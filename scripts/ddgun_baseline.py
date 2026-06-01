#!/usr/bin/env python3
"""
DDGun baseline comparison on the MutCert S2648 val split.

Runs ddgun_seq for each mutation in the val split and computes
Pearson r and Spearman ρ vs experimental ΔΔG.

MSAs are fetched from the ColabFold MMseqs2 API (free, no auth required).
Profile format: DDGun-compatible position-frequency table.

Usage:
    .venv/bin/python scripts/ddgun_baseline.py
"""

import json
import time
import re
import sys
import io
import gzip
import urllib.request
import urllib.parse
import urllib.error
from pathlib import Path
from collections import Counter

import pandas as pd
import numpy as np
from scipy import stats

# ── paths ──────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
SPLIT_JSON  = BASE / "data" / "s2648_split.json"
CIF_DIR     = BASE / "data" / "cif"
PROFILE_DIR = BASE / "data" / "ddgun_profiles"
RESULTS_OUT = BASE / "data" / "ddgun_baseline_results.json"
PROFILE_DIR.mkdir(parents=True, exist_ok=True)

# ── amino-acid helpers ─────────────────────────────────────────────────────────
THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}
AA1 = list("ARNDCQEGHILKMFPSTWYV")

# ── CIF sequence extraction ────────────────────────────────────────────────────
def extract_seq_from_cif(path: str) -> str:
    seq_lines, in_seq = [], False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("_entity_poly.pdbx_seq_one_letter_code_can"):
                in_seq = True
                continue
            if in_seq:
                if line.startswith("_") or line.startswith("#") or line.startswith("loop_"):
                    break
                seq_lines.append(line.replace(";", "").replace(" ", ""))
    return "".join(seq_lines).strip()

# ── ColabFold MSA API ──────────────────────────────────────────────────────────
COLABFOLD_BASE = "https://api.colabfold.com"

def submit_msa(name: str, seq: str) -> str:
    """Submit a sequence to ColabFold MSA API, return ticket ID."""
    fasta = f">{name}\n{seq}"
    body = urllib.parse.urlencode({"q": fasta, "mode": "all"}).encode()
    req = urllib.request.Request(
        f"{COLABFOLD_BASE}/ticket/msa",
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())["id"]


def poll_msa(ticket_id: str, max_wait_s: int = 300) -> bytes:
    """Poll ticket until complete, return raw result bytes."""
    poll_url = f"{COLABFOLD_BASE}/ticket/{ticket_id}"
    deadline = time.time() + max_wait_s
    while time.time() < deadline:
        with urllib.request.urlopen(poll_url, timeout=15) as r:
            info = json.loads(r.read())
        status = info.get("status", "")
        if status == "COMPLETE":
            break
        if status == "ERROR":
            raise RuntimeError(f"ColabFold MSA job failed: {info}")
        print(f"  [{ticket_id[:8]}…] status={status}, waiting…")
        time.sleep(5)
    else:
        raise TimeoutError(f"MSA job {ticket_id} did not complete in {max_wait_s}s")

    dl_url = f"{COLABFOLD_BASE}/result/download/{ticket_id}"
    with urllib.request.urlopen(dl_url, timeout=60) as r:
        return r.read()


def parse_a3m(raw: bytes) -> list[str]:
    """
    Parse A3M from ColabFold result.
    ColabFold returns a gzip-compressed tar archive containing uniref.a3m.
    Returns list of aligned sequences (gaps as '-', inserts stripped).
    """
    import tarfile

    text: str | None = None

    # Case 1: gzip-compressed tar (ColabFold default)
    if raw[:2] == b"\x1f\x8b":
        try:
            tf = tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz")
            for member in tf.getmembers():
                if member.name.endswith(".a3m"):
                    fobj = tf.extractfile(member)
                    if fobj:
                        text = fobj.read().decode(errors="replace")
                        break
        except Exception:
            pass

    # Case 2: zip file
    if text is None and raw[:2] == b"PK":
        import zipfile
        try:
            zf = zipfile.ZipFile(io.BytesIO(raw))
            for name in zf.namelist():
                if name.endswith(".a3m"):
                    text = zf.read(name).decode(errors="replace")
                    break
        except Exception:
            pass

    # Case 3: plain gzip A3M (not tar)
    if text is None and raw[:2] == b"\x1f\x8b":
        try:
            text = gzip.decompress(raw).decode(errors="replace")
        except Exception:
            pass

    # Case 4: plain text
    if text is None:
        text = raw.decode(errors="replace")

    seqs = []
    current: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
            current = []
        elif line:
            # Strip lowercase insert states (A3M convention)
            current.append(re.sub(r"[a-z]", "", line))
    if current:
        seqs.append("".join(current))
    return seqs


def build_profile(seqs: list[str]) -> pd.DataFrame:
    """
    Build a position-frequency profile in DDGun format.
    Columns: POS (1-based index), A-Y (20 AA), SEQ (native), - (gap).
    """
    if not seqs:
        raise ValueError("No sequences in MSA")
    query_seq = seqs[0]
    n_pos = len(query_seq)
    n_seq = len(seqs)

    rows = []
    for i in range(n_pos):
        counts = Counter(s[i] for s in seqs if i < len(s))
        freq = {aa: counts.get(aa, 0) / n_seq for aa in AA1}
        freq["-"] = counts.get("-", 0) / n_seq
        freq["SEQ"] = query_seq[i]
        rows.append(freq)

    df = pd.DataFrame(rows, index=pd.RangeIndex(1, n_pos + 1, name="POS"))
    df.index.name = "POS"
    return df[AA1 + ["-", "SEQ"]]


def get_or_fetch_profile(name: str, seq: str) -> pd.DataFrame:
    """Load cached profile or fetch MSA from ColabFold and build profile."""
    cache_path = PROFILE_DIR / f"{name}.prof"
    if cache_path.exists():
        print(f"  Using cached profile: {cache_path}")
        return pd.read_csv(cache_path, sep=r"\s+", index_col=0)

    print(f"  Submitting MSA job for {name} ({len(seq)} aa)…")
    ticket = submit_msa(name, seq)
    print(f"  Ticket: {ticket}")
    raw = poll_msa(ticket)
    print(f"  Downloaded {len(raw)} bytes")

    seqs = parse_a3m(raw)
    print(f"  Parsed {len(seqs)} sequences from MSA")
    if len(seqs) < 2:
        print(f"  WARNING: only {len(seqs)} seq(s) — using single-seq pseudo-profile")

    profile = build_profile(seqs)
    profile.to_csv(cache_path, sep=" ", float_format="%.6f")
    print(f"  Saved profile: {cache_path}")
    return profile

# ── DDGun runner ───────────────────────────────────────────────────────────────
def run_ddgun_seq(wt1: str, pos1: int, mut1: str, profile_df: pd.DataFrame) -> float:
    """
    Call ddgun.ddgun_seq for a single substitution.
    pos1 is 1-based within the profile.
    """
    from ddgun.aa import Substitution
    from ddgun import ddgun as _ddgun
    from ddgun.profile import Profile

    sub = Substitution(aa_from=wt1, aa_pos=pos1, aa_to=mut1)
    prof = Profile()
    prof.profile = profile_df
    return float(_ddgun.ddgun_seq(sub, prof))

# ── Main ───────────────────────────────────────────────────────────────────────
CIF_PATHS = {
    "t4-lysozyme": CIF_DIR / "t4_lysozyme_P00720_2LZM.cif",
    "ci2":         CIF_DIR / "ci2_P01053.cif",
    "barnase":     CIF_DIR / "barnase_P00648.cif",
}
# Sequence offsets applied by MutCert when reading S2648 positions
# (mature protein position in S2648) → (1-based position in CIF/profile)
# offset: cif_pos = s2648_pos + offset
SEQ_OFFSETS = {
    # DDGun profiles are built directly from AlphaFold CIF sequences (1-indexed).
    # S2648 uses the same UniProt canonical numbering → offset = 0 for all families.
    # (MutCert's Neo4j offsets +1 and +47 corrected for internal graph numbering,
    #  not for raw CIF positions.)
    "t4-lysozyme": 0,
    "ci2":         0,
    "barnase":     0,
}

def main():
    with open(SPLIT_JSON) as f:
        split = json.load(f)

    # Step 1: build profiles
    print("\n=== Step 1: Fetch / build DDGun profiles ===")
    profiles = {}
    seqs_by_fam = {}
    for fam, cif_path in CIF_PATHS.items():
        seq = extract_seq_from_cif(str(cif_path))
        seqs_by_fam[fam] = seq
        print(f"\n{fam}: len={len(seq)}")
        profiles[fam] = get_or_fetch_profile(fam, seq)

    # Step 2: run DDGun on val splits
    print("\n=== Step 2: Run DDGun-seq on val mutations ===")
    all_results = {}

    for fam in ["t4-lysozyme", "ci2", "barnase"]:
        val = split[fam]["Val"]
        offset = SEQ_OFFSETS[fam]
        profile_df = profiles[fam]
        results = []

        print(f"\n--- {fam} (n={len(val)}, offset={offset}) ---")
        for record in val:
            mid = record["MutationId"]
            wt3  = record["WildType3"]
            mut3 = record["MutantType3"]
            # Parse position from MutationId: MUT_<uniprot>_<chain>_<pos>_<mut3>
            parts = mid.split("_")
            s2648_pos = int(parts[3])
            cif_pos = s2648_pos + offset   # 1-based position in CIF/profile
            exp_ddg = record["ExperimentalDdg"]

            wt1  = THREE_TO_ONE.get(wt3)
            mut1 = THREE_TO_ONE.get(mut3)
            if wt1 is None or mut1 is None:
                print(f"  SKIP {mid}: unknown aa {wt3}/{mut3}")
                continue

            # Verify position in range
            n_pos = len(profile_df)
            if cif_pos < 1 or cif_pos > n_pos:
                print(f"  SKIP {mid}: cif_pos={cif_pos} out of range [1,{n_pos}]")
                continue

            # Verify wild-type matches profile SEQ column
            profile_wt = profile_df.at[cif_pos, "SEQ"] if "SEQ" in profile_df.columns else "?"
            if profile_wt != wt1:
                print(f"  WARN {mid}: profile SEQ={profile_wt}, expected {wt1} at pos {cif_pos}")

            try:
                pred_ddg = run_ddgun_seq(wt1, cif_pos, mut1, profile_df)
            except Exception as e:
                print(f"  ERROR {mid}: {e}")
                continue

            results.append({
                "MutationId": mid,
                "CifPos": cif_pos,
                "WT": wt1, "Mut": mut1,
                "DDGun_seq": pred_ddg,
                "Experimental": exp_ddg,
            })

        all_results[fam] = results
        if len(results) >= 2:
            pred  = [r["DDGun_seq"] for r in results]
            exper = [r["Experimental"] for r in results]
            r_p, _ = stats.pearsonr(pred, exper)
            r_s, _ = stats.spearmanr(pred, exper)
            print(f"  n={len(results)}  Pearson r={r_p:.3f}  Spearman ρ={r_s:.3f}")
        else:
            print(f"  Only {len(results)} result(s) — cannot compute correlation")

    # Step 3: overall summary
    print("\n=== Step 3: Summary ===")
    print(f"{'Family':<15} {'n':>4}  {'Pearson r':>10}  {'Spearman ρ':>10}")
    print("-" * 45)
    all_pred, all_exp = [], []
    summary = {}
    for fam in ["t4-lysozyme", "ci2", "barnase"]:
        results = all_results.get(fam, [])
        if len(results) < 2:
            print(f"{fam:<15} {'?':>4}")
            continue
        pred  = [r["DDGun_seq"] for r in results]
        exper = [r["Experimental"] for r in results]
        r_p, _ = stats.pearsonr(pred, exper)
        r_s, _ = stats.spearmanr(pred, exper)
        print(f"{fam:<15} {len(results):>4}  {r_p:>10.3f}  {r_s:>10.3f}")
        all_pred.extend(pred)
        all_exp.extend(exper)
        summary[fam] = {"n": len(results), "pearson_r": round(r_p, 4), "spearman_rho": round(r_s, 4)}

    if len(all_pred) >= 2:
        r_p_all, _ = stats.pearsonr(all_pred, all_exp)
        r_s_all, _ = stats.spearmanr(all_pred, all_exp)
        print(f"{'ALL (pooled)':<15} {len(all_pred):>4}  {r_p_all:>10.3f}  {r_s_all:>10.3f}")
        summary["all_pooled"] = {"n": len(all_pred), "pearson_r": round(r_p_all, 4), "spearman_rho": round(r_s_all, 4)}

    # Save
    output = {"summary": summary, "details": all_results}
    with open(RESULTS_OUT, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to: {RESULTS_OUT}")


if __name__ == "__main__":
    main()
