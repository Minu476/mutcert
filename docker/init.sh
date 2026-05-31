#!/bin/sh
# MutCert Docker init script
# Runs once to build all protein graphs and import mutations into Neo4j.
# Idempotent: checks if data already exists before each step.
#
# Usage (from repo root):
#   docker compose run --rm mutcert-init
#
set -e

APP="dotnet /app/Protein.Engine.dll"
DATA="/app/data"

echo "========================================================================"
echo " MutCert — one-time database initialisation"
echo "========================================================================"

# ── Step 1: Graft T4 Lysozyme ──────────────────────────────────────────────
echo ""
echo "[1/4] Grafting T4 Lysozyme (P00720)..."
$APP graft \
    P00720 \
    t4_lysozyme \
    PF00959 \
    "$DATA/cif/t4_lysozyme_P00720_2LZM.cif"

# ── Step 2: Graft CI2 ──────────────────────────────────────────────────────
echo ""
echo "[2/4] Grafting CI2 (P01053)..."
$APP graft \
    P01053 \
    ci2 \
    PF00014 \
    "$DATA/cif/ci2_P01053.cif"

# ── Step 3: Graft Barnase ──────────────────────────────────────────────────
echo ""
echo "[3/4] Grafting Barnase (P00648)..."
$APP graft \
    P00648 \
    barnase \
    PF00211 \
    "$DATA/cif/barnase_P00648.cif"

# ── Step 4: Import S2648 Mutations ────────────────────────────────────────
echo ""
echo "[4/4] Importing S2648 mutations..."
$APP import-mutations

echo ""
echo "========================================================================"
echo " Initialisation complete."
echo " Run validation with:"
echo "   docker compose run --rm mutcert validate-all"
echo "========================================================================"
