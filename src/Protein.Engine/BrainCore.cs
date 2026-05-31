using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Protein.Engine;

public static class BrainCore
{
    public static bool ValidateStericClash(List<Residue> residues)
    {
        int n = residues.Count;
        for (int i = 0; i < n; i++)
        {
            var r1 = residues[i];
            if (!r1.CA.HasValue) continue;

            for (int j = i + 1; j < n; j++)
            {
                var r2 = residues[j];
                // Non-adjacent only (sequence distance > 1)
                if (Math.Abs(r1.SeqPos - r2.SeqPos) <= 1 && r1.ChainId == r2.ChainId)
                    continue;

                if (!r2.CA.HasValue) continue;

                float dist = Vector3.Distance(r1.CA.Value, r2.CA.Value);
                if (dist < 3.0f)
                {
                    Console.WriteLine($"Steric Clash Detected: {r1.ResidueName}{r1.SeqPos} & {r2.ResidueName}{r2.SeqPos} (Dist: {dist:F2} Å)");
                    return false; // Clash found
                }
            }
        }
        return true;
    }

    public static bool ValidateBackbonePlanarity(List<Residue> residues)
    {
        var chains = residues.GroupBy(r => r.ChainId).ToList();

        foreach (var chain in chains)
        {
            var chainResidues = chain.OrderBy(r => r.SeqPos).ToList();
            int n = chainResidues.Count;

            for (int i = 0; i < n - 1; i++)
            {
                var r1 = chainResidues[i];
                var r2 = chainResidues[i + 1];

                if (r1.C.HasValue && r2.N.HasValue)
                {
                    float dist = Vector3.Distance(r1.C.Value, r2.N.Value);
                    if (dist < 1.17f || dist > 1.47f) // 1.32 +/- 0.15
                    {
                        Console.WriteLine($"Backbone Planarity Violation: {r1.ResidueName}{r1.SeqPos}-C to {r2.ResidueName}{r2.SeqPos}-N (Dist: {dist:F2} Å)");
                        return false;
                    }
                }
            }
        }
        return true;
    }

    public static bool ValidateLocalChargeNeutrality(List<Residue> residues, Vector3 center, float radius = 10.0f)
    {
        double totalCharge = 0.0;

        foreach (var r in residues)
        {
            if (r.SideChainCentroid.HasValue)
            {
                float dist = Vector3.Distance(center, r.SideChainCentroid.Value);
                if (dist <= radius)
                {
                    totalCharge += MiyazawaJernigan.GetCharge(r.ResidueName);
                }
            }
            else if (r.CA.HasValue)
            {
                 float dist = Vector3.Distance(center, r.CA.Value);
                 if (dist <= radius)
                 {
                     totalCharge += MiyazawaJernigan.GetCharge(r.ResidueName);
                 }
            }
        }

        if (Math.Abs(totalCharge) > 3.0)
        {
            Console.WriteLine($"Charge Neutrality Warning: High local charge within {radius}Å shell ({totalCharge:F1})");
            if (Math.Abs(totalCharge) > 5.0) return false;
        }

        return true;
    }
}
