using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Protein.Engine;

public class Residue
{
    public int SeqPos { get; set; }
    public string ResidueName { get; set; } = "";
    public string ChainId { get; set; } = "";
    public double PLDDT { get; set; }
    
    public Vector3? N { get; set; }
    public Vector3? CA { get; set; }
    public Vector3? C { get; set; }
    public Vector3? O { get; set; }
    public Vector3? CB { get; set; }

    public List<AtomSite> AllAtoms { get; } = new();

    public double? Phi { get; set; }
    public double? Psi { get; set; }

    public Vector3? SideChainCentroid { get; set; }

    /// <summary>
    /// Computes the centroid of the side chain heavy atoms (excluding N, CA, C, O backbone atoms).
    /// If no side chain heavy atoms are found (e.g. GLY), falls back to CA coordinate.
    /// </summary>
    public void ComputeSideChainCentroid()
    {
        var heavyAtoms = AllAtoms.Where(a => 
            !a.AtomName.Trim().Equals("N", StringComparison.OrdinalIgnoreCase) &&
            !a.AtomName.Trim().Equals("CA", StringComparison.OrdinalIgnoreCase) &&
            !a.AtomName.Trim().Equals("C", StringComparison.OrdinalIgnoreCase) &&
            !a.AtomName.Trim().Equals("O", StringComparison.OrdinalIgnoreCase)
        ).ToList();

        if (heavyAtoms.Count > 0)
        {
            Vector3 sum = Vector3.Zero;
            foreach (var atom in heavyAtoms)
            {
                sum += new Vector3((float)atom.X, (float)atom.Y, (float)atom.Z);
            }
            SideChainCentroid = sum / heavyAtoms.Count;
        }
        else
        {
            if (CA.HasValue)
            {
                SideChainCentroid = CA.Value;
            }
        }
    }
}
