using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace Protein.Engine;

public record AtomSite(
    string GroupPdb,
    int AtomId,
    string TypeSymbol,
    string AtomName,
    string ResidueName,
    string ChainId,
    int SeqPos,
    double X,
    double Y,
    double Z,
    double Occupancy,
    double TempFactor
);

public static class MmcifParser
{
    /// <summary>
    /// Memory-efficient streaming mmCIF parser that reads a .cif file line-by-line,
    /// dynamically maps the _atom_site column indices, and yields AtomSite records on the fly.
    /// </summary>
    public static IEnumerable<AtomSite> StreamAtoms(string filePath)
    {
        if (!File.Exists(filePath))
        {
            throw new FileNotFoundException($"The specified .mmCIF file was not found: '{filePath}'");
        }

        using var reader = new StreamReader(filePath);
        string? line;
        bool insideAtomSiteLoop = false;
        var colIndex = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        int headerIndex = 0;

        while ((line = reader.ReadLine()) != null)
        {
            string trimmed = line.Trim();
            if (string.IsNullOrEmpty(trimmed)) continue;

            // Handle loop boundary or comments
            if (trimmed.StartsWith('#'))
            {
                // Comment lines can demarcate the end of an atom_site loop or act as separators
                continue;
            }

            if (trimmed.Equals("loop_", StringComparison.OrdinalIgnoreCase))
            {
                insideAtomSiteLoop = false;
                colIndex.Clear();
                headerIndex = 0;
                continue;
            }

            if (trimmed.StartsWith('_'))
            {
                // If this is a loop header tag
                if (trimmed.StartsWith("_atom_site.", StringComparison.OrdinalIgnoreCase))
                {
                    insideAtomSiteLoop = true;
                    colIndex[trimmed] = headerIndex++;
                }
                else
                {
                    // Any other loop header tag means we are not in the atom_site loop
                    insideAtomSiteLoop = false;
                }
                continue;
            }

            // If we get here, this is a data row. Let's see if we are in the atom_site loop.
            if (insideAtomSiteLoop && colIndex.Count > 0)
            {
                if (trimmed.StartsWith("ATOM", StringComparison.OrdinalIgnoreCase) || 
                    trimmed.StartsWith("HETATM", StringComparison.OrdinalIgnoreCase))
                {
                    List<string> tokens = TokenizeLine(trimmed);
                    if (tokens.Count >= colIndex.Count)
                    {
                        AtomSite? atom = ParseAtomSite(tokens, colIndex);
                        if (atom != null)
                        {
                            yield return atom;
                        }
                    }
                }
            }
        }
    }

    /// <summary>
    /// Standard quoted-string character-by-character tokenization to handle single/double quotes and whitespace.
    /// </summary>
    public static List<string> TokenizeLine(string line)
    {
        var tokens = new List<string>();
        int i = 0;
        int len = line.Length;

        while (i < len)
        {
            // Skip leading whitespace
            while (i < len && char.IsWhiteSpace(line[i]))
            {
                i++;
            }
            if (i >= len) break;

            char c = line[i];
            if (c == '"')
            {
                i++;
                int start = i;
                while (i < len && line[i] != '"')
                {
                    i++;
                }
                tokens.Add(line.Substring(start, i - start));
                if (i < len) i++; // skip closing quote
            }
            else if (c == '\'')
            {
                i++;
                int start = i;
                while (i < len && line[i] != '\'')
                {
                    i++;
                }
                tokens.Add(line.Substring(start, i - start));
                if (i < len) i++; // skip closing quote
            }
            else
            {
                int start = i;
                while (i < len && !char.IsWhiteSpace(line[i]))
                {
                    i++;
                }
                tokens.Add(line.Substring(start, i - start));
            }
        }

        return tokens;
    }

    private static AtomSite? ParseAtomSite(List<string> tokens, Dictionary<string, int> colIndex)
    {
        try
        {
            string? groupPdb = GetVal(tokens, colIndex, "_atom_site.group_PDB");
            string? atomIdStr = GetVal(tokens, colIndex, "_atom_site.id");
            string? typeSymbol = GetVal(tokens, colIndex, "_atom_site.type_symbol");
            string? atomName = GetVal(tokens, colIndex, "_atom_site.label_atom_id", "_atom_site.auth_atom_id");
            string? residueName = GetVal(tokens, colIndex, "_atom_site.label_comp_id", "_atom_site.auth_comp_id");
            string? chainId = GetVal(tokens, colIndex, "_atom_site.label_asym_id", "_atom_site.auth_asym_id");
            string? seqPosStr = GetVal(tokens, colIndex, "_atom_site.label_seq_id", "_atom_site.auth_seq_id");
            string? xStr = GetVal(tokens, colIndex, "_atom_site.Cartn_x");
            string? yStr = GetVal(tokens, colIndex, "_atom_site.Cartn_y");
            string? zStr = GetVal(tokens, colIndex, "_atom_site.Cartn_z");
            string? occupancyStr = GetVal(tokens, colIndex, "_atom_site.occupancy");
            string? tempFactorStr = GetVal(tokens, colIndex, "_atom_site.B_iso_or_equiv");

            if (groupPdb == null || atomIdStr == null || typeSymbol == null || atomName == null ||
                residueName == null || chainId == null || seqPosStr == null || xStr == null || yStr == null || zStr == null)
            {
                return null;
            }

            if (!int.TryParse(atomIdStr, NumberStyles.Integer, CultureInfo.InvariantCulture, out int atomId)) return null;
            if (!int.TryParse(seqPosStr, NumberStyles.Integer, CultureInfo.InvariantCulture, out int seqPos)) return null;
            if (!double.TryParse(xStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double x)) return null;
            if (!double.TryParse(yStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double y)) return null;
            if (!double.TryParse(zStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double z)) return null;

            double occupancy = 1.0;
            if (occupancyStr != null)
            {
                double.TryParse(occupancyStr, NumberStyles.Float, CultureInfo.InvariantCulture, out occupancy);
            }

            double tempFactor = 0.0;
            if (tempFactorStr != null)
            {
                double.TryParse(tempFactorStr, NumberStyles.Float, CultureInfo.InvariantCulture, out tempFactor);
            }

            return new AtomSite(
                groupPdb,
                atomId,
                typeSymbol,
                atomName,
                residueName,
                chainId,
                seqPos,
                x,
                y,
                z,
                occupancy,
                tempFactor
            );
        }
        catch
        {
            return null;
        }
    }

    private static string? GetVal(List<string> tokens, Dictionary<string, int> colIndex, string primaryKey, string? fallbackKey = null)
    {
        if (colIndex.TryGetValue(primaryKey, out int idx) && idx < tokens.Count)
        {
            return tokens[idx];
        }
        if (fallbackKey != null && colIndex.TryGetValue(fallbackKey, out idx) && idx < tokens.Count)
        {
            return tokens[idx];
        }
        return null;
    }
}
