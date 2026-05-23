"""
trim_template_pdbs.py v4 — Extract motif-relevant residues from template PDB files.

Changes from v3:
  - FIX: Helix templates now have CORRECT residue names matching the actual RNA
    sequence from the MTF, not the universal helix's default A/U sequence.
    This prevents RiNALMo from getting wrong sequence embeddings.

Usage:
  python trim_template_pdbs.py \
      --templates_tsv results/R1116_manual/R1116_vfold3D.mtf.templates.tsv \
      --templates_dir results/R1116_manual/templates/ \
      --mtf results/R1116_manual/R1116_vfold3D_clean.mtf \
      --helix_pdb /path/to/VfoldPipeline_standalone/Vfold3DLA/Data/data/vfold3D/helix.pdb
"""

import os
import re
import argparse
from collections import defaultdict


# ======================== HELIX EXTRACTION ========================

def parse_helix_motifs_from_mtf(mtf_path: str) -> list[dict]:
    """
    Parse HELIX lines from the clean MTF file.
    
    Example line:
      HELIX 13 chain  A 28 A 40 A 112 A 124 sequence GCGGCUAAAACAG UGUUUUAGCCGC
    """
    helix_motifs = []
    helix_idx = 0
    
    with open(mtf_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if not line.startswith('HELIX'):
                continue
            
            parts = line.split()
            num_bp = int(parts[1])
            
            seq_idx = parts.index('sequence')
            seqs = parts[seq_idx + 1:]
            strand1_seq = seqs[0]
            strand2_seq = seqs[1] if len(seqs) > 1 else ''
            
            helix_motifs.append({
                'helix_idx': helix_idx,
                'num_bp': num_bp,
                'strand1_len': len(strand1_seq),
                'strand2_len': len(strand2_seq),
                'strand1_seq': strand1_seq,
                'strand2_seq': strand2_seq,
                'motif_line': line,
            })
            helix_idx += 1
    
    return helix_motifs


def extract_helix_from_universal(helix_pdb_path: str, num_bp: int,
                                  strand1_seq: str = '', strand2_seq: str = '') -> list[str]:
    """
    Extract N base pairs from the universal helix.pdb, CENTER coordinates,
    and REWRITE residue names to match the actual RNA sequence.
    
    The universal helix has 258 residues in chain A (129 base pairs):
      - Strand 1: residues 1 to 129
      - Strand 2: residues 130 to 258 (pairs with strand 1 in reverse)
    
    For a helix of N base pairs:
      - Take residues 1 to N (strand 1) → rewrite as strand1_seq
      - Take residues (258 - N + 1) to 258 (strand 2) → rewrite as strand2_seq
    """
    strand1_residues = set(range(1, num_bp + 1))
    strand2_start = 258 - num_bp + 1
    strand2_residues = set(range(strand2_start, 259))
    keep_residues = strand1_residues | strand2_residues
    
    # Build residue number → new nucleotide name mapping
    resnum_to_newname = {}
    for i, resnum in enumerate(range(1, num_bp + 1)):
        if i < len(strand1_seq):
            resnum_to_newname[resnum] = strand1_seq[i]  # e.g., 'G', 'C', 'A', 'U'
    for i, resnum in enumerate(range(strand2_start, 259)):
        if i < len(strand2_seq):
            resnum_to_newname[resnum] = strand2_seq[i]
    
    kept_lines = []
    coords = []
    with open(helix_pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                if resnum in keep_residues:
                    # Rewrite residue name if we have the correct one
                    if resnum in resnum_to_newname:
                        new_name = resnum_to_newname[resnum]
                        # PDB format: residue name is at columns 17-19 (0-indexed)
                        # For RNA: single letter like 'A', 'G', 'C', 'U' padded with spaces
                        line = line[:17] + f"  {new_name}" + line[20:]
                    kept_lines.append(line)
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
    
    # Center coordinates around origin
    if coords:
        cx = sum(c[0] for c in coords) / len(coords)
        cy = sum(c[1] for c in coords) / len(coords)
        cz = sum(c[2] for c in coords) / len(coords)
        
        centered_lines = []
        for line in kept_lines:
            x = float(line[30:38]) - cx
            y = float(line[38:46]) - cy
            z = float(line[46:54]) - cz
            new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            centered_lines.append(new_line)
        kept_lines = centered_lines
    
    kept_lines.append("END\n")
    return kept_lines


def generate_helix_templates(mtf_path: str, helix_pdb_path: str, 
                              templates_dir: str, dry_run: bool = False) -> int:
    """
    Generate template PDB files for all HELIX motifs in the MTF.
    Creates files like: HELIX_0_universal.pdb, HELIX_1_universal.pdb, etc.
    """
    helix_motifs = parse_helix_motifs_from_mtf(mtf_path)
    
    if not helix_motifs:
        return 0
    
    print(f"\n  --- Generating HELIX templates from universal helix.pdb ---")
    created = 0
    
    for hm in helix_motifs:
        idx = hm['helix_idx']
        num_bp = hm['num_bp']
        out_filename = f"HELIX_{idx}_universal.pdb"
        out_path = os.path.join(templates_dir, out_filename)
        
        kept_lines = extract_helix_from_universal(
            helix_pdb_path, num_bp,
            strand1_seq=hm['strand1_seq'],
            strand2_seq=hm['strand2_seq']
        )
        kept_atoms = sum(1 for l in kept_lines if l.startswith("ATOM"))
        
        expected_residues = num_bp * 2
        actual_residues = len(set(
            int(l[22:26].strip()) for l in kept_lines if l.startswith("ATOM")
        ))
        
        # Verify sequence was rewritten correctly
        actual_seq = ''
        seen_res = {}
        for l in kept_lines:
            if l.startswith("ATOM"):
                rnum = int(l[22:26].strip())
                if rnum not in seen_res:
                    seen_res[rnum] = l[17:20].strip()
                    actual_seq += l[17:20].strip()
        
        status = "OK" if actual_residues == expected_residues else f"WARN({actual_residues}!={expected_residues})"
        print(f"  {out_filename:45s} | {num_bp:>3d} bp → {kept_atoms:>5d} atoms, "
              f"{actual_residues} residues | seq: {actual_seq} | {status}")
        
        if not dry_run:
            with open(out_path, 'w') as f:
                f.writelines(kept_lines)
            created += 1
    
    return created


# ======================== TEMPLATE TRIMMING ========================

def parse_chain_field(chain_str: str) -> list[tuple[str, int, int]]:
    """Parse the chain field from templates.tsv into (chain_id, start, end) ranges."""
    tokens = chain_str.strip().split()
    if len(tokens) % 2 != 0:
        print(f"  [WARN] Odd number of tokens in chain field: {chain_str}")
        return []
    
    pairs = []
    for i in range(0, len(tokens), 2):
        chain_id = tokens[i]
        try:
            resnum = int(tokens[i + 1])
        except ValueError:
            print(f"  [WARN] Cannot parse: {tokens[i]} {tokens[i+1]}")
            continue
        pairs.append((chain_id, resnum))
    
    ranges = []
    for i in range(0, len(pairs), 2):
        if i + 1 < len(pairs):
            chain1, start = pairs[i]
            chain2, end = pairs[i + 1]
            if chain1 != chain2:
                print(f"  [WARN] Chain mismatch in range: {chain1} vs {chain2}")
            ranges.append((chain1, start, end))
    
    return ranges


def extract_residues_from_pdb(pdb_path: str, ranges: list[tuple[str, int, int]]) -> list[str]:
    """Extract ATOM lines from a PDB file matching the given chain+residue ranges."""
    keep_residues = set()
    for chain_id, start, end in ranges:
        for resnum in range(start, end + 1):
            keep_residues.add((chain_id, resnum))
    
    kept_lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21].strip()
                try:
                    resnum = int(line[22:26].strip())
                except ValueError:
                    continue
                if (chain, resnum) in keep_residues:
                    kept_lines.append(line)
    
    kept_lines.append("END\n")
    return kept_lines


def parse_templates_tsv(tsv_path: str) -> dict:
    """Parse templates.tsv. Keep ONLY the first entry per (motif_idx, pdb_id)."""
    pdb_ranges = {}
    
    with open(tsv_path) as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            motif_idx = parts[0]
            pdb_id = parts[2]
            chain_field = parts[4]
            
            if pdb_id == 'NA':
                continue
            
            key = (motif_idx, pdb_id)
            if key in pdb_ranges:
                continue
            
            ranges = parse_chain_field(chain_field)
            if not ranges:
                continue
            
            pdb_ranges[key] = ranges
    
    return pdb_ranges


def find_template_file(templates_dir: str, motif_idx: str, pdb_id: str) -> str | None:
    """Find the template PDB file matching a motif_idx and pdb_id."""
    for f in os.listdir(templates_dir):
        if not f.endswith('.pdb'):
            continue
        if f.startswith(f"{motif_idx}_") and f.endswith(f"_{pdb_id}.pdb"):
            return os.path.join(templates_dir, f)
    return None


# ======================== MAIN ========================

def main():
    parser = argparse.ArgumentParser(
        description="Trim template PDB files and generate helix templates")
    parser.add_argument('--templates_tsv', required=True,
                        help='Path to templates.tsv (from vfold3D_template_finder)')
    parser.add_argument('--templates_dir', required=True,
                        help='Path to templates/ folder with PDB files')
    parser.add_argument('--mtf', default=None,
                        help='Path to clean MTF file (needed for helix extraction)')
    parser.add_argument('--helix_pdb', default=None,
                        help='Path to universal helix.pdb from Vfold data')
    parser.add_argument('--dry_run', action='store_true',
                        help='Show what would be done without modifying files')
    args = parser.parse_args()
    
    # --- Part 1: Trim non-HELIX template PDBs ---
    print(f"Reading templates TSV: {args.templates_tsv}")
    pdb_ranges = parse_templates_tsv(args.templates_tsv)
    print(f"Found {len(pdb_ranges)} template entries (first chain only per motif)\n")
    
    trimmed = 0
    for (motif_idx, pdb_id), ranges in sorted(pdb_ranges.items()):
        pdb_file = find_template_file(args.templates_dir, motif_idx, pdb_id)
        if pdb_file is None:
            continue
        
        filename = os.path.basename(pdb_file)
        orig_lines = sum(1 for line in open(pdb_file) if line.startswith("ATOM"))
        
        kept_lines = extract_residues_from_pdb(pdb_file, ranges)
        kept_atoms = sum(1 for l in kept_lines if l.startswith("ATOM"))
        
        ranges_str = ", ".join(f"chain {c} {s}-{e}" for c, s, e in ranges)
        print(f"  {filename:45s} | {orig_lines:>6d} → {kept_atoms:>5d} atoms | {ranges_str}")
        
        if not args.dry_run:
            with open(pdb_file, 'w') as f:
                f.writelines(kept_lines)
            trimmed += 1
    
    print(f"\nTrimmed {trimmed} PDB files in {args.templates_dir}")
    
    # --- Part 2: Generate HELIX templates from universal helix.pdb ---
    if args.mtf and args.helix_pdb:
        if not os.path.exists(args.helix_pdb):
            print(f"\n  [ERROR] helix.pdb not found: {args.helix_pdb}")
        elif not os.path.exists(args.mtf):
            print(f"\n  [ERROR] MTF not found: {args.mtf}")
        else:
            created = generate_helix_templates(
                args.mtf, args.helix_pdb, args.templates_dir, args.dry_run
            )
            print(f"\n  Created {created} HELIX template PDBs from universal helix.pdb")
    elif args.mtf or args.helix_pdb:
        print("\n  [NOTE] Both --mtf and --helix_pdb are needed for helix template generation.")
    else:
        print("\n  [NOTE] No --mtf/--helix_pdb provided. Skipping helix template generation.")


if __name__ == "__main__":
    main()
