#!/usr/bin/env python3
"""
Visualize RNA structure comparison: reference vs baseline (no MSA) vs with-MSA.
Outputs a 4-panel PNG: reference alone, baseline overlay on ref, MSA overlay on ref,
and all three overlaid.

Uses py3Dmol if available, else falls back to matplotlib + Bio.PDB.
"""
import argparse
import sys
from pathlib import Path

# Use Bio.PDB and matplotlib for portability
import numpy as np
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa
except ImportError:
    print("Need matplotlib: pip install matplotlib --break-system-packages")
    sys.exit(1)

try:
    from Bio.PDB import PDBParser, MMCIFParser, Superimposer
except ImportError:
    print("Need biopython: pip install biopython --break-system-packages")
    sys.exit(1)


def get_c1prime_coords(structure):
    """Extract C1' atom coordinates from RNA structure."""
    coords = []
    resids = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "C1'" in residue:
                    coords.append(residue["C1'"].coord)
                    resids.append(residue.get_id()[1])
        break  # only first model
    return np.array(coords), resids


def load_structure(path):
    path = Path(path)
    if path.suffix.lower() == '.cif':
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def kabsch_align(P, Q):
    """Align P onto Q using Kabsch algorithm. Returns aligned P and RMSD."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    Pa = (Pc @ R.T) + Q.mean(axis=0)
    rmsd = np.sqrt(((Pa - Q) ** 2).sum() / len(Pa))
    return Pa, rmsd


def plot_3d(ax, coords_list, labels, colors, title):
    for coords, label, color in zip(coords_list, labels, colors):
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
                '-', color=color, label=label, linewidth=2, alpha=0.8)
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                   c=color, s=15, alpha=0.6)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rna', required=True, help='RNA name for title (e.g. 1NBS)')
    p.add_argument('--ref', required=True, help='Reference PDB file')
    p.add_argument('--baseline', required=True, help='Baseline (no MSA) CIF')
    p.add_argument('--msa', required=True, help='With-MSA CIF')
    p.add_argument('--out', required=True, help='Output PNG path')
    args = p.parse_args()

    print(f"Loading reference: {args.ref}")
    ref_struct = load_structure(args.ref)
    ref_coords, ref_ids = get_c1prime_coords(ref_struct)

    print(f"Loading baseline: {args.baseline}")
    base_struct = load_structure(args.baseline)
    base_coords, base_ids = get_c1prime_coords(base_struct)

    print(f"Loading MSA: {args.msa}")
    msa_struct = load_structure(args.msa)
    msa_coords, msa_ids = get_c1prime_coords(msa_struct)

    # Match lengths to common residues
    common_n = min(len(ref_coords), len(base_coords), len(msa_coords))
    ref_coords = ref_coords[:common_n]
    base_coords = base_coords[:common_n]
    msa_coords = msa_coords[:common_n]

    # Align baseline to reference, MSA to reference
    base_aligned, base_rmsd = kabsch_align(base_coords, ref_coords)
    msa_aligned, msa_rmsd = kabsch_align(msa_coords, ref_coords)

    print(f"  Baseline RMSD: {base_rmsd:.2f} Å")
    print(f"  MSA RMSD:      {msa_rmsd:.2f} Å")
    print(f"  Improvement:   {base_rmsd - msa_rmsd:+.2f} Å")

    # Plot 4-panel figure
    fig = plt.figure(figsize=(16, 12))
    
    # Panel 1: Reference alone
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    plot_3d(ax1,
            [ref_coords],
            [f'Reference (native)'],
            ['black'],
            f'{args.rna}: Reference (n={common_n} residues)')

    # Panel 2: Reference vs baseline
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    plot_3d(ax2,
            [ref_coords, base_aligned],
            ['Reference', f'Baseline (no MSA)'],
            ['black', 'red'],
            f'{args.rna}: Baseline vs Reference (RMSD={base_rmsd:.2f} Å)')

    # Panel 3: Reference vs MSA
    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    plot_3d(ax3,
            [ref_coords, msa_aligned],
            ['Reference', f'With Synthetic MSA'],
            ['black', 'green'],
            f'{args.rna}: With MSA vs Reference (RMSD={msa_rmsd:.2f} Å)')

    # Panel 4: All three
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    plot_3d(ax4,
            [ref_coords, base_aligned, msa_aligned],
            ['Reference', 'Baseline', 'With MSA'],
            ['black', 'red', 'green'],
            f'{args.rna}: All three (improvement = {base_rmsd - msa_rmsd:+.2f} Å)')

    plt.suptitle(f'{args.rna}: Synthetic MSA Effect on RNA Structure Prediction',
                 fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(args.out, dpi=120, bbox_inches='tight')
    print(f"Saved: {args.out}")


if __name__ == '__main__':
    main()
