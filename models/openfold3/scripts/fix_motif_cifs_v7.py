#!/usr/bin/env python3
"""
fix_motif_cifs_v7.py

v6 + key fix: write label_seq_id = auth_seq_id (renumbered values like 907-920),
NOT 1-based position-within-entity (1-14).

Why: OF3's parse_mmcif calls update_author_to_pdb_labels which sets
atom_array.res_id from label_seq_id. The cache's idx_map[:, 1] contains
renumbered values like [907..920]. So label_seq_id MUST match those values
or np.isin will return 0 atoms after filtering.

Usage:
    python fix_motif_cifs_v7.py <of3_structures_dir> <sto_path>
"""
import argparse
import re
import shutil
import sys
from pathlib import Path
from collections import OrderedDict

import biotite.structure.io.pdbx as pdbx


RNA_RESIDUES_1 = {"A": "A", "G": "G", "C": "C", "U": "U"}


def to_one_letter(res3: str) -> str:
    if res3 in RNA_RESIDUES_1:
        return RNA_RESIDUES_1[res3]
    return "N"


def parse_sto_starts(sto_path: Path) -> dict:
    starts = {}
    pattern = re.compile(r"^#=GS\s+([^/]+)/(\d+)-(\d+)")
    for line in sto_path.read_text().splitlines():
        m = pattern.match(line)
        if m:
            full_id = m.group(1)
            entry_id = full_id.split("_")[0]
            if entry_id == "query":
                continue
            starts[entry_id] = int(m.group(2))
    return starts


def fix_one(cif_path: Path, start_residue: int) -> bool:
    print(f"\n--- {cif_path.name} (renumber start={start_residue}) ---")

    # Restore from earliest backup
    for suffix in [".cif.before_fix", ".cif.before_fix_v2", ".cif.before_v3",
                   ".cif.before_v4", ".cif.before_v5", ".cif.before_v6"]:
        backup = cif_path.with_suffix(suffix)
        if backup.is_file():
            shutil.copy2(backup, cif_path)
            print(f"  Restored from {backup.name}")
            break
    else:
        v7_backup = cif_path.with_suffix(".cif.before_v7")
        if not v7_backup.exists():
            shutil.copy2(cif_path, v7_backup)

    cif = pdbx.CIFFile.read(cif_path)
    block_name = list(cif.keys())[0]
    block = cif[block_name]

    if "atom_site" not in block:
        return False

    atom_site = block["atom_site"]
    n_atoms = len(atom_site["group_PDB"].as_array())

    group_PDB    = atom_site["group_PDB"].as_array(dtype=str)
    type_symbol  = atom_site["type_symbol"].as_array(dtype=str)
    label_atom   = atom_site["label_atom_id"].as_array(dtype=str)
    label_comp   = atom_site["label_comp_id"].as_array(dtype=str)
    auth_asym    = atom_site["auth_asym_id"].as_array(dtype=str)
    orig_seq     = atom_site["auth_seq_id"].as_array(dtype=int)
    cartn_x      = atom_site["Cartn_x"].as_array(dtype=float)
    cartn_y      = atom_site["Cartn_y"].as_array(dtype=float)
    cartn_z      = atom_site["Cartn_z"].as_array(dtype=float)

    chains = OrderedDict()
    for asym, seqid, comp in zip(auth_asym, orig_seq, label_comp):
        if asym not in chains:
            chains[asym] = OrderedDict()
        if seqid not in chains[asym]:
            chains[asym][seqid] = comp

    # Renumber residues contiguously starting from start_residue
    seqid_remap = {}
    for chain, residues in chains.items():
        for new_idx, orig_seqid in enumerate(residues.keys()):
            seqid_remap[(chain, orig_seqid)] = start_residue + new_idx

    new_seq_arr = []
    for c, oseq in zip(auth_asym, orig_seq):
        new_seq_arr.append(seqid_remap[(c, oseq)])

    print(f"  New residue IDs: {sorted(set(new_seq_arr))}")

    all_residue_names = sorted(set(label_comp))

    chain_to_entity = {c: i for i, c in enumerate(chains.keys(), start=1)}

    # Rebuild new_chains with new seqids preserved
    new_chains = OrderedDict()
    for c, residues in chains.items():
        new_chains[c] = OrderedDict()
        for new_idx, (oseq, comp) in enumerate(residues.items()):
            new_chains[c][start_residue + new_idx] = comp

    new_label_entity = []
    new_label_asym = []
    for i in range(n_atoms):
        c = auth_asym[i]
        new_label_entity.append(str(chain_to_entity[c]))
        new_label_asym.append(c)

    # Build CIF text
    lines = []
    lines.append(f"data_{block_name}")
    lines.append("#")

    lines.append("loop_")
    lines.append("_entry.id")
    lines.append(block_name)
    lines.append("#")

    # chem_comp
    lines.append("loop_")
    lines.append("_chem_comp.id")
    lines.append("_chem_comp.type")
    lines.append("_chem_comp.mon_nstd_flag")
    lines.append("_chem_comp.name")
    lines.append("_chem_comp.formula")
    for resn in all_residue_names:
        if resn in {"A", "G", "C", "U"}:
            lines.append(f"{resn} 'RNA linking' y . ?")
        else:
            lines.append(f"{resn} 'non-polymer' n . ?")
    lines.append("#")

    # entity
    lines.append("loop_")
    lines.append("_entity.id")
    lines.append("_entity.type")
    lines.append("_entity.pdbx_description")
    for chain, eid in chain_to_entity.items():
        lines.append(f"{eid} polymer ?")
    lines.append("#")

    # entity_poly
    lines.append("loop_")
    lines.append("_entity_poly.entity_id")
    lines.append("_entity_poly.type")
    lines.append("_entity_poly.pdbx_seq_one_letter_code")
    lines.append("_entity_poly.pdbx_seq_one_letter_code_can")
    lines.append("_entity_poly.pdbx_strand_id")
    for chain, residues in new_chains.items():
        eid = chain_to_entity[chain]
        seq3 = list(residues.values())
        seq1 = "".join(to_one_letter(r) for r in seq3)
        lines.append(f"{eid} polyribonucleotide  {seq1}  {seq1}  {chain}")
    lines.append("#")

    # entity_poly_seq — _num field is 1-based position within entity
    lines.append("loop_")
    lines.append("_entity_poly_seq.entity_id")
    lines.append("_entity_poly_seq.num")
    lines.append("_entity_poly_seq.mon_id")
    for chain, residues in new_chains.items():
        eid = chain_to_entity[chain]
        for i, (nseq, comp) in enumerate(residues.items(), start=1):
            lines.append(f"{eid} {i} {comp}")
    lines.append("#")

    # pdbx_poly_seq_scheme — seq_id is 1-based position within entity
    lines.append("loop_")
    lines.append("_pdbx_poly_seq_scheme.asym_id")
    lines.append("_pdbx_poly_seq_scheme.entity_id")
    lines.append("_pdbx_poly_seq_scheme.seq_id")
    lines.append("_pdbx_poly_seq_scheme.mon_id")
    lines.append("_pdbx_poly_seq_scheme.pdb_strand_id")
    lines.append("_pdbx_poly_seq_scheme.pdb_seq_num")
    lines.append("_pdbx_poly_seq_scheme.auth_seq_num")
    for chain, residues in new_chains.items():
        eid = chain_to_entity[chain]
        for i, (nseq, comp) in enumerate(residues.items(), start=1):
            lines.append(f"{chain} {eid} {i} {comp} {chain} {nseq} {nseq}")
    lines.append("#")

    # pdbx_audit_revision_history
    lines.append("loop_")
    lines.append("_pdbx_audit_revision_history.ordinal")
    lines.append("_pdbx_audit_revision_history.data_content_type")
    lines.append("_pdbx_audit_revision_history.major_revision")
    lines.append("_pdbx_audit_revision_history.minor_revision")
    lines.append("_pdbx_audit_revision_history.revision_date")
    lines.append("1 'Structure model' 1 0 2000-01-01")
    lines.append("#")

    # atom_site
    # KEY CHANGE FROM v6: label_seq_id = NEW renumbered residue ID (e.g. 907-920),
    # NOT position-within-entity (1-14). OF3's parse_mmcif uses label_seq_id as
    # res_id, and the cache's idx_map[:, 1] contains the renumbered values.
    lines.append("loop_")
    lines.append("_atom_site.group_PDB")
    lines.append("_atom_site.id")
    lines.append("_atom_site.type_symbol")
    lines.append("_atom_site.label_atom_id")
    lines.append("_atom_site.label_alt_id")
    lines.append("_atom_site.label_comp_id")
    lines.append("_atom_site.label_asym_id")
    lines.append("_atom_site.label_entity_id")
    lines.append("_atom_site.label_seq_id")
    lines.append("_atom_site.pdbx_PDB_ins_code")
    lines.append("_atom_site.Cartn_x")
    lines.append("_atom_site.Cartn_y")
    lines.append("_atom_site.Cartn_z")
    lines.append("_atom_site.occupancy")
    lines.append("_atom_site.B_iso_or_equiv")
    lines.append("_atom_site.pdbx_formal_charge")
    lines.append("_atom_site.auth_seq_id")
    lines.append("_atom_site.auth_comp_id")
    lines.append("_atom_site.auth_asym_id")
    lines.append("_atom_site.auth_atom_id")
    lines.append("_atom_site.pdbx_PDB_model_num")
    for i in range(n_atoms):
        atom_name = label_atom[i]
        if "'" in atom_name or '"' in atom_name:
            atom_name_quoted = f'"{atom_name}"'
        else:
            atom_name_quoted = atom_name
        # KEY: label_seq_id and auth_seq_id BOTH = new_seq_arr[i] (e.g. 907-920)
        lines.append(
            f"{group_PDB[i]} {i+1} {type_symbol[i]} {atom_name_quoted} . {label_comp[i]} "
            f"{new_label_asym[i]} {new_label_entity[i]} {new_seq_arr[i]} ? "
            f"{cartn_x[i]:.3f} {cartn_y[i]:.3f} {cartn_z[i]:.3f} 1.00 0.00 ? "
            f"{new_seq_arr[i]} {label_comp[i]} {auth_asym[i]} {atom_name_quoted} 1"
        )
    lines.append("#")

    new_text = "\n".join(lines) + "\n"
    cif_path.write_text(new_text)
    print(f"  Rewrote {n_atoms} atoms — label_seq_id = auth_seq_id = {start_residue}..{start_residue + len(new_chains[list(new_chains.keys())[0]]) - 1}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rna_dir")
    ap.add_argument("sto_path")
    args = ap.parse_args()

    base = Path(args.rna_dir)
    sto = Path(args.sto_path)
    if not base.is_dir() or not sto.is_file():
        sys.exit(1)

    starts = parse_sto_starts(sto)
    print(f"Parsed {len(starts)} starts: {starts}")

    cifs = sorted(p for p in base.glob("*.cif") if p.is_file())
    fixed = 0
    for c in cifs:
        entry_id = c.stem
        if entry_id not in starts:
            print(f"  ⚠ no start for {entry_id}, skipping")
            continue
        if fix_one(c, starts[entry_id]):
            fixed += 1
    print(f"\n[OK] Fixed {fixed}/{len(cifs)}")


if __name__ == "__main__":
    main()
