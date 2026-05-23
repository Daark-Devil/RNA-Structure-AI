#!/usr/bin/env python3
"""
generate_template_yaml.py -- build a Boltz-2 TEMPLATE YAML for one RNA.

Mirrors generate_boltz_yaml.py (the SS-constraint generator) but produces a
template-based YAML instead of constraints.

Requires our PATCHED Boltz-2 (stock Boltz-2 rejects RNA template chains).

Usage:
  python generate_template_yaml.py \
      --rna 1A1T \
      --pdb /path/to/1A1T.pdb \
      --out /path/to/1A1T_template.yaml \
      --cif_out /path/to/1A1T.cif
"""
import argparse
import gemmi


def pdb_to_cif(pdb_path, cif_path):
    """Convert PDB to CIF with full_sequence populated (needed by Boltz-2)."""
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    for entity in st.entities:
        if entity.polymer_type == gemmi.PolymerType.Unknown:
            continue
        for chain in st[0]:
            polymer = chain.get_polymer()
            if not polymer:
                continue
            if not any(r.subchain in entity.subchains for r in chain):
                continue
            entity.full_sequence = [r.name for r in polymer]
            break
    st.make_mmcif_document().write_file(cif_path)


def seq_from_pdb(pdb_path):
    seen, seq = set(), []
    m = {"A": "A", "G": "G", "C": "C", "U": "U"}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            key = (line[21], line[22:26].strip())
            if key in seen:
                continue
            seen.add(key)
            seq.append(m.get(line[17:20].strip(), "N"))
    return "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--pdb", required=True, help="template structure (PDB)")
    ap.add_argument("--out", required=True, help="output YAML path")
    ap.add_argument("--cif_out", required=True, help="output CIF path (YAML points here)")
    ap.add_argument("--chain_id", default="A")
    args = ap.parse_args()

    seq = seq_from_pdb(args.pdb)
    pdb_to_cif(args.pdb, args.cif_out)

    # gemmi renames the chain during PDB->CIF (A -> Axp). Boltz-2 needs that
    # renamed id as template_id. It is the subchain id with "xp" suffix.
    yaml = (
        "version: 1\n"
        "sequences:\n"
        "  - rna:\n"
        f"      id: {args.chain_id}\n"
        f"      sequence: {seq}\n"
        "templates:\n"
        f"  - cif: {args.cif_out}\n"
        f"    chain_id: [{args.chain_id}]\n"
        f"    template_id: [{args.chain_id}xp]\n"
    )
    with open(args.out, "w") as f:
        f.write(yaml)

    print(f"[info] RNA      : {args.rna}")
    print(f"[info] length   : {len(seq)}")
    print(f"[info] cif      : {args.cif_out}")
    print(f"[done] wrote    : {args.out}")


if __name__ == "__main__":
    main()
