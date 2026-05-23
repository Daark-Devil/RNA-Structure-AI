#!/usr/bin/env python3
import re
import json
import shutil
from pathlib import Path
import gemmi

RNA_MAP = {"A":"A","ADE":"A","U":"U","URA":"U","G":"G","GUA":"G","C":"C","CYT":"C"}

def parse_motif_positions_and_seqs(motif_line: str):
    m = re.search(r'chain\s+(.*?)\s+sequence\s+(.+)$', motif_line)
    if not m:
        return None, None
    chain_part = m.group(1).strip()
    seq_part = m.group(2).strip()

    toks = chain_part.split()
    pos = []
    i = 0
    while i < len(toks):
        if i + 1 < len(toks):
            try:
                int(toks[i+1])
                pos.append(int(toks[i+1]))
                i += 2
                continue
            except:
                pass
        i += 1

    seqs = seq_part.split()
    if len(pos) != 2 * len(seqs):
        return None, None

    pieces = []
    for j, s in enumerate(seqs):
        start = pos[2*j]
        end = pos[2*j+1]
        pieces.append((start, end, s))
    return pieces, seqs

def make_row(query_len: int, pieces):
    arr = ["-"] * query_len
    for start, end, seq in pieces:
        start0 = start - 1
        arr[start0:start0+len(seq)] = list(seq)
    return "".join(arr)

def extract_rna_chains(cif_path: Path):
    st = gemmi.read_structure(str(cif_path))
    chains = []
    for model in st:
        for chain in model:
            seq = []
            nums = []
            for res in chain:
                name = res.name.strip().upper()
                if name in RNA_MAP:
                    seq.append(RNA_MAP[name])
                    nums.append(res.seqid.num)
            if seq:
                chains.append((chain.name, "".join(seq), nums[0], nums[-1]))
        break
    return chains

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-name", required=True)
    ap.add_argument("--seq-file", required=True)
    ap.add_argument("--templates-tsv", required=True)
    ap.add_argument("--templates-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    rna = args.rna_name
    seq = Path(args.seq_file).read_text().splitlines()[0].strip()
    qlen = len(seq)

    templates_tsv = Path(args.templates_tsv)
    templates_dir = Path(args.templates_dir)
    out = Path(args.out_dir)

    of3_templates = out / "of3_templates"
    of3_templates_cif = out / "of3_templates_cif"
    of3_structures = out / "of3_structures"
    template_data = out / "template_data"

    of3_templates.mkdir(parents=True, exist_ok=True)
    of3_templates_cif.mkdir(parents=True, exist_ok=True)
    of3_structures.mkdir(parents=True, exist_ok=True)
    (template_data / "template_cache").mkdir(parents=True, exist_ok=True)
    (template_data / "logs").mkdir(parents=True, exist_ok=True)

    kept = []
    seen = set()

    with open(templates_tsv) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            motif_idx = parts[0]
            motif = parts[1]
            pdb = parts[2]

            if motif_idx in seen:
                continue
            seen.add(motif_idx)

            # skip helices and tails for now
            if motif.startswith("HELIX") or motif.startswith("TAIL"):
                continue

            motif_name = motif.split(" chain ")[0].replace(" ", "_")
            cand = list(templates_dir.glob(f"{motif_idx}_{motif_name}_{pdb}.pdb"))
            if not cand:
                cand = list(templates_dir.glob(f"{motif_idx}_*_{pdb}.pdb"))
            if not cand:
                continue

            src_pdb = cand[0]
            dst_pdb = of3_templates / src_pdb.name
            shutil.copy2(src_pdb, dst_pdb)

            st = gemmi.read_structure(str(dst_pdb))
            cif_path = of3_templates_cif / f"{dst_pdb.stem}.cif"
            st.make_mmcif_document().write_file(str(cif_path))

            chains = extract_rna_chains(cif_path)
            if len(chains) != 1:
                continue

            chain_name, chain_seq, start_num, end_num = chains[0]

            pieces, seqs = parse_motif_positions_and_seqs(motif)
            if pieces is None:
                continue

            target_concat = "".join(s for _, _, s in pieces)
            if chain_seq != target_concat:
                continue

            synth_id = f"m{motif_idx}"
            row_name = f"{synth_id}_{chain_name}/{start_num}-{end_num}"
            row = make_row(qlen, pieces)

            d = of3_structures / synth_id
            d.mkdir(exist_ok=True)
            shutil.copy2(cif_path, d / f"{synth_id}_{chain_name}.cif")

            kept.append({
                "motif_idx": motif_idx,
                "motif": motif,
                "pdb": pdb,
                "chain": chain_name,
                "start": start_num,
                "end": end_num,
                "row_name": row_name,
                "entry_chain_id": f"{synth_id}_{chain_name}",
                "row": row
            })

    if not kept:
        raise SystemExit("No usable motif templates found")

    sto = out / f"{rna}_templates.sto"
    with open(sto, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n")
        for k in kept:
            fh.write(f"#=GS {k['row_name']} mol:rna\n")
        fh.write(f"{'query':<22}{seq}\n")
        for k in kept:
            fh.write(f"{k['row_name']:<22}{k['row']}\n")
        fh.write("//\n")

    qjson = out / f"{rna}_query_with_templates.json"
    with open(qjson, "w") as fh:
        json.dump({
            "queries": {
                f"{rna}_motif_test": {
                    "use_msas": False,
                    "chains": [{
                        "molecule_type": "rna",
                        "chain_ids": "A",
                        "sequence": seq,
                        "template_alignment_file_path": str(sto),
                        "template_entry_chain_ids": [k["entry_chain_id"] for k in kept]
                    }]
                }
            }
        }, fh, indent=2)

    runner = out / f"runner_templates_{rna}.yml"
    runner.write_text(
f"""template_preprocessor_settings:
  mode: predict
  moltypes:
    - 1
  max_sequences_parse: 200
  max_templates: 50
  fetch_missing_structures: false
  create_precache: false
  preparse_structures: false
  create_logs: true
  n_processes: 1
  chunksize: 1
  structure_directory: {of3_structures}
  structure_file_format: cif
  output_directory: {template_data}
  precache_directory: null
  structure_array_directory: null
  cache_directory: {template_data / 'template_cache'}
  log_directory: {template_data / 'logs'}
""")

    manifest = out / "template_manifest.tsv"
    with open(manifest, "w") as fh:
        fh.write("motif_idx\tpdb\tentry_chain_id\trow_name\tmotif\n")
        for k in kept:
            fh.write(f"{k['motif_idx']}\t{k['pdb']}\t{k['entry_chain_id']}\t{k['row_name']}\t{k['motif']}\n")

    print(f"[OK] Prepared {len(kept)} motif templates for {rna}")

if __name__ == "__main__":
    main()
