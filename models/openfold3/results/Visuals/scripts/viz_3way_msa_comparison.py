#!/usr/bin/env python3
"""
viz_3way_msa_comparison.py  (CLI version)

For each RNA, generates a 4-panel PNG:
  1. Baseline (no MSA)      best prediction
  2. With Synthetic MSA     best prediction
  3. Reference (ground truth)
  4. All three overlaid     (ref=black, baseline=red, MSA=green)

Best structure per condition is chosen by lowest C1' Kabsch RMSD vs the reference.
Requires PyMOL (cartoon rendering) and PIL (panel assembly).

Usage:
  python viz_3way_msa_comparison.py \
      --rna 1E7K \
      --workdirs   /path/to/rna_workdirs \
      --baseline-base /path/to/batch_runs/withtemplate \
      --msa-base   /path/to/ss_msa_test \
      --out-dir    /path/to/output \
      --rmsd-tool  /path/to/compare_rna_kabsch.py

Multiple RNAs: --rna 1E7K 1P5P 1P5O
"""
import os
import csv
import argparse
import subprocess
import pymol
from pymol import cmd

MOTIF_COLORS = {
    "HELIX": "marine", "HAIRPIN": "firebrick", "2-WAY": "tv_orange",
    "3-WAY": "purple", "4-WAY": "tv_green", "TAIL": "gray60",
    "OPEN": "hotpink", "OPENLOOP": "hotpink", "KISS": "cyan", "PK": "tv_yellow",
}
PANEL_W, PANEL_H = 1000, 1100
DPI = 150
ZOOM_BUFFER = 8

# These get filled from CLI args in main()
CFG = {}


def setup_pymol_aesthetics():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.3)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_ladder_mode", 1)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("cartoon_tube_radius", 0.5)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.1)


def read_manifest(rna):
    path = f"{CFG['workdirs']}/{rna}/markers/{rna}_manifest.tsv"
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def read_marker(rna, fasta_file):
    path = f"{CFG['workdirs']}/{rna}/markers/{fasta_file}"
    if not os.path.exists(path):
        return ""
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(">"):
                seq.append(line)
    return "".join(seq)


def get_motif_residues(marker):
    return [i + 1 for i, ch in enumerate(marker) if ch != '-']


def motif_color_for(name, type_):
    text = (name + " " + type_).upper()
    for key, color in MOTIF_COLORS.items():
        if key in text:
            return color
    return "white"


def get_rmsd(ref, pred):
    if not os.path.exists(CFG['rmsd_tool']) or not os.path.exists(pred):
        return None
    try:
        out = subprocess.run(
            ["python", CFG['rmsd_tool'], "--ref", ref, "--pred", pred],
            capture_output=True, text=True, timeout=60)
        for line in out.stdout.splitlines():
            if "C1'_KABSCH_RMSD" in line:
                return float(line.split("=")[1].strip().split()[0])
    except Exception:
        pass
    return None


def find_best_cif(search_dir, ref):
    best_rmsd = float('inf')
    best_cif = None
    if not os.path.isdir(search_dir):
        return None, None
    for root, _, files in os.walk(search_dir):
        for f in files:
            if f.endswith('_model.cif'):
                path = os.path.join(root, f)
                r = get_rmsd(ref, path)
                if r is not None and r < best_rmsd:
                    best_rmsd = r
                    best_cif = path
    return best_cif, best_rmsd


def render_single_panel(struct_path, struct_name, motif_data, out_png, ref_for_align=None):
    cmd.reinitialize()
    setup_pymol_aesthetics()
    cmd.load(struct_path, struct_name)
    if ref_for_align is not None:
        cmd.load(ref_for_align, "_align_ref")
        try:
            cmd.align(f"{struct_name} and name C1'", "_align_ref and name C1'", cycles=0)
        except Exception:
            pass
        cmd.delete("_align_ref")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.color("white", struct_name)
    for name, type_, residues, color in motif_data:
        if not residues:
            continue
        resi_str = "+".join(str(r) for r in residues)
        sel = f"sel_{name}".replace("-", "_").replace(" ", "_")
        cmd.select(sel, f"{struct_name} and resi {resi_str}")
        cmd.color(color, sel)
        cmd.deselect()
    cmd.zoom(struct_name, buffer=ZOOM_BUFFER)
    cmd.ray(PANEL_W, PANEL_H)
    cmd.png(out_png, dpi=DPI)


def render_overlay_panel(baseline_path, msa_path, ref_path, motif_data, out_png):
    cmd.reinitialize()
    setup_pymol_aesthetics()
    cmd.load(ref_path, "ref")
    cmd.load(baseline_path, "baseline")
    cmd.load(msa_path, "msa")
    for mob in ("baseline", "msa"):
        try:
            cmd.align(f"{mob} and name C1'", "ref and name C1'", cycles=0)
        except Exception:
            pass
    cmd.hide("everything")
    for s in ("ref", "baseline", "msa"):
        cmd.show("cartoon", s)
    cmd.color("black", "ref")
    cmd.color("red", "baseline")
    cmd.color("green", "msa")
    cmd.set("cartoon_transparency", 0.0, "ref")
    cmd.set("cartoon_transparency", 0.4, "baseline")
    cmd.set("cartoon_transparency", 0.4, "msa")
    cmd.zoom("ref", buffer=ZOOM_BUFFER)
    cmd.ray(PANEL_W, PANEL_H)
    cmd.png(out_png, dpi=DPI)


def visualize_one(rna):
    print(f"\n{'-' * 60}\n  {rna}\n{'-' * 60}")

    ref_pdb = f"{CFG['workdirs']}/{rna}/{rna}.pdb"
    if not os.path.exists(ref_pdb):
        print("  REF MISSING"); return False

    baseline_dir = f"{CFG['baseline_base']}/{rna}/output_of3-p2-145k/{rna}_motif_test"
    baseline_cif, baseline_rmsd = find_best_cif(baseline_dir, ref_pdb)
    if baseline_cif is None:
        print("  BASELINE NOT FOUND"); return False
    print(f"  Baseline best: {baseline_rmsd:.3f} A")

    msa_cif = msa_rmsd = None
    for sub in ['of3_output_msa_full', 'of3_output_msa', 'of3_output_msa_with_lma']:
        sd = f"{CFG['msa_base']}/{rna}/{sub}"
        if os.path.isdir(sd):
            msa_cif, msa_rmsd = find_best_cif(sd, ref_pdb)
            if msa_cif:
                print(f"  MSA best:      {msa_rmsd:.3f} A"); break
    if msa_cif is None:
        print("  MSA NOT FOUND"); return False

    print(f"  Improvement:   {baseline_rmsd - msa_rmsd:+.2f} A")

    motif_data = []
    for entry in read_manifest(rna):
        marker = read_marker(rna, entry['fasta_file'])
        residues = get_motif_residues(marker) if marker else []
        if residues:
            motif_data.append((entry['motif_name'], entry['motif_type'], residues,
                               motif_color_for(entry['motif_name'], entry['motif_type'])))

    out_rna = f"{CFG['out_dir']}/{rna}"
    os.makedirs(out_rna, exist_ok=True)
    p1, p2, p3, p4 = (f"{out_rna}/_p{i}.png" for i in (1, 2, 3, 4))

    print("  [1/4] Baseline panel..."); render_single_panel(baseline_cif, "base", motif_data, p1, ref_pdb)
    print("  [2/4] MSA panel...");      render_single_panel(msa_cif, "msa_struct", motif_data, p2, ref_pdb)
    print("  [3/4] Reference panel..."); render_single_panel(ref_pdb, "ref_struct", motif_data, p3, None)
    print("  [4/4] Overlay panel...");  render_overlay_panel(baseline_cif, msa_cif, ref_pdb, motif_data, p4)

    from PIL import Image, ImageDraw, ImageFont
    imgs = [Image.open(p) for p in (p1, p2, p3, p4)]
    pw = max(i.width for i in imgs)
    ph = max(i.height for i in imgs)
    label_h = 60
    combined = Image.new("RGB", (2 * pw, 2 * (ph + label_h)), "white")
    draw = ImageDraw.Draw(combined)
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", 26)
    except Exception:
        font = ImageFont.load_default()
    combined.paste(imgs[0], (0, label_h))
    combined.paste(imgs[1], (pw, label_h))
    combined.paste(imgs[2], (0, ph + 2 * label_h))
    combined.paste(imgs[3], (pw, ph + 2 * label_h))
    labels = [
        (f"{rna} - Baseline (no MSA)  RMSD {baseline_rmsd:.2f} A", 0, pw, 15),
        (f"{rna} - With Synthetic MSA  RMSD {msa_rmsd:.2f} A", pw, 2 * pw, 15),
        (f"{rna} - Reference (ground truth)", 0, pw, ph + label_h + 15),
        (f"{rna} - Overlay: Ref(black) Base(red) MSA(green)", pw, 2 * pw, ph + label_h + 15),
    ]
    for text, xs, xe, y in labels:
        try:
            tw = draw.textlength(text, font=font)
        except Exception:
            tw = len(text) * 14
        draw.text((xs + ((xe - xs) - tw) // 2, y), text, fill="black", font=font)
    out_path = f"{out_rna}/3way_msa_comparison.png"
    combined.save(out_path)
    print(f"  saved {out_path}")
    for p in (p1, p2, p3, p4):
        try: os.remove(p)
        except Exception: pass
    return True


def main():
    ap = argparse.ArgumentParser(description="4-panel baseline/MSA/reference/overlay RNA comparison")
    ap.add_argument("--rna", nargs="+", required=True, help="one or more RNA names")
    ap.add_argument("--workdirs", required=True, help="dir with <RNA>/<RNA>.pdb + markers/")
    ap.add_argument("--baseline-base", required=True, help="batch_runs/withtemplate base dir")
    ap.add_argument("--msa-base", required=True, help="ss_msa_test base dir")
    ap.add_argument("--out-dir", required=True, help="output dir (per-RNA subfolders created)")
    ap.add_argument("--rmsd-tool", required=True, help="path to compare_rna_kabsch.py")
    args = ap.parse_args()

    CFG.update(workdirs=args.workdirs, baseline_base=args.baseline_base,
               msa_base=args.msa_base, out_dir=args.out_dir, rmsd_tool=args.rmsd_tool)
    os.makedirs(args.out_dir, exist_ok=True)

    pymol.finish_launching(['pymol', '-cq'])
    print("4-panel: Baseline | MSA | Reference | Overlay")

    ok = [rna for rna in args.rna if visualize_one(rna)]
    print(f"\nDONE - {len(ok)}/{len(args.rna)} created")
    for rna in ok:
        print(f"  {args.out_dir}/{rna}/3way_msa_comparison.png")


if __name__ == "__main__":
    main()
