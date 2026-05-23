#!/usr/bin/env python3
"""
viz_boltz_comparison.py  (CLI)

3-panel PNG comparing a Boltz-2 prediction against the native structure:
  1. Boltz prediction
  2. Native (reference)
  3. Overlay (native=black, prediction=green), aligned on C1'

Requires PyMOL and PIL.

Usage:
  python viz_boltz_comparison.py \
      --rna 3G8T \
      --pred /path/to/3G8T_model_0.cif \
      --native /path/to/3G8T.pdb \
      --out-dir /path/to/output \
      [--rmsd-tool /path/to/compare_rna_kabsch.py]
"""
import os
import argparse
import subprocess
import pymol
from pymol import cmd

PANEL_W, PANEL_H = 1000, 1100
DPI = 150
ZOOM_BUFFER = 8


def setup():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.3)
    cmd.set("antialias", 2)
    cmd.set("cartoon_ring_mode", 3)
    cmd.set("cartoon_ring_finder", 1)
    cmd.set("cartoon_nucleic_acid_mode", 4)
    cmd.set("cartoon_tube_radius", 0.5)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.1)


def get_rmsd(tool, ref, pred):
    if not tool or not os.path.exists(tool):
        return None
    try:
        out = subprocess.run(["python", tool, "--ref", ref, "--pred", pred],
                             capture_output=True, text=True, timeout=60)
        for line in out.stdout.splitlines():
            if "C1'_KABSCH_RMSD" in line:
                return float(line.split("=")[1].strip().split()[0])
    except Exception:
        pass
    return None


def single_panel(struct_path, name, out_png, color, ref_for_align=None):
    cmd.reinitialize()
    setup()
    cmd.load(struct_path, name)
    if ref_for_align is not None:
        cmd.load(ref_for_align, "_ref")
        try:
            cmd.align(f"{name} and name C1'", "_ref and name C1'", cycles=0)
        except Exception:
            pass
        cmd.delete("_ref")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.color(color, name)
    cmd.zoom(name, buffer=ZOOM_BUFFER)
    cmd.ray(PANEL_W, PANEL_H)
    cmd.png(out_png, dpi=DPI)


def overlay_panel(pred_path, native_path, out_png):
    cmd.reinitialize()
    setup()
    cmd.load(native_path, "native")
    cmd.load(pred_path, "pred")
    try:
        cmd.align("pred and name C1'", "native and name C1'", cycles=0)
    except Exception:
        pass
    cmd.hide("everything")
    cmd.show("cartoon", "native")
    cmd.show("cartoon", "pred")
    cmd.color("black", "native")
    cmd.color("green", "pred")
    cmd.set("cartoon_transparency", 0.0, "native")
    cmd.set("cartoon_transparency", 0.3, "pred")
    cmd.zoom("native", buffer=ZOOM_BUFFER)
    cmd.ray(PANEL_W, PANEL_H)
    cmd.png(out_png, dpi=DPI)


def main():
    ap = argparse.ArgumentParser(description="3-panel Boltz prediction vs native comparison")
    ap.add_argument("--rna", required=True)
    ap.add_argument("--pred", required=True, help="Boltz predicted CIF")
    ap.add_argument("--native", required=True, help="native PDB")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--rmsd-tool", default=None, help="optional compare_rna_kabsch.py")
    args = ap.parse_args()

    out_rna = os.path.join(args.out_dir, args.rna)
    os.makedirs(out_rna, exist_ok=True)

    rmsd = get_rmsd(args.rmsd_tool, args.native, args.pred)
    rmsd_str = f"{rmsd:.2f} A" if rmsd is not None else "n/a"
    print(f"{args.rna}: RMSD vs native = {rmsd_str}")

    pymol.finish_launching(['pymol', '-cq'])

    p1 = f"{out_rna}/_p1.png"
    p2 = f"{out_rna}/_p2.png"
    p3 = f"{out_rna}/_p3.png"
    print("  [1/3] prediction panel..."); single_panel(args.pred, "pred", p1, "green", ref_for_align=args.native)
    print("  [2/3] native panel...");     single_panel(args.native, "native", p2, "gray70", None)
    print("  [3/3] overlay panel...");    overlay_panel(args.pred, args.native, p3)

    from PIL import Image, ImageDraw, ImageFont
    imgs = [Image.open(p) for p in (p1, p2, p3)]
    pw = max(i.width for i in imgs)
    ph = max(i.height for i in imgs)
    label_h = 60
    combined = Image.new("RGB", (3 * pw, ph + label_h), "white")
    draw = ImageDraw.Draw(combined)
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", 26)
    except Exception:
        font = ImageFont.load_default()
    for idx, im in enumerate(imgs):
        combined.paste(im, (idx * pw, label_h))
    labels = [
        (f"{args.rna} - Boltz-2 prediction  RMSD {rmsd_str}", 0, pw),
        (f"{args.rna} - Native (reference)", pw, 2 * pw),
        (f"{args.rna} - Overlay: native(black) pred(green)", 2 * pw, 3 * pw),
    ]
    for text, xs, xe in labels:
        try:
            tw = draw.textlength(text, font=font)
        except Exception:
            tw = len(text) * 14
        draw.text((xs + ((xe - xs) - tw) // 2, 15), text, fill="black", font=font)
    out_path = f"{out_rna}/boltz_vs_native.png"
    combined.save(out_path)
    print(f"  saved {out_path}")
    for p in (p1, p2, p3):
        try: os.remove(p)
        except Exception: pass


if __name__ == "__main__":
    main()
