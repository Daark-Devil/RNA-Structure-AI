#!/usr/bin/env python3
"""
generate_boltz_yaml.py — Boltz-2 YAML generator for full RNAs.

Matches Imdad's pipeline structure:
  - Single chain, sequence only
  - SS base pairs as contact constraints
  - NO templates (Boltz-2 doesn't support RNA templates)

Usage:
  python generate_boltz_yaml.py \
      --rna 1NBS \
      --seq /path/to/1NBS.seq \
      --ss  /path/to/1NBS.2d \
      --out /path/to/output.yaml
"""

import argparse
from pathlib import Path

BRACKET_LEVELS = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]


def parse_dot_bracket(db):
    open_chars = {ob: lvl for lvl, (ob, cb) in enumerate(BRACKET_LEVELS)}
    close_chars = {cb: lvl for lvl, (ob, cb) in enumerate(BRACKET_LEVELS)}
    stacks = [[] for _ in BRACKET_LEVELS]
    pairs = []
    for i, ch in enumerate(db):
        if ch == '.':
            continue
        if ch in open_chars:
            stacks[open_chars[ch]].append(i)
        elif ch in close_chars:
            lvl = close_chars[ch]
            if not stacks[lvl]:
                raise ValueError(f"Unmatched '{ch}' at position {i}")
            j = stacks[lvl].pop()
            pairs.append((min(i, j), max(i, j)))
    for lvl, st in enumerate(stacks):
        if st:
            ob, cb = BRACKET_LEVELS[lvl]
            raise ValueError(f"Unmatched '{ob}'")
    return sorted(pairs)


def classify_pair(nt_i, nt_j):
    wc = {('A','U'),('U','A'),('G','C'),('C','G')}
    wobble = {('G','U'),('U','G')}
    pair = (nt_i.upper(), nt_j.upper())
    if pair in wc: return 'WC'
    if pair in wobble: return 'wobble'
    return 'noncanonical'


def get_max_distance(pair_type):
    # Same defaults as Imdad's pipeline (4.0 Å contact distance)
    return {'WC': 4.0, 'wobble': 4.5, 'noncanonical': 6.0}.get(pair_type, 5.0)


def should_force(pair_type):
    return pair_type in ('WC', 'wobble')


def write_yaml(out_path, sequence, pairs, pair_types, chain_id='A'):
    lines = ["version: 1", "", "sequences:",
             f"  - rna:",
             f"      id: {chain_id}",
             f"      sequence: {sequence}",
             ""]
    
    if pairs:
        lines.append("constraints:")
        for (i, j), pt in zip(pairs, pair_types):
            max_dist = get_max_distance(pt)
            force = should_force(pt)
            lines.append("  - contact:")
            lines.append(f"      token1: [{chain_id}, {i+1}]")
            lines.append(f"      token2: [{chain_id}, {j+1}]")
            lines.append(f"      max_distance: {max_dist}")
            lines.append(f"      force: {str(force).lower()}")
    
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--seq", required=True)
    ap.add_argument("--ss", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain_id", default="A")
    ap.add_argument("--no_constraints", action="store_true",
                    help="Generate YAML without constraints (for baseline)")
    args = ap.parse_args()

    sequence = Path(args.seq).read_text().strip().split()[0].upper()
    db_text = Path(args.ss).read_text().strip().split('\n')
    db = db_text[1] if len(db_text) >= 2 else db_text[0]
    
    if len(db) != len(sequence):
        print(f"[warn] SS length {len(db)} != seq length {len(sequence)}")
    
    if args.no_constraints:
        pairs = []
        pair_types = []
    else:
        pairs = parse_dot_bracket(db)
        pair_types = [classify_pair(sequence[i], sequence[j]) for (i, j) in pairs]
    
    pair_counts = {}
    for pt in pair_types:
        pair_counts[pt] = pair_counts.get(pt, 0) + 1
    
    print(f"[info] RNA          : {args.rna}")
    print(f"[info] length       : {len(sequence)}")
    print(f"[info] base pairs   : {len(pairs)}")
    if pair_counts:
        print(f"[info] pair types   : {pair_counts}")
    print(f"[info] mode         : {'no_constraints' if args.no_constraints else 'with_constraints'}")
    
    write_yaml(args.out, sequence, pairs, pair_types, args.chain_id)
    file_size = Path(args.out).stat().st_size
    print(f"[done] wrote: {args.out}")
    print(f"[done] file size: {file_size} bytes")


if __name__ == "__main__":
    main()
