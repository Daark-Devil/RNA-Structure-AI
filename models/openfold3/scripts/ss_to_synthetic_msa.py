#!/usr/bin/env python3
"""
ss_to_synthetic_msa.py  v5 — Structure-aware synthetic MSA generator.

Adds 4 awareness features over v4:
  1. Length-tier: n_seqs auto-set from sequence length
  2. Stem-length aware: longer stems are more conserved
  3. Stem-position aware: stem middles conserved, ends variable
  4. Loop-type aware: hairpin/internal/bulge/junction/tail get different rates

Backward compatible: pass --legacy for original static-parameter behavior.
"""

import argparse
import random
from pathlib import Path


# ============================================================================
# CONFIG — biology-driven parameters from Rfam/CMfinder studies
# ============================================================================

PAIR_PROBS = {
    'A': {'U': 0.85, 'G': 0.05, 'C': 0.05, 'A': 0.05},
    'U': {'A': 0.70, 'G': 0.25, 'C': 0.03, 'U': 0.02},
    'G': {'C': 0.75, 'U': 0.20, 'A': 0.03, 'G': 0.02},
    'C': {'G': 0.85, 'A': 0.05, 'U': 0.05, 'C': 0.05},
}

# Legacy static rates (used with --legacy flag)
LEGACY_MUT_PAIRED   = 0.15
LEGACY_MUT_UNPAIRED = 0.25
LEGACY_N_SEQS       = 256

# Gap rates (kept stable for both modes)
GAP_RATE_LOOP = 0.05
GAP_RATE_STEM = 0.005

BRACKET_LEVELS = [('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]
NAME_PAD = 40


# ============================================================================
# Length-tier (sets n_seqs)
# ============================================================================
def get_n_seqs_for_length(length: int) -> int:
    if length < 30:    return 256
    elif length < 50:  return 384
    elif length < 85:  return 512
    elif length < 100: return 768
    elif length < 140: return 1024
    else:              return 1500


# ============================================================================
# Stem detection — find consecutive base pairs grouped into stems
# ============================================================================
def detect_stems(pairs):
    """Group consecutive base pairs (i,j), (i+1,j-1), ... into stems.
    Returns list of stems, each stem is list of (i,j) tuples in order."""
    if not pairs:
        return []
    pairs_sorted = sorted(pairs)
    stems = []
    current_stem = [pairs_sorted[0]]
    for prev, cur in zip(pairs_sorted, pairs_sorted[1:]):
        # consecutive in stem: i increments by 1, j decrements by 1
        if cur[0] == prev[0] + 1 and cur[1] == prev[1] - 1:
            current_stem.append(cur)
        else:
            stems.append(current_stem)
            current_stem = [cur]
    stems.append(current_stem)
    return stems


def get_pair_paired_mut_rate(pair, stem):
    """Stem-length × stem-position aware mutation rate for a paired position."""
    stem_len = len(stem)
    # Stem-length base rate
    if stem_len >= 5:   base = 0.08
    elif stem_len >= 3: base = 0.15
    else:               base = 0.25

    # Stem-position multiplier (closing pairs flex more)
    idx = stem.index(pair)
    dist_from_end = min(idx, stem_len - 1 - idx)
    if dist_from_end == 0:   mult = 2.0  # closing pair (very end)
    elif dist_from_end <= 1: mult = 1.5  # near-end
    else:                    mult = 1.0  # middle

    return min(0.50, base * mult)  # cap at 50% to avoid extreme


# ============================================================================
# Loop classification
# ============================================================================
def classify_unpaired_positions(db, pairs):
    """Classify each unpaired position into a loop type:
       'tail5', 'tail3', 'hairpin', 'internal', 'bulge', 'junction'.
    Returns dict {position: loop_type}."""
    n = len(db)
    paired_set = set()
    for i, j in pairs:
        paired_set.add(i)
        paired_set.add(j)

    # Find leftmost and rightmost paired positions
    if not paired_set:
        return {i: 'tail5' for i in range(n)}

    first_paired = min(paired_set)
    last_paired = max(paired_set)

    classification = {}
    for i in range(n):
        if i in paired_set:
            continue
        if i < first_paired:
            classification[i] = 'tail5'
            continue
        if i > last_paired:
            classification[i] = 'tail3'
            continue
        # Internal unpaired - determine loop type by context
        classification[i] = _classify_internal_loop(i, db, pairs)
    return classification


def _classify_internal_loop(pos, db, pairs):
    """Determine if unpaired position is in hairpin, internal, bulge, or junction."""
    pair_dict = {}
    for i, j in pairs:
        pair_dict[i] = j
        pair_dict[j] = i

    # Walk left until we find a paired position
    left_pair_pos = None
    for k in range(pos - 1, -1, -1):
        if k in pair_dict:
            left_pair_pos = k
            break

    # Walk right until we find a paired position
    right_pair_pos = None
    for k in range(pos + 1, len(db)):
        if k in pair_dict:
            right_pair_pos = k
            break

    if left_pair_pos is None or right_pair_pos is None:
        return 'tail5'  # safety

    left_partner = pair_dict[left_pair_pos]
    right_partner = pair_dict[right_pair_pos]

    # Hairpin: left and right paired positions are partners of each other
    if left_partner == right_pair_pos:
        return 'hairpin'

    # Internal/bulge: left_partner is just before right_partner in pair_dict
    # i.e. the two stems on either side of this loop close to form a 2-way junction
    # Check for symmetric internal loop, asymmetric, or bulge
    # by checking how many unpaired residues are between the two stems
    left_unpaired_count = right_pair_pos - left_pair_pos - 1

    # Right side: count unpaired between left_partner and right_partner
    if abs(left_partner - right_partner) <= 1:
        # adjacent stems = single junction point, this is bulge
        return 'bulge' if left_unpaired_count <= 3 else 'internal'

    right_unpaired_count = abs(left_partner - right_partner) - 1

    # Now we have left_unpaired and right_unpaired counts
    if right_unpaired_count == 0:
        return 'bulge' if left_unpaired_count <= 3 else 'internal'
    if left_unpaired_count <= 3 and right_unpaired_count <= 3:
        # short loop on both sides
        if abs(left_unpaired_count - right_unpaired_count) == 0:
            return 'internal'  # symmetric
        return 'internal'
    if left_unpaired_count > 0 and right_unpaired_count > 0:
        return 'internal'

    return 'junction'  # fallback for complex cases


def get_loop_unpaired_mut_rate(loop_type, surrounding_length=None):
    """Loop-type aware mutation rate for an unpaired position."""
    rates = {
        'hairpin':  0.22,   # average of 0.20 (tetraloop) and 0.25 (longer)
        'internal': 0.32,   # average of symmetric 0.30 and asymmetric 0.35
        'bulge':    0.40,
        'junction': 0.20,
        'tail5':    0.45,
        'tail3':    0.45,
    }
    return rates.get(loop_type, 0.30)


# ============================================================================
# SS parsing
# ============================================================================
def parse_dot_bracket(db: str):
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
            pairs.append((j, i))
    for lvl, st in enumerate(stacks):
        if st:
            ob, cb = BRACKET_LEVELS[lvl]
            raise ValueError(f"Unmatched '{ob}' at positions {st}")
    return sorted(pairs)


def sample_complement(nt, rng):
    probs = PAIR_PROBS.get(nt.upper())
    if probs is None:
        return rng.choice(['A', 'U', 'G', 'C'])
    r = rng.random()
    cum = 0
    for partner, p in probs.items():
        cum += p
        if r <= cum:
            return partner
    return list(probs.keys())[-1]


# ============================================================================
# Generate one synthetic sequence
# ============================================================================
def generate_one_sequence(query, pairs, stems, unpaired_class, paired_pos, rng, legacy=False):
    """Generate one synthetic sequence. Position-by-position decision."""
    n = len(query)
    seq = list(query)

    # Build pair-to-stem lookup
    pair_to_stem = {}
    for stem in stems:
        for pair in stem:
            pair_to_stem[pair] = stem

    processed = set()
    for i, j in pairs:
        if i in processed or j in processed:
            continue
        processed.add(i); processed.add(j)

        # Decide mutation rate
        if legacy:
            mut_rate = LEGACY_MUT_PAIRED
        else:
            stem = pair_to_stem[(i, j)]
            mut_rate = get_pair_paired_mut_rate((i, j), stem)

        # Mutate i and pick complementary j
        if rng.random() < mut_rate:
            new_i = rng.choice(['A', 'U', 'G', 'C'])
            new_j = sample_complement(new_i, rng)
            seq[i] = new_i
            seq[j] = new_j

    # Handle unpaired positions
    for k in range(n):
        if k in paired_pos:
            continue
        if legacy:
            mut_rate = LEGACY_MUT_UNPAIRED
        else:
            loop_type = unpaired_class.get(k, 'internal')
            mut_rate = get_loop_unpaired_mut_rate(loop_type)

        if rng.random() < mut_rate:
            seq[k] = rng.choice(['A', 'U', 'G', 'C'])

    # Apply gaps
    for k in range(n):
        gap_prob = GAP_RATE_STEM if k in paired_pos else GAP_RATE_LOOP
        if rng.random() < gap_prob:
            seq[k] = '-'

    return ''.join(seq)


# ============================================================================
# Stockholm output (matches OF3 format)
# ============================================================================
def write_stockholm(query_name, query_seq, synth_seqs, out_path):
    n_seq = len(synth_seqs)
    pp_string = '9' * len(query_seq)
    pp_cons_string = '9' * len(query_seq)
    rf_string = 'x' * len(query_seq)

    name_query = query_name.ljust(NAME_PAD)
    name_pp_cons = '#=GC PP_cons'.ljust(NAME_PAD)
    name_rf = '#=GC RF'.ljust(NAME_PAD)

    lines = ['# STOCKHOLM 1.0', '']
    lines.append(f"#=GS {query_name.ljust(NAME_PAD - 5)}AC {query_name}")
    lines.append('')
    lines.append(f"#=GS {query_name.ljust(NAME_PAD - 5)}DE {query_name}")
    lines.append('')
    lines.append(f"{name_query} {query_seq}")

    for idx, sseq in enumerate(synth_seqs, 1):
        synth_name = f"synth{idx:04d}/1-{len(query_seq)}"
        name_padded = synth_name.ljust(NAME_PAD)
        lines.append(f"{name_padded} {sseq}")
        pp_name = f"#=GR {synth_name} PP".ljust(NAME_PAD)
        lines.append(f"{pp_name} {pp_string}")

    lines.append(f"{name_pp_cons} {pp_cons_string}")
    lines.append(f"{name_rf} {rf_string}")
    lines.append('//')

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ============================================================================
# Main
# ============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq", required=True, help="Path to .seq file")
    ap.add_argument("--ss", required=True, help="Path to .2d file (dot-bracket)")
    ap.add_argument("--out", required=True, help="Output Stockholm path")
    ap.add_argument("--query_name", required=True, help="Query name (e.g. 1NBS_A)")
    ap.add_argument("--n_seqs", type=int, default=None,
                    help="Number of synthetic seqs (default: auto from length)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--legacy", action="store_true",
                    help="Use original static parameters (no awareness)")
    args = ap.parse_args()

    # Read query
    query_seq = Path(args.seq).read_text().strip().split()[0].upper()
    db_text = Path(args.ss).read_text().strip().split('\n')
    db = db_text[1] if len(db_text) >= 2 else db_text[0]

    if len(db) != len(query_seq):
        print(f"[warn] SS length {len(db)} != seq length {len(query_seq)}")

    pairs = parse_dot_bracket(db)
    paired_pos = set()
    for i, j in pairs:
        paired_pos.add(i); paired_pos.add(j)

    # Decide n_seqs
    if args.n_seqs is not None:
        n_seqs = args.n_seqs
    elif args.legacy:
        n_seqs = LEGACY_N_SEQS
    else:
        n_seqs = get_n_seqs_for_length(len(query_seq))

    # Build structure context (only if not legacy)
    stems = []
    unpaired_class = {}
    if not args.legacy:
        stems = detect_stems(pairs)
        unpaired_class = classify_unpaired_positions(db, pairs)
        # Print summary
        loop_counts = {}
        for v in unpaired_class.values():
            loop_counts[v] = loop_counts.get(v, 0) + 1
        stem_lens = [len(s) for s in stems]
        print(f"[info] sequence length : {len(query_seq)}")
        print(f"[info] base pairs      : {len(pairs)}")
        print(f"[info] stems detected  : {len(stems)}  (lengths: {stem_lens})")
        print(f"[info] loop classes    : {loop_counts}")
        print(f"[info] n_seqs (auto)   : {n_seqs}")
    else:
        print(f"[info] LEGACY mode — static parameters")
        print(f"[info] sequence length : {len(query_seq)}")
        print(f"[info] base pairs      : {len(pairs)}")
        print(f"[info] n_seqs          : {n_seqs}")

    print(f"[info] generating {n_seqs} synthetic sequences (seed={args.seed})")
    rng = random.Random(args.seed)
    synth = [generate_one_sequence(query_seq, pairs, stems, unpaired_class,
                                    paired_pos, rng, legacy=args.legacy)
             for _ in range(n_seqs)]

    write_stockholm(args.query_name, query_seq, synth, args.out)
    file_size = Path(args.out).stat().st_size
    print(f"[done] wrote: {args.out}")
    print(f"[done] file size: {file_size} bytes")


if __name__ == "__main__":
    main()
