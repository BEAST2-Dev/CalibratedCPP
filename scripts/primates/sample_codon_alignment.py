#!/usr/bin/env python3
"""Sample N protein-coding loci at random and concatenate them into an alignment
that is safe to run with codon partitioning.

This follows the Vanderpool et al. (2020) "molecular dating" recipe -- draw loci
without replacement from their 1730 and concatenate -- but keeps the result IN
READING FRAME, which their published dating datasets are not.

Why it matters
--------------
Charsets like "1-.\\3" select columns by index mod 3. That equals codon position
only if every locus in the concatenation is a whole number of codons AND starts
on a 1st position. In the Dryad locus alignments neither holds universally:

  * 246 of the 1730 loci have a length that is not a multiple of 3;
  * ~13% do not start on a 1st codon position -- the "NoNcol.Noambig" filtering
    dropped individual columns, not whole codons;
  * a minority change frame *inside* the locus, for the same reason.

Concatenate those unchanged and each codon charset becomes a ~1/3 mixture of all
three codon positions, so the partitions come back with near-identical rates
(r1 ~ r2 ~ r3) and near-identical base frequencies. So every locus is screened
before it is accepted:

  1. its start phase is detected and the leading 0-2 bases dropped,
  2. the trailing partial codon is dropped,
  3. loci that change frame mid-sequence are rejected and replaced.

Frame is called from the data: 3rd codon positions are mostly synonymous and so
are much more variable than 1st and 2nd, and whichever offset maximises mean
column variability marks the 3rd positions. The majority start codon (ATG) is
used only as a fallback for loci too invariant to call -- it is not reliable on
its own (15 of 167 ATG-starting loci in this pool are out of frame anyway).

The output is verified before the script exits, by two independent checks:
translating in all three frames and counting internal stop codons, and the
variability-by-offset test above.

Input
-----
The per-locus FASTA alignments from the Dryad deposit, either as the tarball or
as an already-extracted directory:
    doi_10_5061_dryad_rfj6q577d / 1730_Alignments_FINAL.tar.gz

Usage
-----
    ./sample_codon_alignment.py <loci.tar.gz | locus_dir> <n_loci> <out.nex>
                                [--seed S] [--min-taxa K] [--partitions out.txt]

Example (400 loci, reproducible):
    ./sample_codon_alignment.py ~/Downloads/doi_10_5061_dryad_rfj6q577d__v20201103/\\
1730_Alignments_FINAL.tar.gz 400 ../../calibratedcpp-lphy/examples/data/primates_400loci.nex \\
        --seed 1 --partitions ../../calibratedcpp-lphy/examples/data/primates_400loci_partitions.txt
"""
import argparse
import glob
import os
import random
import sys
import tarfile
from collections import Counter

import numpy as np

MISSING = set("-N?")
STOPS = {"TAA", "TAG", "TGA"}
MIN_RATIO = 1.5        # variability ratio below which a locus cannot be called
QC_WINDOW = 300        # window for the internal frame-consistency check


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def parse_fasta(lines):
    seqs, name = {}, None
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    seqs = {k: "".join(v).upper() for k, v in seqs.items()}
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f"not aligned: lengths {sorted(lengths)}")
    return seqs, lengths.pop()


def read_pool(src, pattern):
    """Yield (name, seqs, length) from a directory or a .tar.gz of FASTAs."""
    if os.path.isdir(src):
        files = sorted(glob.glob(os.path.join(src, pattern)))
        if not files:
            sys.exit(f"no {pattern} files under {src}")
        for f in files:
            seqs, L = parse_fasta(open(f))
            yield os.path.basename(f), seqs, L
    else:
        with tarfile.open(src) as tf:
            members = sorted((m for m in tf.getmembers()
                              if m.isfile() and m.name.endswith((".fa", ".fasta"))),
                             key=lambda m: m.name)
            if not members:
                sys.exit(f"no FASTA members inside {src}")
            for m in members:
                text = tf.extractfile(m).read().decode().splitlines()
                seqs, L = parse_fasta(text)
                yield os.path.basename(m.name), seqs, L


# --------------------------------------------------------------------------
# frame detection
# --------------------------------------------------------------------------
def column_variability(seqs, length):
    """Per column: 1 - frequency of the commonest non-missing base."""
    mat = np.array([list(s) for s in seqs.values()])
    var = np.zeros(length)
    for j in range(length):
        col = [c for c in mat[:, j] if c not in MISSING]
        var[j] = 0.0 if len(col) < 4 else 1.0 - Counter(col).most_common(1)[0][1] / len(col)
    return var


def call_offset(var):
    """Offset (0/1/2) holding the 3rd codon positions, or None if uncallable."""
    means = [var[k::3].mean() for k in range(3)]
    if min(means) <= 0 or max(means) / min(means) < MIN_RATIO:
        return None
    return int(np.argmax(means))


def screen_locus(seqs, length):
    """Return (leading bases to drop, how) or (None, reason) if unusable."""
    var = column_variability(seqs, length)
    offset = call_offset(var)

    if offset is None:
        # Too little variation to call. Trust the start codon if it is there.
        start = Counter(s[:3] for s in seqs.values()).most_common(1)[0][0]
        return 0, "weak/ATG" if start == "ATG" else "weak"

    # The frame must hold all the way along; column-level filtering can shift it
    # mid-locus, and no amount of end-trimming repairs that.
    calls = [c for c in (call_offset(var[s:s + QC_WINDOW])
                         for s in range(0, length - QC_WINDOW + 1, QC_WINDOW))
             if c is not None]
    if len(set(calls)) > 1:
        return None, "internal frame break"

    # 3rd positions sit at `offset`, so 1st positions sit at offset + 1.
    return (offset + 1) % 3, "variability"


# --------------------------------------------------------------------------
# verification of the finished alignment
# --------------------------------------------------------------------------
def check_stop_codons(names, seqs_out):
    """Translate in all three frames and count internal stops.

    The correct frame sits near 0%; a wrong frame scrambles every codon and
    stops turn up at random, ~4-5%. If all three frames land at a similar
    middling value, the frame CHANGES partway along the alignment.
    """
    ref = "Homo_sapiens" if "Homo_sapiens" in seqs_out else names[0]
    s = seqs_out[ref]
    rates = []
    for off in range(3):
        cods = (s[j:j + 3] for j in range(off, len(s) - 2, 3))
        real = [c for c in cods if set(c) <= set("ACGT")]
        rates.append(100.0 * sum(c in STOPS for c in real) / max(len(real), 1))
    best = int(np.argmin(rates))
    others = min(r for k, r in enumerate(rates) if k != best)
    ok = best == 0 and rates[0] < 1.0 and others > 2 * rates[0]
    print(f"  stop codons in {ref}: " + "  ".join(
        f"frame {k}={r:.3f}%" for k, r in enumerate(rates)))
    return ok, rates


def check_variability(seqs_out, charsets):
    """Variability by offset, globally and per locus."""
    mat = np.array([list(s) for s in seqs_out.values()])
    n = mat.shape[1]
    var = np.zeros(n)
    for j in range(n):
        col = [c for c in mat[:, j] if c not in MISSING]
        var[j] = 0.0 if len(col) < 4 else 1.0 - Counter(col).most_common(1)[0][1] / len(col)

    means = [var[k::3].mean() for k in range(3)]
    ratio = max(means) / min(means) if min(means) > 0 else float("inf")
    print(f"  variability by offset: " + "  ".join(f"{m:.5f}" for m in means)
          + f"   (3rd positions at offset {int(np.argmax(means))}, ratio {ratio:.2f})")

    good = bad = weak = 0
    for _, a, b in charsets:
        c = call_offset(var[a - 1:b])
        if c is None:
            weak += 1
        elif c == 2:
            good += 1
        else:
            bad += 1
    print(f"  per-locus frame: {good} in frame, {bad} out of frame, {weak} uncallable")
    return bad == 0 and int(np.argmax(means)) == 2


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("loci", help="directory of per-locus FASTAs, or the .tar.gz")
    ap.add_argument("n_loci", type=int)
    ap.add_argument("out_nex")
    ap.add_argument("--seed", type=int, default=None,
                    help="RNG seed; set it to make the sample reproducible")
    ap.add_argument("--min-taxa", type=int, default=None,
                    help="skip loci with fewer taxa than this "
                         "(default: only loci with every taxon present)")
    ap.add_argument("--partitions", default=None,
                    help="also write a charset/partition file here")
    ap.add_argument("--pattern", default="*.fa")
    ap.add_argument("--no-screen", action="store_true",
                    help="concatenate as-is, skipping the frame screening "
                         "(reproduces the original, broken behaviour)")
    args = ap.parse_args()

    pool, all_taxa = [], set()
    for name, seqs, L in read_pool(args.loci, args.pattern):
        pool.append((name, seqs, L))
        all_taxa |= set(seqs)
    all_taxa = sorted(all_taxa)
    min_taxa = args.min_taxa if args.min_taxa is not None else len(all_taxa)

    eligible = [p for p in pool if len(p[1]) >= min_taxa]
    print(f"pool: {len(pool)} loci, {len(all_taxa)} taxa total")
    print(f"  loci with >= {min_taxa} taxa: {len(eligible)}")
    if len(eligible) < args.n_loci:
        sys.exit(f"ERROR: only {len(eligible)} eligible loci, asked for {args.n_loci}")

    # --- sample without replacement, screening as we go ---------------------
    rng = random.Random(args.seed)
    order = eligible[:]
    rng.shuffle(order)

    chosen, how_counts, rejected = [], Counter(), Counter()
    for name, seqs, L in order:
        if len(chosen) == args.n_loci:
            break
        head, how = (0, "unscreened") if args.no_screen else screen_locus(seqs, L)
        if head is None:
            rejected[how] += 1
            continue
        how_counts[how] += 1
        chosen.append((name, seqs, L, head))

    if len(chosen) < args.n_loci:
        sys.exit(f"ERROR: only {len(chosen)} loci passed screening, "
                 f"asked for {args.n_loci}")
    chosen.sort(key=lambda p: p[0])          # deterministic concatenation order

    # --- concatenate --------------------------------------------------------
    parts = {t: [] for t in all_taxa}
    charsets, pos = [], 0
    shifted, trimmed, padded = Counter(), Counter(), Counter()
    for name, seqs, L, head in chosen:
        if head:
            shifted[head] += 1
        keep = (L - head) - ((L - head) % 3)
        if keep != L - head:
            trimmed[(L - head) % 3] += 1
        for t in all_taxa:
            s = seqs.get(t)
            if s is None:
                s = "-" * L                  # taxon absent from this locus
                padded[t] += 1
            parts[t].append(s[head:head + keep])
        charsets.append((name.split(".")[0], pos + 1, pos + keep))
        pos += keep

    nchar = pos
    assert nchar % 3 == 0
    seqs_out = {t: "".join(parts[t]) for t in all_taxa}
    assert all(len(s) == nchar for s in seqs_out.values())

    print(f"\nsampled {args.n_loci} loci (seed={args.seed}) -> "
          f"{nchar} bp x {len(all_taxa)} taxa")
    print(f"  loci screened to fill the sample: {len(chosen) + sum(rejected.values())}")
    for reason, k in rejected.items():
        print(f"    rejected, {reason}: {k}")
    print(f"  start phase called by: {dict(how_counts)}")
    print(f"  re-phased at the start: {sum(shifted.values())} loci "
          f"({sum(k * v for k, v in shifted.items())} leading bases dropped)")
    print(f"  trimmed at the end: {sum(trimmed.values())} loci "
          f"({sum(k * v for k, v in trimmed.items())} trailing bases dropped)")
    if padded:
        print(f"  gap-padded taxon slots: {sum(padded.values())} "
              f"across {len(padded)} taxa")

    # --- write --------------------------------------------------------------
    width = max(len(t) for t in all_taxa) + 2
    with open(args.out_nex, "w") as fh:
        fh.write("#NEXUS\n\nBEGIN DATA;\n")
        fh.write(f"\tDimensions ntax={len(all_taxa)} nchar={nchar};\n")
        fh.write("\tFormat datatype=nucleotide gap=-;\n\tMatrix\n")
        for t in all_taxa:
            fh.write(f"\t{t.ljust(width)}{seqs_out[t]}\n")
        fh.write(";\nEND;\n")
    print(f"  wrote {args.out_nex}")

    if args.partitions:
        with open(args.partitions, "w") as fh:
            fh.write(f"# {args.n_loci} loci sampled from {args.loci}, seed={args.seed}\n")
            fh.write("# Every locus is a whole number of codons and starts on a\n")
            fh.write("# 1st position, so codon position == column index mod 3:\n")
            fh.write(f"#   1st positions: 1-{nchar}\\3\n")
            fh.write(f"#   2nd positions: 2-{nchar}\\3\n")
            fh.write(f"#   3rd positions: 3-{nchar}\\3\n\n")
            for name, a, b in charsets:
                fh.write(f"charset {name} = {a}-{b};\n")
        print(f"  wrote {args.partitions}")

    # --- verify -------------------------------------------------------------
    print("\nverifying the output is safe for codon partitioning:")
    ok_stop, _ = check_stop_codons(all_taxa, seqs_out)
    ok_var = check_variability(seqs_out, charsets)
    if ok_stop and ok_var:
        print("  OK -- charsets \"1-.\\3\", \"2-.\\3\", \"3-.\\3\" are 1st, 2nd, "
              "3rd codon positions.")
        return 0
    print("  FAILED -- do NOT run this with codon partitioning.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
