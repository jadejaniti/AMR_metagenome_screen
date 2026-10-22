#!/usr/bin/env python3
"""
Screen assembled metagenomic contigs for resistance determinants.

Aligns reference gene sequences (ABRicate-bundled databases) against contigs
with minimap2 (via mappy), applies identity/coverage thresholds comparable to
ABRicate/ResFinder defaults, and resolves overlapping hits.

Usage:
  python3 -I argscan.py --contigs CONTIGS.fsa_nt.gz --dbdir ABRICATE/db \
      --dbs card ncbi resfinder argannot megares bacmet2 plasmidfinder \
      --minid 80 --mincov 60 --out OUTDIR
"""
import argparse
import gzip
import json
import os
import sys
from collections import defaultdict

import mappy


def open_maybe_gzip(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def read_fasta(path):
    name, buf = None, []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def assembly_stats(path):
    lengths, gc, n_amb, total = [], 0, 0, 0
    for _, seq in read_fasta(path):
        s = seq.upper()
        lengths.append(len(s))
        gc += s.count("G") + s.count("C")
        n_amb += s.count("N")
        total += len(s)
    lengths.sort(reverse=True)
    half, run, n50, l50 = total / 2, 0, 0, 0
    for i, L in enumerate(lengths, 1):
        run += L
        if run >= half:
            n50, l50 = L, i
            break
    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "longest_bp": lengths[0] if lengths else 0,
        "shortest_bp": lengths[-1] if lengths else 0,
        "mean_bp": round(total / len(lengths), 1) if lengths else 0,
        "n50_bp": n50,
        "l50": l50,
        "gc_percent": round(100 * gc / max(1, total - n_amb), 2),
        "ambiguous_bases": n_amb,
        "contigs_ge_1kb": sum(1 for L in lengths if L >= 1000),
        "contigs_ge_5kb": sum(1 for L in lengths if L >= 5000),
    }


def parse_header(h):
    """ABRicate format: DB~~~GENE~~~ACCESSION~~~PRODUCT"""
    parts = h.split("~~~")
    if len(parts) >= 4:
        return parts[0], parts[1], parts[2], parts[3].strip()
    return "?", h.split()[0], "", h


def is_protein(path, n=20):
    """crude alphabet check on the first n records"""
    nt = set("ACGTUNRYKMSWBDHV")
    seen, non_nt = 0, 0
    for _, seq in read_fasta(path):
        s = set(seq.upper())
        if s and len(s - nt) / len(s) > 0.3:
            non_nt += 1
        seen += 1
        if seen >= n:
            break
    return non_nt > seen / 2


def overlap(a, b):
    lo = max(a[0], b[0])
    hi = min(a[1], b[1])
    return max(0, hi - lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--dbdir", required=True)
    ap.add_argument("--dbs", nargs="+", required=True)
    ap.add_argument("--minid", type=float, default=80.0)
    ap.add_argument("--mincov", type=float, default=60.0)
    ap.add_argument("--preset", default="asm20")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)

    stats = assembly_stats(a.contigs)
    print("=== assembly ===")
    for k, v in stats.items():
        print(f"{k:22s} {v}")
    with open(os.path.join(a.out, "assembly_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)

    print("\nindexing contigs ...", flush=True)
    aln = mappy.Aligner(a.contigs, preset=a.preset, best_n=50)
    if not aln:
        sys.exit("failed to index contigs")

    all_rows = []
    per_db_counts = {}
    for db in a.dbs:
        path = os.path.join(a.dbdir, db, "sequences")
        if not os.path.exists(path):
            print(f"  ! {db}: no sequences file, skipped")
            continue
        if is_protein(path):
            print(f"  ! {db}: protein database, skipped "
                  f"(nucleotide alignment only)")
            per_db_counts[db] = "skipped (protein db)"
            continue

        rows, n_ref = [], 0
        for hdr, seq in read_fasta(path):
            n_ref += 1
            if len(seq) < 90:
                continue
            dbn, gene, acc, product = parse_header(hdr)
            for h in aln.map(seq):
                ident = 100.0 * h.mlen / max(1, h.blen)
                cov = 100.0 * (h.q_en - h.q_st) / max(1, len(seq))
                if ident < a.minid or cov < a.mincov:
                    continue
                rows.append({
                    "database": db, "gene": gene, "accession": acc,
                    "product": product, "contig": h.ctg,
                    "c_start": h.r_st, "c_end": h.r_en,
                    "strand": "+" if h.strand > 0 else "-",
                    "ref_len": len(seq),
                    "identity": round(ident, 2), "coverage": round(cov, 2),
                    "score": round(ident * cov / 100.0, 2),
                })

        # resolve overlapping hits within a database: greedy best-first
        kept, by_ctg = [], defaultdict(list)
        for r in sorted(rows, key=lambda x: -x["score"]):
            iv = (r["c_start"], r["c_end"])
            if any(overlap(iv, k) > 0.6 * min(iv[1] - iv[0], k[1] - k[0])
                   for k in by_ctg[r["contig"]]):
                continue
            by_ctg[r["contig"]].append(iv)
            kept.append(r)

        per_db_counts[db] = len(kept)
        all_rows.extend(kept)
        print(f"  {db:14s} {n_ref:6d} refs -> {len(rows):5d} raw, "
              f"{len(kept):4d} kept")

    cols = ["database", "gene", "accession", "contig", "c_start", "c_end",
            "strand", "ref_len", "identity", "coverage", "score", "product"]
    out_tsv = os.path.join(a.out, "arg_hits.tsv")
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(all_rows, key=lambda x: (x["database"], -x["score"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\ntotal hits: {len(all_rows)}  ->  {out_tsv}")
    with open(os.path.join(a.out, "per_db_counts.json"), "w") as fh:
        json.dump(per_db_counts, fh, indent=2)


if __name__ == "__main__":
    main()
