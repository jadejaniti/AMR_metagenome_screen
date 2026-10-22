#!/usr/bin/env python3
"""Non-redundant ARG set and contig-level co-location analysis."""
import argparse
import csv
import gzip
from collections import defaultdict

# map database-specific names onto a common label
ALIASES = [
    ("blavim", "blaVIM"), ("vim-", "blaVIM"),
    ("blages", "blaGES"), ("ges-", "blaGES"),
    ("blaoxa", "blaOXA"), ("oxa-", "blaOXA"),
    ("sul1", "sul1"), ("sul2", "sul2"), ("sul3", "sul3"),
    ("arr", "arr (rifampicin)"),
    ("dfra", "dfrA"), ("dfrb", "dfrB"),
    ("stra", "strA/aph(3'')"), ("aph3-dprime", "strA/aph(3'')"),
    ("strb", "strB/aph(6)"), ("aph6", "strB/aph(6)"),
    ("aac6", "aac(6')"), ("aac(6", "aac(6')"),
    ("aac3", "aac(3)"), ("aac(3", "aac(3)"),
    ("ant", "ant (aminoglycoside)"), ("aad", "ant (aminoglycoside)"),
    ("msr", "msr(E) macrolide"), ("mph", "mph macrolide"),
    ("cmla", "cmlA chloramphenicol"), ("flor", "floR"),
    ("cat", "cat chloramphenicol"),
    ("tet", "tet (tetracycline)"),
    ("qac", "qacE (biocide)"),
    ("inti", "intI1 integrase"), ("int1", "intI1 integrase"),
    ("ere", "ere macrolide"), ("erm", "erm macrolide"),
    ("mex", "mex efflux"), ("acr", "acr efflux"), ("oqx", "oqx efflux"),
    ("blaa", "class A beta-lactamase"),
    ("blaz", "blaZ"), ("ampc", "ampC"), ("blatem", "blaTEM"),
    ("blactx", "blaCTX-M"), ("blakpc", "blaKPC"), ("blandm", "blaNDM"),
]


def canon(gene):
    g = gene.lower().lstrip("(").replace(" ", "")
    for pre in ("agly)", "bla)", "mls)", "phe)", "rif)", "sul)", "tmt)",
                "tet)", "col)", "fcyn)", "flq)", "gly)", "ntmdf)", "tric)"):
        if g.startswith(pre):
            g = g[len(pre):]
    for key, label in ALIASES:
        if g.startswith(key) or key in g:
            return label
    return gene


def read_contig_lengths(path):
    lens, name = {}, None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lens[name] = 0
            elif name:
                lens[name] += len(line.strip())
    return lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.hits), delimiter="\t"))
    for r in rows:
        r["canon"] = canon(r["gene"])
        r["identity"] = float(r["identity"])
        r["coverage"] = float(r["coverage"])
        r["c_start"], r["c_end"] = int(r["c_start"]), int(r["c_end"])

    clen = read_contig_lengths(a.contigs)

    # --- non-redundant gene set: cluster by contig + overlapping locus
    loci = defaultdict(list)
    for r in rows:
        loci[r["contig"]].append(r)

    clusters = []
    for ctg, rs in loci.items():
        rs.sort(key=lambda x: x["c_start"])
        cur = []
        for r in rs:
            if cur and r["c_start"] < max(x["c_end"] for x in cur) - 30:
                cur.append(r)
            else:
                if cur:
                    clusters.append((ctg, cur))
                cur = [r]
        if cur:
            clusters.append((ctg, cur))

    print(f"{len(rows)} hits -> {len(clusters)} distinct loci on "
          f"{len(loci)} contigs\n")

    # consensus call per locus: most frequent canonical name, best identity
    summary = []
    for ctg, rs in clusters:
        names = defaultdict(int)
        for r in rs:
            names[r["canon"]] += 1
        best = max(names.items(), key=lambda kv: kv[1])[0]
        top = max(rs, key=lambda x: x["identity"] * x["coverage"])
        summary.append({
            "contig": ctg, "contig_len": clen.get(ctg, 0),
            "start": min(r["c_start"] for r in rs),
            "end": max(r["c_end"] for r in rs),
            "gene": best,
            "databases": ",".join(sorted({r["database"] for r in rs})),
            "n_db": len({r["database"] for r in rs}),
            "best_identity": top["identity"],
            "best_coverage": top["coverage"],
            "raw_names": ";".join(sorted({r["gene"] for r in rs})),
            "product": top["product"][:90],
        })

    summary.sort(key=lambda x: (-x["n_db"], -x["best_identity"]))
    cols = ["contig", "contig_len", "start", "end", "gene", "databases",
            "n_db", "best_identity", "best_coverage", "raw_names", "product"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for s in summary:
            w.writerow(s)

    print("=== non-redundant loci (confirmed by >=2 databases) ===")
    print(f"{'gene':26s} {'contig':22s} {'len':>7s} {'id%':>6s} "
          f"{'cov%':>6s} {'dbs':>4s}")
    for s in summary:
        if s["n_db"] >= 2:
            print(f"{s['gene']:26s} {s['contig']:22s} {s['contig_len']:7d} "
                  f"{s['best_identity']:6.1f} {s['best_coverage']:6.1f} "
                  f"{s['n_db']:4d}")

    single = [s for s in summary if s["n_db"] == 1]
    if single:
        print(f"\n=== single-database loci ({len(single)}) ===")
        for s in single:
            print(f"{s['gene']:26s} {s['contig']:22s} {s['contig_len']:7d} "
                  f"{s['best_identity']:6.1f} {s['best_coverage']:6.1f}  "
                  f"[{s['databases']}]")

    # --- co-location
    multi = defaultdict(list)
    for s in summary:
        multi[s["contig"]].append(s)
    multi = {k: v for k, v in multi.items() if len(v) > 1}
    print(f"\n=== contigs carrying >1 distinct resistance locus "
          f"({len(multi)}) ===")
    for ctg, ss in sorted(multi.items(), key=lambda kv: -kv[1][0]["contig_len"]):
        ss.sort(key=lambda x: x["start"])
        print(f"\n{ctg}  ({ss[0]['contig_len']} bp)")
        for s in ss:
            print(f"   {s['start']:7d}-{s['end']:<7d} {s['gene']:26s} "
                  f"{s['best_identity']:5.1f}% id  {s['best_coverage']:5.1f}% cov")


if __name__ == "__main__":
    main()
