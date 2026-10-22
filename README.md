# Resistance determinants in industrial effluent metagenomes

A reproducible screen of two publicly deposited metagenome assemblies for antimicrobial
resistance genes, biocide resistance determinants, class 1 integron markers and plasmid
replicons.

Both assemblies are from pharmaceutical and industrial effluent treatment in India and were
generated during the work described in Jadeja *et al.* (2014) and Jadeja *et al.* (2019).
This repository re-analyses them with current reference databases and reports what is, and
is not, supportable from assemblies of this fragmentation.

---

## The assemblies

| Accession | Platform | Contigs | Total | N50 | Longest | GC |
|---|---|---:|---:|---:|---:|---:|
| **AERA01** | 454 pyrosequencing | 36,270 | 27.9 Mb | 849 bp | 128.8 kb | 58.9% |
| **JRYH01** | Illumina | 159,908 | 130.8 Mb | 1,503 bp | 398.6 kb | 62.9% |

Both are retrievable from NCBI under their WGS accessions.

---

## Method

`scripts/argscan.py` aligns each contig against reference gene sets with
[`mappy`](https://pypi.org/project/mappy/) (minimap2, `asm20` preset, `best_n=50`), then:

- computes identity as `100 × mlen / blen` and coverage as the fraction of the reference
  gene spanned by the alignment;
- resolves overlapping hits greedily at a 0.6 reciprocal-overlap threshold, so a gene called
  by several databases counts once;
- reports every hit (`*_hits.tsv`) and the resolved non-redundant set (`*_loci.tsv`).

`scripts/context.py` clusters hits into loci by contig position, canonicalises gene names
across databases, and reports which contigs carry more than one distinct determinant.

**Databases.** CARD, NCBI AMRFinder, ResFinder, ARG-ANNOT, MEGARes, BacMet2 (biocide and
metal resistance) and PlasmidFinder, as distributed with ABRicate.

**Thresholds.** 80% identity, 60% coverage for the reported set. A relaxed run at lower
thresholds was used only to test whether sparse results were a threshold artefact.

```bash
python scripts/argscan.py \
  --contigs AERA01.1.fsa_nt.gz \
  --dbdir  abricate/db \
  --dbs    card,ncbi,resfinder,argannot,megares,bacmet2,plasmidfinder \
  --minid  80 --mincov 60 \
  --out    results/AERA01_hits.tsv

python scripts/context.py \
  --hits    results/AERA01_hits.tsv \
  --contigs AERA01.1.fsa_nt.gz \
  --out     results/AERA01_loci.tsv
```

---

## Results

**AERA01** — 134 raw hits resolving to **34 distinct loci**.
**JRYH01** — 15 raw hits resolving to **5 distinct loci**.

Selected loci at full identity and coverage:

| Gene | Contig | Contig length | Identity | Coverage | Databases |
|---|---|---:|---:|---:|---:|
| **blaVIM** | AERA01023697.1 | 1,505 bp | 100.0% | 100.0% | 5 |
| **aac(6′)** | AERA01023697.1 | 1,505 bp | 100.0% | 100.0% | 4 |
| **dfrA1** | AERA01035294.1 | 581 bp | 100.0% | 100.0% | 5 |
| **dfrB5** | AERA01018650.1 | 411 bp | 100.0% | 100.0% | 5 |
| **blaOXA** | AERA01018647.1 | 831 bp | 100.0% | 100.0% | 5 |
| **blaGES** | AERA01001083.1 | 1,020 bp | 100.0% | 100.0% | 5 |
| **sul1** | AERA01018605.1 | 810 bp | 100.0% | 97.8% | 5 |
| **sul2** | AERA01009517.1 | 1,347 bp | 100.0% | 100.0% | 5 |
| **strA / strB** | AERA01018606.1 | 2,440 bp | 100.0% | 100.0% | 5 |
| **msr(E) / mph** | AERA01017201.1 | 2,781 bp | 100.0% | 100.0% | 5 |
| **arr** (rifampicin) | AERA01014631.1 | 609 bp | 100.0% | 100.0% | 5 |
| **IncP-1 replicon** | AERA01000503.1 | 2,007 bp | 99.1% | 100.0% | 1 |

The resolved set spans β-lactamases (VIM, GES, OXA, a class A enzyme), aminoglycoside
modifying enzymes, sulphonamide and trimethoprim resistance, macrolide, phenicol,
tetracycline and rifampicin determinants, biocide efflux (*qacE*, four loci), RND efflux,
a tigecycline efflux component (*toprJ1*) and a broad-host-range IncP-1 plasmid replicon.

### The one genuine co-location

**blaVIM and aac(6′) are adjacent on AERA01023697.1**, both at 100% identity and coverage,
within a 1,505 bp contig — a carbapenemase and an aminoglycoside acetyltransferase in
tandem, the arrangement expected of an integron cassette array.

The other multi-locus contigs carry gene pairs that are near-universally linked anyway
(*strA–strB*, *msr(E)–mph*, *toprJ1–*RND efflux) and so are not independent evidence of
cassette structure.

---

## Limitations

These are stated plainly because they bound what the data can be used to claim.

**Contig length.** Most loci sit on contigs of 400–1,500 bp — barely longer than the genes
themselves. **No claim of physical linkage between genes on separate contigs is supportable
from these assemblies.** Determining whether *sul1*, *qacEΔ1*, *dfrA1* and *blaVIM* share a
genetic element requires long reads.

**Integron structure is inferred, not observed.** `intI1` was not detected at either strict
or relaxed thresholds. The class 1 integron is inferred from the presence of *sul1* and
*qacEΔ1*, which constitute its 3′ conserved segment — so their co-occurrence is a property
of the element's definition and is not independent evidence of biocide co-selection.

**No metal resistance determinants were found.** *mer*, *ars*, *czc*, *pco*, *sil*, *cop*,
*chr*, *znt* and *cad* loci were absent from both assemblies in all databases, BacMet2
included.

**Assembly asymmetry.** The larger Illumina assembly (JRYH01, 131 Mb) yields far fewer loci
than the smaller 454 assembly (AERA01, 28 Mb). A self-mapping control and a relaxed-threshold
re-run confirmed this is not a detection artefact; the two samples differ.

---

## Files

```
scripts/argscan.py        contig-vs-database alignment and hit resolution
scripts/context.py        locus clustering, name canonicalisation, co-location
results/AERA01_hits.tsv   all alignments passing threshold
results/AERA01_loci.tsv   resolved non-redundant loci
results/JRYH01_hits.tsv
results/JRYH01_loci.tsv
```

Column definitions are in the header of each file.

---

## References

Jadeja NB, More RP, Purohit HJ, Kapley A (2014). Metagenomic analysis of oxygenases from
activated sludge. *Bioresource Technology* 165: 250–256.

Jadeja NB, Purohit HJ, Kapley A (2019). Decoding microbial community intelligence through
metagenomics for efficient wastewater treatment. *Functional & Integrative Genomics* 19:
839–851.

---

## License

MIT. Reference databases retain their own licences.
