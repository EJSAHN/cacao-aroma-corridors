# Analysis definitions

This document describes implemented operations. It does not establish phenotype
effects, favorable alleles, haplotype conservation or an optimal breeding window.
All source-specific names, chromosomes, windows, sample sensitivity IDs and QC
thresholds are explicit in `config/study.json`; no private drive paths are embedded.

## Source units and identity

Association worksheets are retained separately. Only the declared master is the
primary association-record source. Unique event keys distinguish chromosome, peak,
trait, reported model and marker filter. A distinct chromosome/peak pair is a
physical point, not an independent association discovery. The fruity subset is a
separately named point set. Source-defined intervals are used only when boundaries
are supplied, ordered, in range and compatible with the peak. Missing boundaries
are not replaced with an assumed LD block. Potential rotated fields and trait-name
aliases are reported but do not silently edit primary records. Qualitative candidate
flags (`Near`, `?`) remain distinct from explicit membership. Repeated gene IDs with
consistent coordinates are consolidated; conflicting intervals are quarantined.

Raw and narrowly canonicalized marker identifiers are compared independently from
coordinates. Canonicalization recognizes configured DArT naming conventions, not
arbitrary punctuation removal or bare RefSeq-version digits. Identifier collisions
and coordinate disagreements are documented. Gene-body recovery uses a 1-based
inclusive interval; peak exact recovery requires the identical position on the same
chromosome. A flank of 50,000 means 50 kb on **each** side, not a 50-kb total width.

## Reference placement

Both allelic tags are reconstructed from five-prime flank, directional base
substitution and three-prime flank. A reference-sequence key is supported only when
source tag identity is established. The exact policy searches all reference contigs
in both orientations using an SNP-free seed and full-tag equality. Its uniqueness
claim is limited to exact full-tag matches.

The bounded alignment policy searches all eligible tags, including exact-only
successes, using Bowtie 2 2.5.5, end-to-end, very-sensitive, FASTA, quality-ignored,
seed 1729, four threads, `-k 10`, `-L 18`, `-N 1`, `-D 25`, `-R 4`, interval
`S,1,0.5`, constant minimum score -30, mismatch penalty 6, gap open/extend 5/3,
N penalty 1 and gap barrier 4. Candidate reporting is heuristic and bounded, not
an exhaustive proof of global mapping uniqueness. Reaching the reported-hit cap
causes conservative exclusion. The two allelic tags must support the same focal
position and direction. CIGAR operations locate the focal base; soft clips,
nonfocal edits above 3, gap bases above 2, a gap at/within the 3-base focal guard,
and insufficient best-to-alternative score separation (6) are excluded. The parser
also checks reference bases and score-floor censoring. Full decision and top-two
allelic alignment evidence are preserved; all SAM records remain in the cache.

Diversity is mapped from its own tags. Nacional coordinates are projected only
through a unique supported reference-sequence key. This projection is **not an
independent sequence validation** of Nacional or of an aroma-associated allele.
Original association peaks are not moved through positional coincidence alone.
All candidate gene IDs and boundaries are checked against the pinned gene catalog.

## Fixed A/B/C cohorts and coverage

A contains source-eligible primary-chromosome markers whose original position is
in range. B is the retained mapping-policy subset with the same original positions.
C consists of the exact same marker records as B with mapped SNP positions. A
record that cannot enter B cannot appear only in C. `B-A` is subset selection;
`C-B` isolates coordinate change on fixed records; `C-A` combines both. Position
counts for coverage use unique coordinates. Exact-only versus bounded-alignment
comparisons are separate mapping-policy sensitivity analyses.

Distances are nearest-target-marker distances on the same chromosome. For an
interval, distance is zero for a gene-body hit and otherwise is the distance to
the nearest boundary. Missing target chromosomes are misses, not zero distances.
All configured windows are evaluated without optimizing a cutoff. Query sets are
all unique candidates, explicit pathway memberships, unique master peaks, fruity
subset peaks, valid source intervals and a coding-only candidate sensitivity set.
Overlap among pathways does not create independent gene counts.

## Marker counts and annotated-gene backgrounds

Chromosome-specific minimum unique marker counts across the comparison resources
set the quotas. Sampling without replacement uses 200 deterministic replicates.
A and C have separate quotas/draws and do not constitute a paired correction
experiment; use fixed B/C for that interpretation. A resource that determines all
minimum counts may have no variability across count-matched replicates.

All original candidate IDs are excluded from the background catalog. Noncandidate
gene sets are matched on chromosome, biotype and gene-length bins (five bins),
optionally adding local annotated-gene-density bins (three bins, 250-kb radius).
Background genes are sampled without replacement, 1,000 draws per evaluable model.
The same draw is applied across resource/coordinate arms. No depleted stratum is
silently widened; incomplete models are flagged without dropping target genes.
The full candidate and protein-coding-only sets are distinct analyses. Background
quantiles describe repeated draws, not biological confidence intervals. No p-value,
FDR or claim of enrichment against the original tested-GWAS universe is produced.
Count matching and gene-background matching are separate diagnostic controls, not
joint adjustment for all ascertainment, relatedness or spatial dependence.

## Conditional dosage correlations

The declared Nacional coding is A/M/C = 0/1/2, conditional on export symbols. It
is not biological reference/alternate or favorable-allele orientation. No genotype
is imputed, changed or swapped. Unsupported symbols/declarations and duplicate
corrected positions are excluded explicitly. Sample and marker call rates are at
least 0.90, MAF at least 0.05, and each pair requires at least 50 complete samples.
A fixed marker set passes QC in both inclusive and prespecified sample-omission
conditions. `SNA604` and `L17H53` are omitted only in the sensitivity group because
between-export call differences concentrate there, not because error is established.

Pairwise-complete Pearson dosage correlations squared are computed in bounded
batches. All pairs within 1 Mb under either original or corrected coordinates are
retained. The same genotype r2 is assigned to each coordinate representation. Low
and zero values are not filtered. Separately sampled same-chromosome pairs >=5 Mb
and interchromosomal pairs (up to 5,000 each) provide descriptive context. Shared
G7/MAF5 comparisons fix the identical marker keys, corrected positions, samples and
pair inventory. Per-SNP neighbor summaries exclude self and refer only to other
observed SNPs. They do not test tagging of an absent causal variant.

These are unphased, unadjusted dosage correlations, not haplotype-frequency r2,
structure-corrected LD, genome-wide association validation or breeding predictions.
