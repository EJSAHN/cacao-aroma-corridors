# Validation and reproducibility

`selftest` exercises parsing, source duplicates, candidate conflicts, marker naming,
exact placement, CIGAR/allele checks, cache integrity, cohort immutability, query
bounds, nearest distances, deterministic subsampling/backgrounds, dosage r2,
shared collections, workbook export, location safety and numerical comparisons.
Fixtures in tests are explicitly synthetic software tests, not research data.
A raw-input miniature orchestration test substitutes known native placements and
requires no pre-existing result folders. It does not claim to test the production
native binary or the whole cacao reference. Native executables separately pass a
seven-case software smoke test before processing cacao tags; its verified cache
key includes executable hashes and settings.

The development acceptance comparison re-parsed eight actual source workbooks and
recomputed downstream results using preserved reference/catalog/mapping decisions
as fixtures when whole reference/native execution was unavailable in that runtime.
This establishes the integrated numerical and source-reader interface behavior,
not an independent repetition of whole-genome alignment. The local raw-input run
and external reference comparison provide the final execution check on the target
Windows environment. Only a passing local agreement report should support a claim
that the complete public-source run reproduces the designated reference outputs.

The optional `verify-results` operation compares every non-README worksheet of six
coverage workbooks and five dosage workbooks. Record keys, sample membership,
counts, distance values, exclusions, iteration identities and strings must agree;
finite floating values have an absolute tolerance of 1e-12. The `engine` reader
provenance column is the sole ignored column. A short explicit compatibility map
renames old source-eligibility headers and one exact-policy comparison phrase; it
does not change biological labels or numbers. Metadata and workbook styling are
not used as proxies for numerical equality. All rows are compared, not a sample.

Comparison references are optional, externally supplied and hash checked. They
are never used to fill scientific outputs. A missing workbook, changed reference,
column mismatch or altered value is a failed comparison. No reference fixture is
included as a hidden analysis input in the public source.

Known limitations are substantive, not software failures: alternative alignments
outside bounded search are not excluded; Nacional sequence coordinates are linked
through a shared resource; genotype coding remains conditional; Amazonia has no
approved encoding/sequence bridge here; no causal effect, flavor phenotype,
haplotype, optimal window, ancestry-controlled tagging or breeding outcome is tested.
