# Validation and reproducibility

`selftest` exercises parsing, source duplicates, candidate conflicts, marker naming,
exact placement, CIGAR/allele checks, cache integrity, cohort immutability, query
bounds, nearest distances, deterministic subsampling/backgrounds, dosage r2,
shared collections, workbook export, location safety and numerical comparisons.
Fixtures in tests are explicitly synthetic software tests, not research data.
A raw-input miniature orchestration test substitutes known native placements and
requires no pre-existing result folders. This test excludes the production
native binary and the whole cacao reference. Native executables separately pass a
seven-case software smoke test before processing cacao tags; its verified cache
key includes executable hashes and settings.

The development acceptance comparison re-parsed eight actual source workbooks and
recomputed downstream results using preserved reference/catalog/mapping decisions
as fixtures when whole reference/native execution was unavailable in that runtime.
This establishes the integrated numerical and source-reader interface behavior,
not an independent repetition of whole-genome alignment. Full numerical
reproduction is evaluated by a raw-input run followed by comparison with the
designated reference outputs on the target Windows environment.

The optional `verify-results` operation compares every non-README worksheet of six
coverage workbooks and five dosage workbooks. Record keys, sample membership,
counts, distance values, exclusions, iteration identities and strings must agree;
finite floating values have an absolute tolerance of 1e-12. The `engine` reader
provenance column is the sole ignored column. A short explicit compatibility map
renames old source-eligibility headers and one exact-policy comparison phrase; it
does not change biological labels or numbers. Metadata and workbook styling are
not used as proxies for numerical equality. All rows are compared, not a sample.

Comparison references are optional, externally supplied and hash checked. They
are used to compare completed outputs, rather than as analysis inputs. A missing
workbook, changed reference, column mismatch or altered value is a failed comparison.

Interpretive limitations remain: alternative alignments
outside bounded search are not excluded; Nacional sequence coordinates are linked
through a shared resource; genotype coding remains conditional; Amazonia has no
approved encoding/sequence bridge here; no causal effect, flavor phenotype,
haplotype, optimal window, ancestry-controlled tagging or breeding outcome is tested.
