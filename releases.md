## Changelog for releases
### 08.09.2021: v0.5.0
Name change from anonymizeBAM to BAMboozle.
### 11.01.2021: v0.4.5
Add correct installation dependencies.
### 11.01.2021: v0.4.4
Fix version verbose.
### 11.01.2021: v0.4.3
Added `@PG` header line to mark output .bam files with the anonymizeBAM call.
### 07.01.2021: v0.4.2
Changes for upload to PyPI
### 06.01.2021: v0.4.1
Added a switch to keep/discard secondary alignments.
Correctly deal with unmapped reads that inhert the chromosome & mapping position of their mate.
Sanitize bwa-mem auxiliary tags with information on mismatches.
Speed up indexing by passing threads argument.
### 05.01.2021: v0.4
Initial release to GitHub
