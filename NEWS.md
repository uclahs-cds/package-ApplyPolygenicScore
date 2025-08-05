# Unreleased

# ApplyPolygenicScore unreleased

# ApplyPolygenicScore 4.0.0

## Changed
* Refactored all main functions for gains in RAM efficiency and runtime
* Introduced a breaking change to the output of `import.vcf`. The outputed list object has a different naming scheme and different data formats. Previous data formats are still supported by setting `long.format` to `TRUE`, however the naming scheme is still changed.
* Introduced a breaking change to `apply.polygenic.score`. The expected default `vcf.data` input format has changed. The previous input format is still supported by setting `vcf.long.format` to `TRUE` from the default `FALSE`.

## Added
* Added support for more efficient storage and manipulation of imported VCF data. The default output of `import.vcf` now returns VCF data in a split format. A `data.table` object contains VCF data from fixed fields (CHROM, POS, ID, REF, ALT). A `matrix` object contains sample-specific genotypes in allele-format in a sample (columns) by variant (rows) matrix.

# ApplyPolygenicScore 3.1.0

## Changed
* Fixed regression of combine.vcf.with.pgs() function that prevented it from handling multiple rsIDs on the same line.
* Fixed bug caused by the case of a sample-specific missing variant at a multiallelic site
* Use updated `R CMD check` CI/CD action with `renv` dependency management

## Added
* Added new contributor
* Added minimum sample size check for grouped density curves
* Added new plotting function `create.pgs.boxplot`
* Added option for user to provide custom PGS source column(s) for plotting functions
* Added option to `assess.pgs.vcf.allele.match` to condition the handling of ambiguous strand flips on the total number of unambiguous strand flips.
* Added new function `analyze.pgs.binary.predictiveness` which given a PGS and phenotypes runs a logistic regression and returns statistics (OR, p-value, AUC) describing how well the PGS predicts the phenotype. It also automatically plots a receiver-operator-curve.

# ApplyPolygenicScore 3.0.2

## Changed
* ApplyPolygenicScore released on CRAN! Updated README with CRAN links.

# ApplyPolygenicScore 3.0.1

## Added
* Added hemizygous allele handling to dosage calculation
* Added toggle to hexbinplot at sample size threshold in `create.pgs.with.continuous.phenotype.plot`

## Changed
* Updated INDEL effect switch reporting by strand flip checker
* Updated data structuring for automated statistical analysis in `apply.polygenic.score`

# ApplyPolygenicScore 3.0.0 (2024-12-02)

## Added
* Added handling of overlapping deletion allele notation
* Added secondary PGS/VCF variant matching method using rsID after first attempt with genomic coordinates
* Added checks for rsID as an optional column in input PGS weight files
* Added functionality to assess allele matches and correct strand flips

# ApplyPolygenicScore 2.0.0 (2024-07-31)

## Changed
* Renamed functions starting with reserved vocabulary for S3 generic methods merge. -> combine.

# ApplyPolygenicScore 1.0.0 (2024-07-25)

* First release

# ApplyPolygenicScore 0.1.0 (2024-07-25)

* INITIAL FEATURES
