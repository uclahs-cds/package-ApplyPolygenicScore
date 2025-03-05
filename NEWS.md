# Unreleased

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
