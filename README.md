# package-ApplyPolygenicScore

An R package that provides utilities for the application of an existing polygenic score (PGS) to genotype data in VCF form.
Matches coordinates of a provided polygenic score weight file to VCF a input and calculates a weighted sum of dosages in each individual.

## Description
This tool is intended to simply and transparently parse genotype/dosage data from an input VCF, match genotype coordinates to the component SNPs of an existing polygenic score, and apply SNP weights to dosages to calculate a polygenic score for each individual in accordance with the additive weighted sum of dosages model.

## Installation
To install the latest development version from GitHub:

```
# install.packages("devtools")

devtools::install_github("uclahs-cds/package-ApplyPolygenicScore")
```

## Getting Started
This package is a fantastic resource for people in the following scenario: you have just received some genetic data from a DNA sequencing or genotyping experiment, or perhaps you have just downloaded one such dataset from a public repository. You noticed that a research group or genetic consortium has recently published a new polygenic score for the trait you are interested in studying. You would now like to apply that PGS to your data. ApplyPolygenicScore provides all the functions you need to import all required input data into R, perform the necessary calculations, and transition smoothly into analysis.

Below is an overview of the data you will require to get started. For even more details, check out our discussion on [What is a PGS?](https://github.com/uclahs-cds/package-ApplyPolygenicScore/discussions/2).

### Input Data
You will need only two pieces of data to get started: 
1. A VCF file: Genotype data of the individuals upon which you wish to apply a polygenic score.
2. A PGS weight file: Coordinates of each SNP that compose the polygenic score you with to apply, and their associated weights.

#### Genotype data
Genotype data should be provided in the form of a VCF file. Others have done a great job of describing Variant Call Format. For those with a basic understanding of genetic nomenclature, we recommend the GATK [resource](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format). For those who need a refresher on genomics and genomic data, we recommend starting with the [fact sheets](https://www.genome.gov/about-genomics/fact-sheets) curated by the National Human Genome Research Institute (NHGRI).

If you wish to apply a PGS to a cohort, we recommend that genotypes for the whole cohort be aggregated in one VCF file, either through a regenotyping process, or through VCF merging with an external tool designed for manipulating VCF files. VCF files can be very large, causing memory-related complications in the R environment. To reduce memory usage and improve speed of PGS application, we recommend pre-filtering the input VCF for only the coordinates that compose the PGS you wish to apply. This action can be performed using a coordinate BED file and tools such as bcftools or bedtools. To facilitate this process, ApplyPolygenicScore provides a function that outputs a BED file containing coordinates for any number of PGS weight files provided as input.

#### PGS weight file
The PGS weight file describes a PGS by providing a list of component SNPs, their genomic coordinates, and their respective weights. [The PGS Catalog](https://www.pgscatalog.org/) is a public database of PGSs and their weight files, and a great first stop for acquiring a PGS weight file. The functions of ApplyPolygenicScore have been designed to operate on weight files that have been formatted according to the standards established by The PGS Catalog. There are very well documented [here](https://www.pgscatalog.org/downloads/#:~:text=the%20different%20samples.-,PGS%20Scoring%20Files,-Formatted%20Files). However, you could easily create your own compatible PGS weight file, simply by formatting all required columns by Catalog standards. When in doubt, use our `check.pgs.weight.columns()` function to make sure any data table you import into R contains the reqiuired columns for downstream functions.

### Recommended Workflow

## Resources

## Getting Help

## Contributors

## License

Author: Nicole Zeltser(nzeltser@mednet.ucla.edu)

package-ApplyPolygenicScore is licensed under the GNU General Public License version 2. See the file LICENSE.md for the terms of the GNU GPL license.

A package that provides utilities for the application of an existing polygenic score to a VCF.

Copyright (C) 2021 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
