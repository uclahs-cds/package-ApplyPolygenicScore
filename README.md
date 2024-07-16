# package-ApplyPolygenicScore

1. [Description](#description)
2. [Installation](#installation)
3. [Getting Started](#getting-started)
4. [Resources](#resources)
5. [Contributors](#contributors)
6. [License](#license)


## Description
This R package provides a set of utilities to simply and transparently parse genotype/dosage data from an input VCF, match genotype coordinates to the component SNPs of an existing polygenic score, and apply SNP weights to dosages to calculate a polygenic score for each individual in accordance with the additive weighted sum of dosages model.

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
- Genotype data should be provided in the form of a VCF (Variant Call Format) file.
- Others have done a great job of describing Variant Call Format. For those with a basic understanding of genetic nomenclature, we recommend the GATK [resource](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format).
- For those who need a refresher on genomics and genomic data, we recommend starting with the [fact sheets](https://www.genome.gov/about-genomics/fact-sheets) curated by the National Human Genome Research Institute (NHGRI).

If you wish to apply a PGS to a cohort, we recommend that genotypes for the whole cohort be aggregated in one VCF file, either through a regenotyping process, or through VCF merging with an external tool designed for manipulating VCF files. VCF files can be very large, causing memory-related complications in the R environment. To reduce memory usage and improve speed of PGS application, we recommend pre-filtering the input VCF for only the coordinates that compose the PGS you wish to apply. This action can be performed using a coordinate BED file and tools such as bcftools or bedtools. To facilitate this process, ApplyPolygenicScore provides a function that outputs a BED file containing coordinates for any number of PGS weight files provided as input.

#### PGS weight file
- The PGS weight file describes a PGS by providing a list of component SNPs, their genomic coordinates, and their respective weights.
- [The PGS Catalog](https://www.pgscatalog.org/) is a public database of PGSs and their weight files, and a great first stop for acquiring a PGS weight file.
- The functions of ApplyPolygenicScore have been designed to operate on weight files that have been formatted according to the standards established by The PGS Catalog. These are very well documented [here](https://www.pgscatalog.org/downloads/#:~:text=the%20different%20samples.-,PGS%20Scoring%20Files,-Formatted%20Files).
- You could easily create your own compatible PGS weight file, simply by formatting all required columns by Catalog standards.
- When in doubt, use our `check.pgs.weight.columns()` function to make sure any data table you import into R contains the reqiuired columns for downstream functions.

### Recommended Workflow


1. Convert PGS weight files to BED coordinate files

    We recommend starting by filtering your input VCF for just the variants in your PGS weight files. Several software tools are available to do this, and most all require a coordinate BED file. A description of BED format can be found [here](https://bedtools.readthedocs.io/en/latest/content/general-usage.html).

    The function `import.pgs.weight.file` can be used to import your PGS weight files into R.
    The functions `convert.pgs.to.bed` and `merge.pgs.bed` can be used to make the conversion, and merge several BED dataframes into one, respectively.

2. Import your VCF file.

    Once you have filtered down your VCF, simply import it into R using `import.vcf`. This function is a wrapper of `vcfR::vcfR2tidy` that ensures all required fields are imported.


3. Apply your PGS.

    Provide your imported VCF and PGS weight files to `apply.polygenic.score`. It's as simple as that.
    Under the hood, this function begins by calling `merge.vcf.with.pgs`. The merge function also outputs a list of variants in your PGS that could not be found in your VCf data, which you can obtain by calling the function independently.
    `apply.polygenic.score` outputs lots of useful information along with the score and provides various customizeable options, such as methods for handling missing sites (see [this discussion](https://github.com/uclahs-cds/package-ApplyPolygenicScore/discussions/17) for more) and basic analyses with phenotype data.

4. Create summary plots.

    ApplyPolygenicScore comes with several plotting functions designed to operate on the results of `apply.polygenic.score`. Display PGS density curves with `create.pgs.density.plot` and PGS percentile ranks with `create.pgs.rank.plot`. If you provided phenotype data in step 3, you can incorporate categorical data into the density plots and categorical and continuous phenotype data into the rank plots, and use `create.pgs.with.continuous.phenotype.plot` to make scatterplots of your PGS against any continuous phenotype data.

For more step-by-step instructions, check out our vignettes.

## Resources

## Getting Help
Looking for guidance or support with ApplyPolygenicScore? Check out our [Discussions](https://github.com/uclahs-cds/package-ApplyPolygenicScore/discussions) page.

Submit bugs, suggest new features or see current works at our [Issues](https://github.com/uclahs-cds/package-ApplyPolygenicScore/issues) page.

## Contributors
For lists of contributors please visit [here](https://github.com/uclahs-cds/package-ApplyPolygenicScore/graphs/contributors) at GitHub.

## License

Author: Nicole Zeltser(nzeltser@mednet.ucla.edu)

package-ApplyPolygenicScore is licensed under the GNU General Public License version 2. See the file LICENSE.md for the terms of the GNU GPL license.

A package that provides utilities for the application of an existing polygenic score to a VCF.

Copyright (C) 2021 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
