# package-ApplyPolygenicScore

A package that provides utilities for the application of an existing polygenic score to a VCF.
Matches coordinates of a provided polygenic score weight file to VCF input and calculates weighted sum of dosages in each individual.

## Description
This tool is intended to simply and transparently parse genotype/dosage data from an input VCF, match genotype coordinates to the component SNPs of an existing polygenic score, and apply SNP weights to dosages to calculate a polygenic score for each individual in accordance with the additive weighted sum of dosages model.

### Inputs
ApplyPolygenicScore needs only two inputs: 
1. Genotype data of the individuals upon which you wish to apply a polygenic score.
2. Coordinates of each SNP that compose the polygenic score you with to apply, and their associated weights (weight file).

#### Genotype data
Genotype data should be provided in the form of a VCF file. If you wish to apply a PGS to a cohort, we recommend that genotypes for the whole cohort be aggregated in one VCF file, ideally through a regenotyping process, or through VCF merging with a tool like bcftools. To reduce memory usage and improve speed of PGS application, we recommend pre-filtering the input VCF for only the coordinates that compose the PGS you wish to apply. This action can be performed using a coordinate BED file and tools such as bcftools or bedtools. To facilitate this process, ApplyPolygenicScore provides a function that outputs a BED file containing coordinates for any number of PGS weight files provided as input.

#### PGS weight file
[Discuss non-standardized sources including PGS catalog. Discuss importance of matching reference genome, related harmonization, liftover, and VCF formatting differences like "chr" prefix. Discuss weight format e.g. Odds Ratio vs beta.]

## License

Author: Nicole Zeltser(nzeltser@mednet.ucla.edu)

package-ApplyPolygenicScore is licensed under the GNU General Public License version 2. See the file LICENSE.md for the terms of the GNU GPL license.

A package that provides utilities for the application of an existing polygenic score to a VCF.

Copyright (C) 2021 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
