# utility function to identify entries of multiallelic sites that should be excluded from
# pgs weight calculation due to mismatch of effect allele and genotype
# Input is a data.frame with VCF and pgs weight data, requires columns CHROM, POS, Indiv, gt_GT_alleles, effect_allele, and beta
# Output is the input data.frame filtered for non-risk multiallellic site entries
get.non.risk.multiallelic.site.row <- function(merged.vcf.with.pgs.data, merged.vcf.allele.matrix = NULL, vcf.long.format, original.df.row.index = NULL, current.chrom = NULL, current.pos = NULL) {
    if (vcf.long.format) {
        # handle case where there are no multiallelic sites
        if (nrow(merged.vcf.with.pgs.data) < 2) {
            return(data.frame());
            }

        # handle cases where this variant is missing in this individual sample
        if (all(merged.vcf.with.pgs.data$gt_GT_alleles == '.') | all(is.na(merged.vcf.with.pgs.data$gt_GT_alleles))) { #standard importation marks missing genotypes with '.'
            # if all alleles are missing (variant not called for this sample), the first entry of the multiallelic site
            # is arbitrarily chosen to represent the sample (any choice would be NA regardless)
            risk.allele.site.row.index <- 1;
            return(merged.vcf.with.pgs.data[-risk.allele.site.row.index, ]);
            }

        # split genotyped alleles into separate columns
        GT.alleles <- data.table::tstrsplit(merged.vcf.with.pgs.data$gt_GT_alleles, split = c('/|\\|'), keep = c(1,2));
        names(GT.alleles) <- c('called.allele.a', 'called.allele.b');
        GT.alleles <- data.frame(do.call(cbind, GT.alleles));

        # one of these alleles must be chosen to represent the sample, the rest will not be counted in pgs calculation
        # identify which multiallelic entry represents the effect allele that is present in the sample's genotype
        risk.allele.to.gt.matching <- GT.alleles == merged.vcf.with.pgs.data$effect_allele
        risk.allele.site.index <- which(risk.allele.to.gt.matching, arr.ind = TRUE);
        if (length(risk.allele.site.index) == 0) {
            # if neither effect allele is found in the sample genotype, the sample is homozygous for the non-risk allele
            # in this case, the first entry of the multiallelic site is arbitrarily chosen to represent the sample
            risk.allele.site.row.index <- 1;
            } else if (length(unique(risk.allele.site.index[ ,'row'])) > 1) {
                # if both effect alleles are found in the sample genotype, the sample is heterozygous for the effect alleles
                # this tool is not designed to handle multiple betas for one genotype
                # in this case, the entry with the higher effect alle beta is chosen to represent the sample
                risk.allele.site.row.index <- which.max(merged.vcf.with.pgs.data$beta);
                warning(paste0(
                    'Multiple effect alleles found in a sample, choosing effect allele with highest beta for dosage calculation. Check coordinates ',
                    merged.vcf.with.pgs.data$CHROM, ':', merged.vcf.with.pgs.data$POS, '\n'
                    ))
            } else {
                risk.allele.site.row.index <- risk.allele.site.index[1, 'row'];
            }
        non.risk.allele.entries <- merged.vcf.with.pgs.data[-risk.allele.site.row.index, ];
        return(non.risk.allele.entries);
        } else {
            # handle case where there are no multiallelic sites
            if (nrow(merged.vcf.with.pgs.data) < 2) {
                # return a matrix of input allele matrix dimensions with just TRUE (all alleles are non-risk)
                output.matrix <- NULL;
                return(output.matrix);
                }

            # Extract effect alleles from the fixed data
            effect.alleles <- merged.vcf.with.pgs.data$effect_allele;
            beta.values <- merged.vcf.with.pgs.data$beta;

            # Get the corresponding columns from the allele matrix (the external matrix)
            allele.matrix.subset <- merged.vcf.allele.matrix[original.df.row.index, , drop = FALSE];

            # Split alleles and find matches
            alleles.split <- data.table::tstrsplit(as.vector(allele.matrix.subset), split = '/|\\|');
            alleles.a <- matrix(alleles.split[[1]], nrow = nrow(allele.matrix.subset));
            alleles.b <- matrix(alleles.split[[2]], nrow = nrow(allele.matrix.subset));


            # Replicate effect alleles across columns
            effect.allele.matrix <- matrix(
                effect.alleles,
                nrow = length(effect.alleles),
                ncol = ncol(allele.matrix.subset)
                );

            # Determine if the effect allele is present in a vectorized way
            is.effect.allele.present.matrix <- (alleles.a == effect.allele.matrix) | (alleles.b == effect.allele.matrix);

            # Logic to find the single risk allele entry to keep for each sample
            # browser()
            risk.alleles.to.keep.index <- apply(
                X = is.effect.allele.present.matrix,
                MARGIN = 2, # Apply to each column (sample)
                FUN = function(is.present.col) {
                    # Find which row has a match
                    matched.rows <- which(is.present.col);

                    if (length(matched.rows) == 1) {
                        return(matched.rows);
                    } else if (length(matched.rows) == 0) {
                        return(1); # Arbitrarily choose the first entry
                    } else {
                        # Handle multiple matches: choose the one with the highest beta
                        warning(paste0(
                            'Multiple effect alleles found in a sample, choosing effect allele with highest beta for dosage calculation. ',
                            'Check coordinates: ', current.chrom, ':', current.pos, '\n'
                        ));
                        max.beta.index <- matched.rows[which.max(beta.values[matched.rows])];
                        return(max.beta.index);
                    }
                    }
                );

            # Create a final logical mask
            # This mask identifies which specific allele (row in the current subset)
            # for a given sample (column) is the "risk" allele chosen to be kept.
            risk.allele.selected.mask <- matrix(
                FALSE,
                nrow = nrow(allele.matrix.subset),
                ncol = ncol(allele.matrix.subset)
            );
            for (j in 1:ncol(risk.allele.selected.mask)) {
                risk.allele.selected.mask[risk.alleles.to.keep.index[j], j] <- TRUE;
            }

            # The non-risk mask is the inverse of the selected risk alleles
            final.non.risk.mask.for.site <- !risk.allele.selected.mask;

            # Return a data.table with the original row indices and the logical mask
            # Return a data.table with list columns.
            # This prevents data.table from unlisting the vector/matrix when combining results.
            return(data.table(
                original.df.row.index = list(original.df.row.index), # Wrap the vector of row indices in a list
                mask.matrix = list(final.non.risk.mask.for.site)     # Wrap the logical matrix in a list
                ));

        }
    }
