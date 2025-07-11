# utility function to identify entries of multiallelic sites that should be excluded from
# pgs weight calculation due to mismatch of effect allele and genotype
# Input is a data.frame with VCF and pgs weight data, requires columns CHROM, POS, Indiv, gt_GT_alleles, effect_allele, and beta
# Output is the input data.frame filtered for non-risk multiallellic site entries
get.non.risk.multiallelic.site.row <- function(single.sample.multialellic.pgs.with.vcf.data, vcf.long.format) {
    if (vcf.long.format) {
        # handle case where there are no multiallelic sites
        if (nrow(single.sample.multialellic.pgs.with.vcf.data) < 2) {
            return(data.frame());
            }

        # handle cases where this variant is missing in this individual sample
        if (all(single.sample.multialellic.pgs.with.vcf.data$gt_GT_alleles == '.') | all(is.na(single.sample.multialellic.pgs.with.vcf.data$gt_GT_alleles))) { #standard importation marks missing genotypes with '.'
            # if all alleles are missing (variant not called for this sample), the first entry of the multiallelic site
            # is arbitrarily chosen to represent the sample (any choice would be NA regardless)
            risk.allele.site.row.index <- 1;
            return(single.sample.multialellic.pgs.with.vcf.data[-risk.allele.site.row.index, ]);
            }

        # split genotyped alleles into separate columns
        GT.alleles <- data.table::tstrsplit(single.sample.multialellic.pgs.with.vcf.data$gt_GT_alleles, split = c('/|\\|'), keep = c(1,2));
        names(GT.alleles) <- c('called.allele.a', 'called.allele.b');
        GT.alleles <- data.frame(do.call(cbind, GT.alleles));

        # one of these alleles must be chosen to represent the sample, the rest will not be counted in pgs calculation
        # identify which multiallelic entry represents the effect allele that is present in the sample's genotype
        risk.allele.to.gt.matching <- GT.alleles == single.sample.multialellic.pgs.with.vcf.data$effect_allele
        risk.allele.site.index <- which(risk.allele.to.gt.matching, arr.ind = TRUE);
        if (length(risk.allele.site.index) == 0) {
            # if neither effect allele is found in the sample genotype, the sample is homozygous for the non-risk allele
            # in this case, the first entry of the multiallelic site is arbitrarily chosen to represent the sample
            risk.allele.site.row.index <- 1;
            } else if (length(unique(risk.allele.site.index[ ,'row'])) > 1) {
                # if both effect alleles are found in the sample genotype, the sample is heterozygous for the effect alleles
                # this tool is not designed to handle multiple betas for one genotype
                # in this case, the entry with the higher effect alle beta is chosen to represent the sample
                risk.allele.site.row.index <- which.max(single.sample.multialellic.pgs.with.vcf.data$beta);
                warning(paste0(
                    'Multiple effect alleles found in ',
                    unique(single.sample.multialellic.pgs.with.vcf.data$Indiv),
                    ' genotype, choosing effect allele with highest beta for dosage calculation. Check coordinates ',
                    single.sample.multialellic.pgs.with.vcf.data$CHROM, ':', single.sample.multialellic.pgs.with.vcf.data$POS, '\n'
                    ))
            } else {
                risk.allele.site.row.index <- risk.allele.site.index[1, 'row'];
            }
        non.risk.allele.entries <- single.sample.multialellic.pgs.with.vcf.data[-risk.allele.site.row.index, ];
        return(non.risk.allele.entries);
        } else {
            # handle case where there are no multiallelic sites
            if (nrow(single.sample.multialellic.pgs.with.vcf.data$fixed.data) < 2) {
                # return a matrix of input allele matrix dimensions with just TRUE (all alleles are non-risk)
                output.matrix <- NULL;
                return(output.matrix);
                }
            # # Chat GPT suggestion
            # test.df <- data.frame(
            #     CHROM = c('chr1', 'chr1'),
            #     POS = c(100, 100),
            #     effect_allele = c('A', 'C'),
            #     beta = c(1,2)
            #     );
            # test.allele.matrix <- matrix(
            #     c('C/A', 'C/A', 'T/C', 'T/C'),
            #     nrow = 2,
            #     ncol = 2
            #     );
            # single.sample.multialellic.pgs.with.vcf.data = list(fixed.data = test.df, allele.matrix = test.allele.matrix);

            # Extract effect alleles from the fixed data
            effect.alleles <- single.sample.multialellic.pgs.with.vcf.data$fixed.data$effect_allele;

            # Expand effect_alleles across allele matrix dimensions
            effect.allele.matrix <- matrix(effect.alleles, nrow = nrow(single.sample.multialellic.pgs.with.vcf.data$fixed.data), ncol = ncol(single.sample.multialellic.pgs.with.vcf.data$allele.matrix));
            # Extract alleles from genotypes using regex split
            allele.list <- strsplit(single.sample.multialellic.pgs.with.vcf.data$allele.matrix, split = '/|\\|');

            # Determine whether effect allele is present in each genotype
            is.effect.allele.present.vector <- mapply(
                FUN = function(alleles, effect) {
                    # if an allele is missing or NA, return FALSE
                    if (any(alleles == '.' | is.na(alleles))) return(FALSE);
                    # otherwise return whether the effect allele is present in the corresponding genotype
                    effect %in% alleles;
                    },
                    allele.list,
                    as.vector(effect.allele.matrix)
                    );

            # Reshape logical result to same matrix dims as allele.matrix
            is.effect.allele.present.matrix <- matrix(is.effect.allele.present.vector, nrow = nrow(single.sample.multialellic.pgs.with.vcf.data$allele.matrix), ncol = ncol(single.sample.multialellic.pgs.with.vcf.data$allele.matrix));

            # Apply conditions for various multiallelic scenarios
            risk.allele.site.matrix.index <- apply(
                X = is.effect.allele.present.matrix,
                MARGIN = 2,
                FUN = function(col) {
                    if (sum(col) == 1) {
                        return(col); # one risk/effect allele present, return the column as is
                        }
                    if (sum(col) == 0) {
                        # Both alleles are non-risk, arbitrarily choose the first entry as the risk allele
                        col <- rep(FALSE, length(col)); # Set all to FALSE
                        col[1] <- TRUE; # Set the first entry to TRUE
                        return(col);
                        }
                    if (sum(col) > 1) {
                        # More than one effect allele present, choose the one with the highest beta as the risk allele
                        max.beta.index <- which.max(single.sample.multialellic.pgs.with.vcf.data$fixed.data$beta[col]);
                        col[col] <- FALSE; # Set all to FALSE
                        col[max.beta.index] <- TRUE; # Set the one with max beta to TRUE
                        return(col);
                        }
                    }
                );
            # Return a matrix with the same dimension as the allele matrix, with TRUE for non-risk alleles
            output.matrix <- !risk.allele.site.matrix.index;
            return(output.matrix);

        }
    }
