# utility function to identify entries of multiallelic sites that should be excluded from
# pgs weight calculation due to mismatch of effect allele and genotype
# Input is a data.frame with VCF and pgs weight data, requires columns CHROM, POS, Indiv, gt_GT_alleles, effect_allele, and beta
# Output is the input data.frame filtered for non-risk multiallellic site entries
get.non.risk.multiallelic.site.row <- function(single.sample.multialellic.pgs.with.vcf.data) {
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
        # in this case, the first entry of the multiallelic site is arbitratily chosen to represent the sample
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
    }
