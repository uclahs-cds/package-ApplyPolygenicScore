# utility function to get the coordinates of multiallelic sites in VCF data
# Input is a data.frame with CHROM, POS, and Indiv columns
# Output is the input data frame filtered for the first (of 2 or more) entries of each multiallelic site
get.multiallelic.site.coordinates <- function(vcf.data) {
    single.sample.cohort.representative <- unique(vcf.data$Indiv)[1];
    single.sample.vcf.data <- vcf.data[vcf.data$Indiv == single.sample.cohort.representative, ];
    first.multiallelic.entry.index <- which(duplicated(paste0(single.sample.vcf.data$CHROM, single.sample.vcf.data$POS)));
    first.multiallelic.entries <- single.sample.vcf.data[first.multiallelic.entry.index, ];
    return(first.multiallelic.entries);
    }

# utility function to identify entries of multiallelic sites that should be excluded from
# pgs weight calculation due to mismatch of effect allele and genotype
# Input is a data.frame with VCF and pgs weight data, requires columns CHROM, POS, Indiv, gt_GT_alleles, effect_allele, and beta
# Output is the input data.frame filtered for non-risk multiallellic site entries
get.non.risk.multiallelic.site.row <- function(single.sample.multialellic.pgs.with.vcf.data) {
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
            risk.allele.site.row.index <- which.max(single.sample.multiallelic.vcf.data$beta);
            warning(paste0(
                'Multiple effect alleles found in sample ',
                unique(single.sample.multiallellic.vcf.data$Indiv),
                ' genotype, choosing effect allele with highest beta for dosage calculation. Check coordinates ',
                single.sample.multiallellic.vcf.data$CHROM, ':', single.sample.multiallellic.vcf.data$POS
                ))
        } else {
            risk.allele.site.row.index <- risk.allele.site.index[1, 'row'];
        }
    non.risk.allele.entries <- single.sample.multialellic.pgs.with.vcf.data[-risk.allele.site.row.index, ];
    return(non.risk.allele.entries);
    }
