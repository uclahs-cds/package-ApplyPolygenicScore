#' @title Apply polygenic score to VCF data
#' @description Apply a polygenic score to VCF data.
#' @param vcf.data A data.frame containing VCF genotype data.
#' @param pgs.weight.data A data.frame containing PGS weight data.
#' @return A data.frame containing the PGS per sample.
#' @export
apply.polygenic.score <- function(vcf.data, pgs.weight.data) {
    # check that inputs are data.frames
    if (!is.data.frame(vcf.data)) {
        stop('vcf.data must be a data.frame');
        }
    if (!is.data.frame(pgs.weight.data)) {
        stop('pgs.weight.data must be a data.frame');
        }

    # check that inputs contain required columns for PGS application
    required.vcf.columns <- c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles');
    required.pgs.columns <- c('CHROM', 'POS', 'effect_allele', 'beta');
    if (!all(required.vcf.columns %in% colnames(vcf.data))) {
        stop('vcf.data must contain columns named CHROM, POS, REF, ALT, Indiv, and gt_GT_alleles');
        }
    if (!all(required.pgs.columns %in% colnames(pgs.weight.data))) {
        stop('pgs.weight.data must contain columns named CHROM, POS, effect_allele, and beta');
        }

    # check for duplicate variants in PGS data
    if (any(duplicated(paste0(pgs.weight.data$CHROM, pgs.weight.data$POS)))) {
        stop('Duplicate variants are present in the PGS weight data. Please remove duplicate variants.');
        }

    # check that all samples have variant data represented for all variants
    n.samples <- length(unique(vcf.data$Indiv));
    n.variants <- length(unique(paste0(vcf.data$CHROM, vcf.data$POS, vcf.data$REF, vcf.data$ALT)));
    if (nrow(vcf.data) != n.samples * n.variants) {
        stop('Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.');
        }

    # merge VCF and PGS data
    merged.vcf.with.pgs <- merge.vcf.with.pgs(
        vcf.data = vcf.data,
        pgs.weight.data = pgs.weight.data
        );
    merged.vcf.with.pgs.data <- merged.vcf.with.pgs$merged.vcf.with.pgs.data;

    # identify multiallelic sites
    # if multiallelics are merged?
    # if multiallelics are split, remove the entry with risk allele

    # identify coordinates of multiallelic sites
    indiv.for.filter <- unique(merged.vcf.with.pgs.data$Indiv)[1];
    sites.one.sample <- merged.vcf.with.pgs.data[merged.vcf.with.pgs.data$Indiv == indiv.for.filter, ];
    first.multiallelic.site.index <- which(duplicated(paste0(sites.one.sample$CHROM, sites.one.sample$POS)));
    first.multiallelic.sites <- sites.one.sample[first.multiallelic.site.index, ];

    # identify which multiallelic sites have the risk allele
    for (i in 1:nrow(first.multiallelic.sites)) {
        # extract all of the same multiallelic site
        multiallelic.site.row.index <- which(merged.vcf.with.pgs.data$CHROM == first.multiallelic.sites$CHROM[i] & merged.vcf.with.pgs.data$POS == first.multiallelic.sites$POS[i]);
        multiallelic.site <- merged.vcf.with.pgs.data[multiallelic.site.row.index, ];
        # identify which of the duplicate positions contains the risk allele
        risk.allele.site.index <- which(multiallelic.site$effect_allele == multiallelic.site$REF | multiallelic.site$effect_allele == multiallelic.site$ALT);
        # remove position with non-risk allele
        multiallelic.site <- multiallelic.site[risk.allele.site.index, ];
        # if the reference allele is the risk allele, there can still be multiple entries
        # arbitrarily keep one of them
        # sort by Indiv then REF then ALT
        multiallelic.site <- multiallelic.site[order(multiallelic.site$Indiv, multiallelic.site$REF, multiallelic.site$ALT), ];
        # keep the first n entries, where n is the number of unique samples
        n.samples <- length(unique(multiallelic.site$Indiv));
        multiallelic.site <- multiallelic.site[1:n.samples, ];
        # replace the multiallelic site with the single allele site
        # first remove all multiallelic rows
        merged.vcf.with.pgs.data <- merged.vcf.with.pgs.data[-multiallelic.site.row.index, ];
        # then add the single allele rows
        merged.vcf.with.pgs.data <- rbind(merged.vcf.with.pgs.data, multiallelic.site);

    }

    # calculate dosage
    merged.vcf.with.pgs.data$dosage <- convert.alleles.to.pgs.dosage(
        called.alleles = merged.vcf.with.pgs.data$gt_GT_alleles,
        risk.alleles = merged.vcf.with.pgs.data$effect_allele
        );

    # calculate weighted dosage
    merged.vcf.with.pgs.data$weighted.dosage <- merged.vcf.with.pgs.data$dosage * merged.vcf.with.pgs.data$beta;

    # calculate PGS per sample using base R
    pgs.per.sample <- aggregate(
        x = merged.vcf.with.pgs.data$weighted.dosage,
        by = list(merged.vcf.with.pgs.data$Indiv),
        FUN = sum,
        na.rm = TRUE
        );
    colnames(pgs.per.sample) <- c('sample', 'PGS');

    return(pgs.per.sample);
    }
