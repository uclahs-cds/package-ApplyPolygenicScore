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
        # sort by Indiv then REF then ALT
        multiallelic.site <- multiallelic.site[order(multiallelic.site$Indiv, multiallelic.site$REF, multiallelic.site$ALT), ];

        # merge multiallelic sites, one sample at a time
        merged.multiallelic.sites <- sapply(
            X = unique(multiallelic.site$Indiv),
            FUN = function(x) {
                merge.multiallelic.site(multiallelic.site[multiallelic.site$Indiv == x, ])
                }
            )

        # remove position with non-risk/non-selected allele
        multiallelic.site <- multiallelic.site[risk.allele.site.index, ];

        
        # # keep the first n entries, where n is the number of unique samples
        # n.samples <- length(unique(multiallelic.site$Indiv));
        # multiallelic.site <- multiallelic.site[1:n.samples, ];
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

merge.multiallelic.site <- function(multiallelic.site) {
    # check for required columns
    required.columns <- c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles');
    if (!all(required.columns %in% colnames(multiallelic.site))) {
        stop('multiallelic.site must contain columns named CHROM, POS, REF, ALT, and gt_GT_alleles');
        }

    # check that only two extra alleles are present
    if (length(unique(multiallelic.site$ALT)) != 2) {
        stop('multiallelic.site must contain exactly two unique ALT alleles');
        }

    # check that only one sample is present
    if (length(unique(multiallelic.site$Indiv)) != 1) {
        stop('multiallelic.site must contain exactly one unique sample');
        }

    # check that only two rows are present
    if (nrow(multiallelic.site) != 2) {
        stop('multiallelic.site must contain exactly two rows');
        }

    # split GT alleles
    split.alleles <- data.table::tstrsplit(called.alleles, split = c('/|\\|'), keep = c(1,2));
    names(split.alleles) <- c('called.allele.a', 'called.allele.b');

    # define three possible alleles
    ref.allele <- unique(multiallelic.site$REF);
    alt1.allele <- multiallelic.site$ALT[1];
    alt2.allele <- multiallelic.site$ALT[2];

    all.alleles <- c(split.alleles.called.allele.a, split.alleles.called.allele.b);

    # count alt alleles
    n.alt1.alleles <- sum(alt1.allele %in% all.alleles);
    n.alt2.alleles <- sum(alt2.allele %in% all.alleles);

    if (unique(all.alleles) == ref.allele) {
        # if all alleles are REF, the merged genotype is REF/REF
        merged.gt <- paste0(ref.allele, '/', ref.allele);
        } else if (n.alt1.alleles == 1 & n.alt2.alleles == 0) {
            # if one allele is ALT1 and no alleles are ALT2, the merged genotype is REF/ALT1
            merged.gt <- paste0(ref.allele, '/', alt1.allele);
        } else if (n.alt1.alleles == 0 & n.alt2.alleles == 1) {
            # if no alleles are ALT1 and one allele is ALT2, the merged genotype is REF/ALT2
            merged.gt <- paste0(ref.allele, '/', alt2.allele);
        } else if (n.alt1.alleles == 1 & n.alt2.alleles == 1) {
            # if one allele is ALT1 and one allele is ALT2, the merged genotype is ALT1/ALT2
            merged.gt <- paste0(alt1.allele, '/', alt2.allele);
        } else if (n.alt1.alleles == 2 & n.alt2.alleles == 0) {
            # if two alleles are ALT1 and no alleles are ALT2, the merged genotype is ALT1/ALT1
            merged.gt <- paste0(alt1.allele, '/', alt1.allele);
        } else if (n.alt1.alleles == 0 & n.alt2.alleles == 2) {
            # if no alleles are ALT1 and two alleles are ALT2, the merged genotype is ALT2/ALT2
            merged.gt <- paste0(alt2.allele, '/', alt2.allele);
        } else {
            stop('unrecognized multiallelic site format');
        }

    # assemble merged site
    merged.site <- multiallelic.site[1, ];
    merged.site$ALT <- paste0(alt1.allele, ',', alt2.allele);
    merged.site$gt_GT_alleles <- merged.gt;

    return(merged.site);
    }
