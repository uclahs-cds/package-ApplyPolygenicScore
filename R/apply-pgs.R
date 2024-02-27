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
        warning('Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.');
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

    # calculate dosage
    merged.vcf.with.pgs.data$dosage <- convert.alleles.to.pgs.dosage(
        called.alleles = merged.vcf.with.pgs.data$gt_GT_alleles,
        risk.alleles = merged.vcf.with.pgs.data$effect_allele
        );

    # calculate weighted dosage
    merged.vcf.with.pgs.data$weighted.dosage <- merged.vcf.with.pgs.data$dosage * merged.vcf.with.pgs.data$beta;


    # identify coordinates of multiallelic sites
    indiv.for.filter <- unique(merged.vcf.with.pgs.data$Indiv)[1];
    sites.one.sample <- merged.vcf.with.pgs.data[merged.vcf.with.pgs.data$Indiv == indiv.for.filter, ];
    first.multiallelic.site.index <- which(duplicated(paste0(sites.one.sample$CHROM, sites.one.sample$POS)));
    first.multiallelic.sites <- sites.one.sample[first.multiallelic.site.index, ];

    extracted.non.risk.multiallelic.entries <- list();
    # identify which multiallelic sites have the risk allele
    for (i in 1:nrow(first.multiallelic.sites)) {
        # extract all of the same multiallelic site
        multiallelic.site.row.index <- which(merged.vcf.with.pgs.data$CHROM == first.multiallelic.sites$CHROM[i] & merged.vcf.with.pgs.data$POS == first.multiallelic.sites$POS[i]);
        multiallelic.site <- merged.vcf.with.pgs.data[multiallelic.site.row.index, ];
        # sort by Indiv then REF then ALT
        multiallelic.site <- multiallelic.site[order(multiallelic.site$Indiv, multiallelic.site$REF, multiallelic.site$ALT), ];

        # iterate by sample through the multiallelic site
        for (j in 1:length(unique(multiallelic.site$Indiv))) {
            single.sample.multiallelic.site <- multiallelic.site[multiallelic.site$Indiv == unique(multiallelic.site$Indiv)[j], ];

            # split GT alleles into separate columns
            GT.alleles <- data.table::tstrsplit(single.sample.multiallelic.site$gt_GT_alleles, split = c('/|\\|'), keep = c(1,2));
            names(GT.alleles) <- c('called.allele.a', 'called.allele.b');
            GT.alleles <- data.frame(do.call(cbind, GT.alleles));

            # one of these alleles must be chosen to represent the sample, the rest will not be counted
            risk.allele.to.gt.matching <- GT.alleles == single.sample.multiallelic.site$effect_allele
            risk.allele.site.index <- which(risk.allele.to.gt.matching, arr.ind = TRUE);
            non.risk.allele.entry <- single.sample.multiallelic.site[-risk.allele.site.index[1, 'row'], ];
            extracted.non.risk.multiallelic.entries[[length(extracted.non.risk.multiallelic.entries) + 1]] <- non.risk.allele.entry;

            }
        }
    extracted.non.risk.multiallelic.entries <- do.call(rbind, extracted.non.risk.multiallelic.entries);
    # row.match returns index of first match of each row of x in table, returns NA for no match
    non.risk.multiallelic.entries.index <- prodlim::row.match(x = extracted.non.risk.multiallelic.entries, table = merged.vcf.with.pgs.data);

    merged.vcf.with.pgs.data$multiallelic.weighted.dosage <- merged.vcf.with.pgs.data$weighted.dosage;
    merged.vcf.with.pgs.data$multiallelic.weighted.dosage[non.risk.multiallelic.entries.index] <- NA;

    # calculate PGS per sample using base R
    pgs.per.sample <- aggregate(
        x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage,
        by = list(merged.vcf.with.pgs.data$Indiv),
        FUN = sum,
        na.rm = TRUE
        );
    colnames(pgs.per.sample) <- c('sample', 'PGS');

    return(pgs.per.sample);
    }
