#' @title Apply polygenic score to VCF data
#' @description Apply a polygenic score to VCF data.
#' @param vcf.data A data.frame containing VCF genotype data.
#' @param pgs.weight.data A data.frame containing PGS weight data.
#' @param missing.genotype.method A character string indicating the method to handle missing genotypes. Options are "mean.dosage", "normalize", or "none". Default is "mean.dosage".
#' @return A data.frame containing the PGS per sample.
#' @export
apply.polygenic.score <- function(vcf.data, pgs.weight.data, missing.genotype.method = 'mean.dosage') {
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
    if (any(duplicated(paste0(pgs.weight.data$CHROM, pgs.weight.data$POS, pgs.weight.data$effect_allele)))) {
        stop('Duplicate variants detected in the PGS weight data. Please ensure only unique coordinate:effect allele combinations are present.');
        }

    # check for duplicate coordinates in PGS data
    if (any(duplicated(paste0(pgs.weight.data$CHROM, pgs.weight.data$POS)))) {
        warning('Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.');
        }

    # check that all samples have variant data represented for all variants
    n.samples <- length(unique(vcf.data$Indiv));
    n.variants <- length(unique(paste0(vcf.data$CHROM, vcf.data$POS, vcf.data$REF, vcf.data$ALT)));
    if (nrow(vcf.data) != n.samples * n.variants) {
        stop('Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.');
        }

    # check missing genotype method input
    if (all(missing.genotype.method %in% c('mean.dosage', 'normalize', 'none'))) {
        missing.genotype.method <- missing.genotype.method;
        } else {
        stop('missing.genotype.method must be either "mean.dosage", "normalize", or "none"');
        }
    if ('none' %in% missing.genotype.method && length(missing.genotype.method) > 1) {
        stop('If "none" is included in missing.genotype.method, it must be the only method included');
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

    ### Start Missing Genotype Handling ###

    # create sample by variant dosage matrix
    variant.id <- paste(merged.vcf.with.pgs.data$CHROM, merged.vcf.with.pgs.data$POS, merged.vcf.with.pgs.data$REF, merged.vcf.with.pgs.data$effect_allele, sep = ':');
    dosage.matrix <- get.variant.by.sample.matrix(
        long.data = merged.vcf.with.pgs.data,
        variant.id = variant.id,
        value.var = 'dosage'
        );

    if ('mean.dosage' %in% missing.genotype.method) {
        # calculate dosage to replace missing genotypes
        missing.genotype.dosage <- calculate.missing.genotype.dosage(dosage.matrix = dosage.matrix);

        # identify missing genotypes
        missing.genotype.row.index <- which(is.na(merged.vcf.with.pgs.data$dosage) & !is.na(merged.vcf.with.pgs.data$Indiv));
        # start a column for replaced missing dosages
        merged.vcf.with.pgs.data$dosage.with.replaced.missing <- merged.vcf.with.pgs.data$dosage;
        # assign mean dosage to missing genotypes
        for (i in missing.genotype.row.index) {
            missing.variant.id <- paste(merged.vcf.with.pgs.data[i, 'CHROM'], merged.vcf.with.pgs.data[i, 'POS'], merged.vcf.with.pgs.data[i, 'REF'], merged.vcf.with.pgs.data[i, 'effect_allele'], sep = ':');
            missing.variant.dosage <- missing.genotype.dosage[missing.variant.id];
            merged.vcf.with.pgs.data[i, 'dosage.with.replaced.missing'] <- missing.variant.dosage;
            }
    }

    ### End Misssing Genotype Handling ###

    # calculate weighted dosage
    if ('mean.dosage' %in% missing.genotype.method) {
        merged.vcf.with.pgs.data$weighted.dosage.with.replaced.missing <- merged.vcf.with.pgs.data$dosage.with.replaced.missing * merged.vcf.with.pgs.data$beta;
        }
    if ('normalize' %in% missing.genotype.method || 'none' %in% missing.genotype.method) {
        merged.vcf.with.pgs.data$weighted.dosage <- merged.vcf.with.pgs.data$dosage * merged.vcf.with.pgs.data$beta;
        }


    ### Start Multiallelic Site Handling ###
    # create a dictionary to each unique sample:coordinate combination
    sample.coordinate.to.row.dict.hash <- new.env(hash = TRUE, parent = emptyenv());

    for (i in 1:nrow(merged.vcf.with.pgs.data)) {
        key <- paste(merged.vcf.with.pgs.data[i, 'Indiv'], merged.vcf.with.pgs.data[i, 'CHROM'], merged.vcf.with.pgs.data[i, 'POS'], sep = '_');
        sample.coordinate.to.row.dict.hash[[key]] <- c(sample.coordinate.to.row.dict.hash[[key]], i);
        }

    non.risk.multiallelic.entries.index <- lapply(
        X = ls(sample.coordinate.to.row.dict.hash),
        FUN = function(x) {
            row.index <- sample.coordinate.to.row.dict.hash[[x]];
            single.sample.multialellic.pgs.with.vcf.data <- merged.vcf.with.pgs.data[row.index, ];
            single.sample.multialellic.pgs.with.vcf.data$original.df.row.index <- row.index;
            non.risk.multiallelic.site.rows <- get.non.risk.multiallelic.site.row(
                single.sample.multialellic.pgs.with.vcf.data = single.sample.multialellic.pgs.with.vcf.data
                );
            return(non.risk.multiallelic.site.rows$original.df.row.index);
            }
        );

    non.risk.multiallelic.entries.index <- unlist(non.risk.multiallelic.entries.index);

    if ('mean.dosage' %in% missing.genotype.method) {
        merged.vcf.with.pgs.data$multiallelic.weighted.dosage.with.replaced.missing <- merged.vcf.with.pgs.data$weighted.dosage.with.replaced.missing;
        merged.vcf.with.pgs.data$multiallelic.weighted.dosage.with.replaced.missing[non.risk.multiallelic.entries.index] <- NA;
        }
    if ('normalize' %in% missing.genotype.method || 'none' %in% missing.genotype.method) {
        merged.vcf.with.pgs.data$multiallelic.weighted.dosage <- merged.vcf.with.pgs.data$weighted.dosage;
        merged.vcf.with.pgs.data$multiallelic.weighted.dosage[non.risk.multiallelic.entries.index] <- NA;
        }

    ### End Multiallelic Site Handling ###

    # calculate PGS per sample
    pgs.output.list <- list();

    if ('none' %in% missing.genotype.method) {
        pgs.per.sample <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample) <- c('sample', 'PGS');
        return(pgs.per.sample);
        }

    if ('normalize' %in% missing.genotype.method) {
        pgs.per.sample.with.normalized.missing <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample.with.normalized.missing) <- c('sample', 'PGS');
        bialellic.variant.id <- paste(merged.vcf.with.pgs.data$CHROM, merged.vcf.with.pgs.data$POS, sep = ':');
        biallelic.snp.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = merged.vcf.with.pgs.data,
            variant.id = bialellic.variant.id,
            value.var = 'multiallelic.weighted.dosage'
            );
        per.sample.non.missing.genotype.count <- colSums(!is.na(biallelic.snp.by.sample.matrix));
        pgs.per.sample.with.normalized.missing$PGS <- pgs.per.sample.with.normalized.missing$PGS / per.sample.non.missing.genotype.count;
        # account for division by zero
        pgs.per.sample.with.normalized.missing$PGS[is.nan(pgs.per.sample.with.normalized.missing$PGS)] <- 0;
        pgs.output.list$PGS.with.normalized.missing <- pgs.per.sample.with.normalized.missing;
        }

    if ('mean.dosage' %in% missing.genotype.method) {
        pgs.per.sample <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage.with.replaced.missing,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample) <- c('sample', 'PGS');
        pgs.output.list$PGS.with.replaced.missing <- pgs.per.sample
        }

    # format output
    # bind PGS columns of list components
    PGS.cols <- lapply(pgs.output.list, function(x) x$PGS);
    PGS.cols <- data.frame(do.call(cbind, PGS.cols));
    colnames(PGS.cols) <- names(pgs.output.list);
    # bind sample column of first list component
    sample <- pgs.output.list[[1]]$sample;
    pgs.output <- cbind(sample, PGS.cols);
    return(pgs.output);


    return(pgs.per.sample);
    }
