
validate.vcf.input <- function(vcf.data) {
    # check that inputs are data.frames
    if (!is.data.frame(vcf.data)) {
        stop('vcf.data must be a data.frame');
        }

    # check that vcf.data contains required columns
    required.vcf.columns <- c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles');

    if (!all(required.vcf.columns %in% colnames(vcf.data))) {
        stop('vcf.data must contain columns named CHROM, POS, REF, ALT, Indiv, and gt_GT_alleles');
        }

    # check that all samples have variant data represented for all variants
    n.samples <- length(unique(vcf.data$Indiv));
    n.variants <- length(unique(paste(vcf.data$CHROM, vcf.data$POS, vcf.data$REF, vcf.data$ALT, sep = ':')));
    if (nrow(vcf.data) != n.samples * n.variants) {
        stop('Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.');
        }

    }

validate.pgs.data.input <- function(pgs.weight.data, use.external.effect.allele.frequency, correct.strand.flips, remove.ambiguous.allele.matches, remove.mismatched.indels) {
    if (!is.data.frame(pgs.weight.data)) {
        stop('pgs.weight.data must be a data.frame');
        }

    required.pgs.columns <- c('CHROM', 'POS', 'effect_allele', 'beta');

    if (!all(required.pgs.columns %in% colnames(pgs.weight.data))) {
        stop('pgs.weight.data must contain columns named CHROM, POS, effect_allele, and beta');
        }

    # additional required columns if strand flip correction is enabled
    if (correct.strand.flips || remove.ambiguous.allele.matches || remove.mismatched.indels) {
        if (!('other_allele' %in% colnames(pgs.weight.data))) {
            stop('pgs.weight.data must contain a column named other_allele if correct.strand.flips, remove.ambiguous.allele.matches, or remove.mismatched.indels is TRUE');
            }
        }

    if (use.external.effect.allele.frequency) {
        required.eaf.column <- 'allelefrequency_effect';
        if (!(required.eaf.column %in% colnames(pgs.weight.data))) {
            stop('pgs.weight.data must contain a column named allelefrequency_effect if use.external.effect.allele.frequency is TRUE');
            }
        }

    # check for duplicate variants in PGS data
    if (any(duplicated(paste(pgs.weight.data$CHROM, pgs.weight.data$POS, pgs.weight.data$effect_allele, sep = ':')))) {
        stop('Duplicate variants detected in the PGS weight data. Please ensure only unique coordinate:effect allele combinations are present.');
        }

    # check for duplicate coordinates in PGS data
    if (any(duplicated(paste(pgs.weight.data$CHROM, pgs.weight.data$POS, sep = ':')))) {
        warning('Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.');
        }
    }

validate.phenotype.data.input <- function(phenotype.data, phenotype.analysis.columns, vcf.data) {
    if (!is.null(phenotype.data)) {
        if (!is.data.frame(phenotype.data)) {
            stop('phenotype.data must be a data.frame');
            }

        required.phenotype.columns <- 'Indiv';

        if (!all(required.phenotype.columns %in% colnames(phenotype.data))) {
            stop('phenotype.data must contain columns named Indiv');
            }

        # check for at least one matching Indiv between phenotype.data and vcf.data
        if (length(intersect(phenotype.data$Indiv, vcf.data$Indiv)) == 0) {
            stop('No matching Indiv between phenotype.data and vcf.data');
            }

        # validate phenotype.analysis.columns if provided
        if (!is.null(phenotype.analysis.columns)) {
            if (!all(phenotype.analysis.columns %in% colnames(phenotype.data))) {
                stop('phenotype.analysis.columns must be columns in phenotype.data');
                }
            }

        } else if (!is.null(phenotype.analysis.columns)) {
            stop('phenotype.analysis.columns provided but no phenotype data detected');
            }

    }

#' @title Apply polygenic score to VCF data
#' @description Apply a polygenic score to VCF data.
#' @param vcf.data A data.frame containing VCF genotype data as formatted by \code{import.vcf()}.
#' @param pgs.weight.data A data.frame containing PGS weight data as formatted by \code{import.pgs.weight.file()}.
#' @param phenotype.data A data.frame containing phenotype data. Must have an Indiv column matching vcf.data. Default is \code{NULL}.
#' @param phenotype.analysis.columns A character vector of phenotype columns from phenotype.data to analyze in a regression analsyis. Default is \code{NULL}.
#' Phenotype variables are automatically classified as continuous, binary, or neither based on data type and number of unique values. The calculated PGS is associated
#' with each phenotype variable using linear or logistic regression for continuous or binary phenotypes, respectively. See \code{run.pgs.regression} for more details.
#' If no phenotype.analysis.columns are provided, no regression analysis is performed.
#' @param correct.strand.flips A logical indicating whether to check PGS weight data/VCF genotype data matches for strand flips and correct them. Default is \code{TRUE}.
#' The PGS catalog standard column \code{other_allele} in \code{pgs.weight.data} is required for this check.
#' @param remove.ambiguous.allele.matches A logical indicating whether to remove PGS variants with ambiguous allele matches between PGS weight data and VCF genotype data. Default is \code{FALSE}.
#' The PGS catalog standard column \code{other_allele} in \code{pgs.weight.data} is required for this check.
#' @param remove.mismatched.indels A logical indicating whether to remove indel variants that are mismatched between PGS weight data and VCF genotype data. Default is \code{FALSE}.
#' The PGS catalog standard column \code{other_allele} in \code{pgs.weight.data} is required for this check.
#' @param output.dir A character string indicating the directory to write output files. Separate files are written for per-sample pgs results and optional regression results.
#' Files are tab-separate .txt files. Default is NULL in which case no files are written.
#' @param file.prefix A character string to prepend to the output file names. Default is \code{NULL}.
#' @param missing.genotype.method A character string indicating the method to handle missing genotypes. Options are "mean.dosage", "normalize", or "none". Default is "mean.dosage".
#' @param use.external.effect.allele.frequency A logical indicating whether to use an external effect allele frequency for calculating mean dosage when handling missing genotypes. Default is \code{FALSE}.
#' Provide allele frequency as a column is \code{pgs.weight.data} named \code{allelefrequency_effect}.
#' @param n.percentiles An integer indicating the number of percentiles to calculate for the PGS. Default is \code{NULL}.
#' @param analysis.source.pgs A character string indicating the source PGS for percentile calculation and regression analyses. Options are "mean.dosage", "normalize", or "none".
#' When not specified, defaults to \code{missing.genotype.method} choice and if more than one PGS missing genotype method is chosen, calculation defaults to the first selection.
#' @param validate.inputs.only A logical indicating whether to only perform input data validation checks without running PGS application.
#' If no errors are triggered, a message is printed and TRUE is returned. Default is \code{FALSE}.
#' @return A list containing per-sample PGS output and per-phenotype regression output if phenotype analysis columns are provided.
#'
#' \strong{Output Structure}
#'
#' The outputed list contains the following elements:
#' \itemize{
#' \item pgs.output: A data.frame containing the PGS per sample and optional phenotype data.
#' \item regression.output: A data.frame containing the results of the regression analysis if phenotype.analysis.columns are provided, otherwise \code{NULL}.
#' }
#'
#' pgs.output columns:
#' \itemize{
#' \item \code{Indiv}: A character string indicating the sample ID.
#' \item \code{PGS}: A numeric vector indicating the PGS per sample. (only if missing.genotype.method includes "none")
#' \item \code{PGS.with.normalized.missing}: A numeric vector indicating the PGS per sample with missing genotypes normalized. (only if missing.genotype.method includes "normalize")
#' \item \code{PGS.with.replaced.missing}: A numeric vector indicating the PGS per sample with missing genotypes replaced by mean dosage. (only if missing.genotype.method includes "mean.dosage")
#' \item \code{percentile}: A numeric vector indicating the percentile rank of the PGS.
#' \item \code{decile}: A numeric vector indicating the decile rank of the PGS.
#' \item \code{quartile}: A numeric vector indicating the quartile rank of the PGS.
#' \item \code{percentile.X:} A numeric vector indicating the user-specified percentile rank of the PGS where "X" is substituted by \code{n.percentiles}. (only if \code{n.percentiles} is specified)
#' \item \code{n.missing.genotypes}: A numeric vector indicating the number of missing genotypes per sample.
#' \item \code{percent.missing.genotypes}: A numeric vector indicating the percentage of missing genotypes per sample.
#' \item All columns in \code{phenotype.data} if provided.
#' }
#'
#' regression.output columns:
#' \itemize{
#' \item phenotype: A character vector of phenotype names.
#' \item \code{model}: A character vector indicating the regression model used. One of "logistic.regression" or "linear.regression".
#' \item \code{beta}: A numeric vector indicating the beta coefficient of the regression analysis.
#' \item \code{se}: A numeric vector indicating the standard error of the beta coefficient.
#' \item \code{p.value}: A numeric vector indicating the p-value of the beta coefficient.
#' \item \code{r.squared}: A numeric vector indicating the r-squared value of linear regression analysis. NA for logistic regression.
#' \item \code{AUC}: A numeric vector indicating the area under the curve of logistic regression analysis. NA for linear regression.
#' }
#'
#' \strong{PGS Calculation}
#'
#' PGS for each individual \emph{i} is calculated as the sum of the product of the dosage and beta coefficient for each variant in the PGS:
#' \deqn{PGS_i = \sum_{m=1}^{M} \left( \beta_m \times dosage_{im} \right)}
#' Where \emph{m} is a PGS component variant out of a total \emph{M} variants.
#'
#' \strong{Missing Genotype Handling}
#'
#' VCF genotype data are matched to PGS data by chromosome and position. If a SNP cannot be matched by genomic coordinate,
#' an attempt is made to match by rsID (if available). If a SNP from the PGS weight data is not found in the VCF data after these two matching attempts,
#' it is considered a cohort-wide missing variant.
#'
#' Missing genotypes (in individual samples) among successfully matched variants are handled by three methods:
#'
#' \code{none}: Missing genotype dosages are excluded from the PGS calculation.
#' This is equivalent to assuming that all missing genotypes are homozygous for the non-effect allele, resulting in a dosage of 0.
#'
#' \code{normalize}: Missing genotypes are excluded from score calculation but the final score is normalized by the number of non-missing alleles.
#' The calculation assumes a diploid genome:
#' \deqn{PGS_i = \dfrac{\sum \left( \beta_m \times dosage_{im} \right)}{P_i * M_{non-missing}}}
#' Where \emph{P} is the ploidy and has the value \code{2} and \eqn{M_{non-missing}} is the number of non-missing genotypes.
#'
#' \code{mean.dosage}: Missing genotype dosages are replaced by the mean population dosage of the variant which is calculated as the product of the effect allele frequency \emph{EAF} and the ploidy of a diploid genome:
#' \deqn{\overline{dosage_{k}} = EAF_k * P}
#' where \emph{k} is a PGS component variant that is missing in between 1 and n-1 individuals in the cohort and \emph{P} = ploidy = 2
#' This dosage calculation holds under assumptions of Hardy-Weinberg equilibrium.
#' By default, the effect allele frequency is calculated from the provided VCF data.
#' For variants that are missing in all individuals (cohort-wide), dosage is assumed to be zero (homozygous non-reference) for all individuals.
#' An external allele frequency can be provided in the \code{pgs.weight.data} as a column named \code{allelefrequency_effect} and by setting \code{use.external.effect.allele.frequency} to \code{TRUE}.
#'
#' \strong{Multiallelic Site Handling}
#'
#' If a PGS weight file provides weights for multiple effect alleles, the appropriate dosage is calculated for the alleles that each individual carries.
#' It is assumed that multiallelic variants are encoded in the same row in the VCF data. This is known as "merged" format. Split multiallelic sites are not accepted.
#' VCF data can be formatted to merged format using external tools for VCF file manipulation.
#'
#' \strong{Allele Mismatch Handling}
#' Variants from the PGS weight data are merged with records in the VCF data by genetic coordinate.
#' After the merge is complete, there may be cases where the VCF reference (REF) and alternative (ALT) alleles do not match their conventional counterparts in the
#' PGS weight data (other allele and effect allele, respectively).
#' This is usually caused by a strand flip: the variant in question was called against opposite DNA reference strands in the PGS training data and the VCF data.
#' Strand flips can be detected and corrected by flipping the affected allele to its reverse complement.
#' \code{apply.polygenic.score} uses \code{assess.pgs.vcf.allele.match} to assess allele concordance, and is controlled through the following arguments:
#'
#' \itemize{
#' \item \code{correct.strand.flips}: When \code{TRUE}, detected strand flips are corrected by flipping the affected value in the \code{effect_allele} column prior to dosage calling.
#' \item \code{remove.ambiguous.allele.matches}: Corresponds to the \code{return.ambiguous.as.missing} argument in \code{assess.pgs.vcf.allele.match}. When \code{TRUE}, non-INDEL allele
#' mismatches that cannot be resolved (due to palindromic alleles or causes other than strand flips) are removed by marking the affected value in the \code{effect_allele} column as missing
#' prior to dosage calling and missing genotype handling. The corresponding dosage is set to NA and the variant is handled according to the chosen missing genotype method.
#' \item \code{remove.mismatched.indels}: Corresponds to the \code{return.indels.as.missing} argument in \code{assess.pgs.vcf.allele.match}. When \code{TRUE}, INDEL allele mismatches
#' (which cannot be assessed for strand flips) are removed by marking the affected value in the \code{effect_allele} column as missing prior to dosage calling and missing genotype handling.
#' The corresponding dosage is set to NA and the variant is handled according to the chosen missing genotype method.
#' }
#'
#' Note that an allele match assessment requires the presence of both the \code{other_allele} and \code{effect_allele} in the PGS weight data.
#' The \code{other_allele} column is not required by the PGS Catalog, and so is not always available.
#'
#' @examples
#' # Example VCF
#' vcf.path <- system.file(
#'     'extdata',
#'     'HG001_GIAB.vcf.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' vcf.import <- import.vcf(vcf.path);
#'
#' # Example pgs weight file
#' pgs.weight.path <- system.file(
#'     'extdata',
#'     'PGS000662_hmPOS_GRCh38.txt.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' pgs.import <- import.pgs.weight.file(pgs.weight.path);
#'
#' pgs.data <- apply.polygenic.score(
#'     vcf.data = vcf.import$dat,
#'     pgs.weight.data = pgs.import$pgs.weight.data,
#'     missing.genotype.method = 'none'
#'     );
#'
#' # Specify different methods for handling missing genotypes
#' pgs.import$pgs.weight.data$allelefrequency_effect <- rep(0.5, nrow(pgs.import$pgs.weight.data));
#' pgs.data <- apply.polygenic.score(
#'     vcf.data = vcf.import$dat,
#'     pgs.weight.data = pgs.import$pgs.weight.data,
#'     missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
#'     use.external.effect.allele.frequency = TRUE
#'     );
#'
#' # Specify allele mismatch handling
#' pgs.data <- apply.polygenic.score(
#'    vcf.data = vcf.import$dat,
#'    pgs.weight.data = pgs.import$pgs.weight.data,
#'    correct.strand.flips = TRUE,
#'    remove.ambiguous.allele.matches = TRUE,
#'    remove.mismatched.indels = FALSE
#'    );
#'
#' # Provide phenotype data for basic correlation analysis
#' phenotype.data <- data.frame(
#'     Indiv = unique(vcf.import$dat$Indiv),
#'     continuous.phenotype = rnorm(length(unique(vcf.import$dat$Indiv))),
#'     binary.phenotype = sample(
#'         c('a', 'b'),
#'         length(unique(vcf.import$dat$Indiv)),
#'         replace = TRUE
#'         )
#'     );
#'
#' pgs.data <- apply.polygenic.score(
#'     vcf.data = vcf.import$dat,
#'     pgs.weight.data = pgs.import$pgs.weight.data,
#'     phenotype.data = phenotype.data
#'     );
#'
#' # Only run validation checks on input data and report back
#' apply.polygenic.score(
#'     vcf.data = vcf.import$dat,
#'     pgs.weight.data = pgs.import$pgs.weight.data,
#'     validate.inputs.only = TRUE
#'     );
#' @export
apply.polygenic.score <- function(
    vcf.data,
    pgs.weight.data,
    phenotype.data = NULL,
    phenotype.analysis.columns = NULL,
    correct.strand.flips = TRUE,
    remove.ambiguous.allele.matches = FALSE,
    remove.mismatched.indels = FALSE,
    output.dir = NULL,
    file.prefix = NULL,
    missing.genotype.method = 'mean.dosage',
    use.external.effect.allele.frequency = FALSE,
    n.percentiles = NULL,
    analysis.source.pgs = NULL,
    validate.inputs.only = FALSE
    ) {

    ### Start Input Validation ###

    validate.vcf.input(vcf.data = vcf.data);
    validate.pgs.data.input(
        pgs.weight.data = pgs.weight.data,
        use.external.effect.allele.frequency = use.external.effect.allele.frequency,
        correct.strand.flips = correct.strand.flips,
        remove.ambiguous.allele.matches = remove.ambiguous.allele.matches,
        remove.mismatched.indels = remove.mismatched.indels
        );
    validate.phenotype.data.input(phenotype.data = phenotype.data, phenotype.analysis.columns = phenotype.analysis.columns, vcf.data = vcf.data);

    if (validate.inputs.only) {
        message('Input data passed validation');
        return(TRUE);
        }

    # check missing genotype method input
    if (all(missing.genotype.method %in% c('mean.dosage', 'normalize', 'none'))) {
        missing.genotype.method <- missing.genotype.method;
        } else {
        stop('missing.genotype.method must be either "mean.dosage", "normalize", or "none"');
        }

    # check that n.percentiles is a mathematical integer
    if (!is.null(n.percentiles) && (n.percentiles %% 1 != 0)) {
        stop('n.percentiles must be an integer');
        }

    # check that analysis.source.pgs is NULL or a character string representing a missing genotype method
    if (!is.null(analysis.source.pgs)) {
        if (length(analysis.source.pgs) > 1) {
            stop('analysis.source.pgs must be one of the chosen missing genotype methods');
            }
        if (!(analysis.source.pgs %in% missing.genotype.method)) {
            stop('analysis.source.pgs must be one of the chosen missing genotype methods');
            }
        } else {
        # if no source is provided, set analysis.source.pgs to the first chosen method.
        analysis.source.pgs <- missing.genotype.method[1];
        }

    ### End Input Validation ###

    # merge VCF and PGS data
    merged.vcf.with.pgs.data <- combine.vcf.with.pgs(
        vcf.data = vcf.data,
        pgs.weight.data = pgs.weight.data
        )$merged.vcf.with.pgs.data;

    # free up some memory
    rm(vcf.data);

    ### Start Allele Match Check ###
    if (remove.ambiguous.allele.matches || correct.strand.flips) {
        match.assessment <- ApplyPolygenicScore::assess.pgs.vcf.allele.match(
            vcf.ref.allele = merged.vcf.with.pgs.data$REF,
            vcf.alt.allele = merged.vcf.with.pgs.data$ALT,
            pgs.ref.allele = merged.vcf.with.pgs.data$other_allele,
            pgs.effect.allele = merged.vcf.with.pgs.data$effect_allele,
            return.ambiguous.as.missing = remove.ambiguous.allele.matches,
            return.indels.as.missing = remove.mismatched.indels
            );
        merged.vcf.with.pgs.data$effect_allele <- match.assessment$new.pgs.effect.allele;
        }

    # calculate dosage
    merged.vcf.with.pgs.data$dosage <- convert.alleles.to.pgs.dosage(
        called.alleles = merged.vcf.with.pgs.data$gt_GT_alleles,
        risk.alleles = merged.vcf.with.pgs.data$effect_allele
        );

    ### Start Missing Genotype Handling ###
    variant.id <- paste(merged.vcf.with.pgs.data$CHROM, merged.vcf.with.pgs.data$POS, merged.vcf.with.pgs.data$effect_allele, sep = ':');
    if ('mean.dosage' %in% missing.genotype.method) {
        # calculate dosage to replace missing genotypes
        if (use.external.effect.allele.frequency) {
            missing.genotype.dosage <- data.frame(dosage = convert.allele.frequency.to.dosage(allele.frequency = merged.vcf.with.pgs.data$allelefrequency_effect));
            missing.genotype.dosage$variant.id <- variant.id;
            # remove duplicated rows
            missing.genotype.dosage <- missing.genotype.dosage[!duplicated(missing.genotype.dosage$variant.id), ];
            # convert to named vector
            missing.genotype.dosage <- setNames(missing.genotype.dosage$dosage, missing.genotype.dosage$variant.id);
            # identify missing genotypes
            # this method includes variants that are missing in all Indivs, these are all replaced with the same mean dosage
            missing.genotype.row.index <- which(is.na(merged.vcf.with.pgs.data$dosage));
            } else {
            # create sample by variant dosage matrix

            dosage.matrix <- get.variant.by.sample.matrix(
                long.data = merged.vcf.with.pgs.data,
                variant.id = variant.id,
                value.var = 'dosage'
                );
            missing.genotype.dosage <- calculate.missing.genotype.dosage(dosage.matrix = dosage.matrix);
            # free up some memory
            rm(dosage.matrix);
            # identify missing genotypes
            # this method excludes variants that are missing in all Indivs, these remain as NA
            missing.genotype.row.index <- which(is.na(merged.vcf.with.pgs.data$dosage) & !is.na(merged.vcf.with.pgs.data$Indiv));
            }


        # start a column for replaced missing dosages
        merged.vcf.with.pgs.data$dosage.with.replaced.missing <- merged.vcf.with.pgs.data$dosage;
        # assign mean dosage to missing genotypes
        for (i in missing.genotype.row.index) {
            missing.variant.id <- paste(merged.vcf.with.pgs.data[i, 'CHROM'], merged.vcf.with.pgs.data[i, 'POS'], merged.vcf.with.pgs.data[i, 'effect_allele'], sep = ':');
            missing.variant.dosage <- missing.genotype.dosage[missing.variant.id];
            merged.vcf.with.pgs.data[i, 'dosage.with.replaced.missing'] <- missing.variant.dosage;
            }
    }

    ### End Missing Genotype Handling ###

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
        # skip if the variant is missing from all samples
        if (is.na(merged.vcf.with.pgs.data[i, 'Indiv'])) {
            next;
            }
        key <- paste(merged.vcf.with.pgs.data[i, 'Indiv'], merged.vcf.with.pgs.data[i, 'CHROM'], merged.vcf.with.pgs.data[i, 'POS'], sep = '_');
        # save all row indexes that have the same sample:coordinate combination under one key
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

    ### Start Missing SNP Count ###
    bialellic.variant.id <- paste(merged.vcf.with.pgs.data$CHROM, merged.vcf.with.pgs.data$POS, sep = ':');
    biallelic.snp.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
        long.data = merged.vcf.with.pgs.data,
        variant.id = bialellic.variant.id,
        value.var = 'dosage'
        );
    per.sample.missing.genotype.count <- colSums(is.na(biallelic.snp.by.sample.matrix));

    per.sample.missing.genotype.percent <- round(per.sample.missing.genotype.count / nrow(biallelic.snp.by.sample.matrix), 2);

    ### End Missing SNP Count ###

    ### Start PGS Application ###
    # calculate PGS per sample
    pgs.output.list <- list();

    missing.method.to.colname.ref <- c(
        'mean.dosage' = 'PGS.with.replaced.missing',
        'normalize' = 'PGS.with.normalized.missing',
        'none' = 'PGS'
        );

    if ('none' %in% missing.genotype.method) {
        pgs.per.sample <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample) <- c('Indiv', 'PGS');
        pgs.output.list$PGS <- pgs.per.sample;
        # free up some memory
        rm(pgs.per.sample);
        }

    if ('normalize' %in% missing.genotype.method) {
        pgs.per.sample.with.normalized.missing <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample.with.normalized.missing) <- c('Indiv', 'PGS');
        per.sample.non.missing.genotype.count <- colSums(!is.na(biallelic.snp.by.sample.matrix));
        ploidy <- 2; # hard-coded ploidy for human diploid genome
        per.sample.non.missing.genotype.count.ploidy.adjusted <- ploidy * per.sample.non.missing.genotype.count
        pgs.per.sample.with.normalized.missing$PGS <- pgs.per.sample.with.normalized.missing$PGS / per.sample.non.missing.genotype.count.ploidy.adjusted;
        # account for division by zero
        pgs.per.sample.with.normalized.missing$PGS[is.nan(pgs.per.sample.with.normalized.missing$PGS)] <- NA;
        pgs.output.list$PGS.with.normalized.missing <- pgs.per.sample.with.normalized.missing;
        # free up some memory
        rm(pgs.per.sample.with.normalized.missing);
        }

    if ('mean.dosage' %in% missing.genotype.method) {
        pgs.per.sample <- aggregate(
            x = merged.vcf.with.pgs.data$multiallelic.weighted.dosage.with.replaced.missing,
            by = list(merged.vcf.with.pgs.data$Indiv),
            FUN = sum,
            na.rm = TRUE
            );
        colnames(pgs.per.sample) <- c('Indiv', 'PGS');
        pgs.output.list$PGS.with.replaced.missing <- pgs.per.sample;
        # free up some memory
        rm(pgs.per.sample);
        }

    ### End PGS Application ###

    # format output
    # bind PGS columns of list components
    PGS.cols <- lapply(pgs.output.list, function(x) x$PGS);
    PGS.cols <- data.frame(do.call(cbind, PGS.cols));
    colnames(PGS.cols) <- names(pgs.output.list);
    # bind sample column of first list component
    Indiv <- pgs.output.list[[1]]$Indiv;
    pgs.output <- cbind(Indiv, PGS.cols);

    # retrieve pgs for statisitical analyses
    pgs.for.stats <- pgs.output[ , missing.method.to.colname.ref[analysis.source.pgs]];

    # calculate percentiles
    percentiles <- get.pgs.percentiles(pgs = pgs.for.stats, n.percentiles = n.percentiles);

    pgs.output <- cbind(pgs.output, percentiles);

    # add missing genotype count
    pgs.output$n.missing.genotypes <- per.sample.missing.genotype.count;
    pgs.output$percent.missing.genotypes <- per.sample.missing.genotype.percent;

    # initialize regression output
    regression.output <- NULL;

    # merge PGS data with phenotype data by Indiv column
    if (!is.null(phenotype.data)) {
        pgs.output <- merge(
            x = pgs.output,
            y = phenotype.data,
            by = 'Indiv',
            all.x = TRUE,
            all.y = TRUE
            );

        ### Begin Phenotype Analysis ###

        if (!is.null(phenotype.analysis.columns)) {
            regression.output <- run.pgs.regression(
                pgs = pgs.for.stats,
                phenotype.data = subset(phenotype.data, select = phenotype.analysis.columns)
                );
            }
        ### End Phenotype Analysis ###

        }

    final.output <- list(
        pgs.output = pgs.output,
        regression.output = regression.output
        );

    # If requested, write output to file
    if (!is.null(output.dir)) {

        if (is.null(file.prefix)) {
            file.prefix <- 'ApplyPolygenicScore';
            }

        write.apply.polygenic.score.output.to.file(
            apply.polygenic.score.output = final.output,
            output.dir = output.dir,
            file.prefix = file.prefix
            );
        }

    return(final.output);
    }
