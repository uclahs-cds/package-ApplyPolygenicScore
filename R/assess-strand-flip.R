validate.allele.input <- function(alleles, na.allowed = TRUE) {
    # check that alleles is a character vector
    if (!is.character(alleles)) {
        stop('allele must be a character vector');
        }

    if (!na.allowed) {
        # check for NA
        if (any(is.na(alleles))) {
            stop('allele must not contain NA');
            }
        }

    # Verify acceptable SNP alleles
    is.snp <- nchar(alleles) == 1;
    snp.alleles <- alleles[is.snp & !is.na(alleles)];
    # convert to uppercase
    snp.alleles.upper <- toupper(snp.alleles);
    accepted.alleles <- c('A', 'T', 'C', 'G', '*'); # Note that "*" is a regex wild card and may not be comparable if using regular expressions
    allele.check <- snp.alleles.upper %in% accepted.alleles;

    if (any(!allele.check, na.rm = TRUE)) {
        invalid.snps <- snp.alleles[!allele.check];
        stop('Invalid allele: ', invalid.snps);
        }
    }

#' @title Flip DNA allele
#' @description Flip single base pair DNA alleles to their reverse complement. INDEL flipping is not supported.
#' @param alleles A character vector of DNA alleles.
#' @param return.indels.as.missing A logical value indicating whether to return NA for INDEL alleles. Default is \code{FALSE}.
#' @return A character vector of flipped DNA alleles. INDEL alleles are returned as is unless \code{return.indels.as.missing} is \code{TRUE}.
#' @examples
#' alleles <- c('A', 'T', 'C', 'G', 'ATG', NA);
#' flip.DNA.allele(alleles);
#' @export
flip.DNA.allele <- function(alleles, return.indels.as.missing = FALSE) {
    if (all(is.na(alleles))) {
        return(NA);
        }
    # validate allele input
    validate.allele.input(alleles);

    flipped.alleles <- sapply(
        X = alleles,
        FUN = function(x) {
            # NA handling
            if (is.na(x)) {
                return(NA);
                # INDEL handling
                } else if (nchar(x) > 1) {
                if (return.indels.as.missing) {
                    return(NA);
                    } else {
                        return(x); # no INDEL flipping supported
                        }
                # SNP handling
                } else {
                    flip <- switch(
                        x,
                        'A' = 'T',
                        'T' = 'A',
                        'C' = 'G',
                        'G' = 'C',
                        '*' = '*', # wildcard for missing alleles
                        `NA` = NA
                        );
                }
            }
        );

    names(flipped.alleles) <- NULL;
    return(flipped.alleles);
    }
#' @title Assess PGS allele match to VCF allele
#' @description Assess whether PGS reference and effect alleles match provided VCF reference and alternative alleles.
#' Mismatches are checked for potential switching of effect and reference PGS alleles (cases where the effect allele is the REF VCF allele)
#' and are evaluated for DNA strand flips (by flipping the PGS alleles). INDEL alleles are not supported for strand flip assessment.
#' @param vcf.ref.allele A character vector of singular VCF reference (REF) alleles.
#' @param vcf.alt.allele A character vector of VCF alternative (ALT) alleles. Multiple alleles at a multiallelic site must be separated by commas.
#' @param pgs.ref.allele A character vector of singular PGS reference alleles aka "non-effect" or "other" alleles.
#' @param pgs.effect.allele A character vector of singular PGS effect alleles.
#' @param return.indels.as.missing A logical value indicating whether to return NA for INDEL alleles with detected mismatches. Default is \code{FALSE}.
#' @param return.ambiguous.as.missing A logical value indicating whether to return NA for ambiguous cases where both a strand flip and effect switch are possible,
#' or no strand flip is detected and a mismatch cannot be resolved. Default is \code{FALSE}.
#' @param max.strand.flips An integer indicating the number of non-ambiguous strand flips that must be present to implement the discarding all allele matches labeled "ambiguous_flip". Only applies if {return.ambiguous.as.missing == TRUE}.
#' Defaults to \code{0}, meaning that no strand flips are allowed. Allele matches labeled "unresolved_mismatch" are not affected by this parameter.
#' @return A list containing the match assessment, a new PGS effect allele, and a new PGS other allele.
#'
#' \strong{Output Structure}
#'
#' The outputed list contains the following elements:
#' \itemize{
#' \item \code{match.status}: A character vector indicating the match status for each pair of allele pairs. Possible values are \code{default_match}, \code{effect_switch}, \code{strand_flip}, \code{effect_switch_with_strand_flip}, \code{ambiguous_flip}, \code{indel_mismatch}, and \code{unresolved_mismatch}.
#' \item \code{new.pgs.effect.allele}: A character vector of new PGS effect alleles based on the match status. If the match status is \code{default_match}, \code{effect_switch} or \code{missing_allele}, the original PGS effect allele is returned.
#' If the match status is \code{strand_flip} or \code{effect_switch_with_strand_flip} the flipped PGS effect allele is returned. If the match status is \code{ambiguous_flip}, \code{indel_mismatch}, or \code{unresolved_mismatch},
#' the return value is either the original allele or NA as dictated by the \code{return.indels.as.missing}, \code{return.ambiguous.as.missing}, and \code{max.strand.flips} parameters.
#' \item \code{new.pgs.other.allele}: A character vector of new PGS other alleles based on the match status, following the same logic as \code{new.pgs.effect.allele}.
#' }
#'
#' The match.status output indicates the following:
#' \itemize{
#' \item \code{default_match}: The default PGS reference allele matches the VCF REF allele and the default PGS effect allele matches one of the VCF ALT alleles.
#' \item \code{effect_switch}: The PGS effect allele matches the VCF REF allele and the PGS reference allele matches one of the VCF ALT alleles.
#' \item \code{strand_flip}: The PGS reference and effect alleles match their respective VCF pairs when flipped.
#' \item \code{effect_switch_with_strand_flip}: The PGS effect allele matches the VCF REF allele and the PGS reference allele matches one of the VCF ALT alleles when flipped.
#' \item \code{ambiguous_flip}: Both an effect switch and a strand flip have been detected. This is an ambiguous case caused by palindromic SNPs.
#' \item \code{indel_mismatch}: A mismatch was detected between pairs of alleles where at least one was an INDEL. INDEL alleles are not supported for strand flip assessment.
#' \item \code{unresolved_mismatch}: A mismatch was detected between pairs of non-INDEL alleles that could not be resolved by an effect switch or flipping the PGS alleles.
#' \item \code{missing_allele}: One of the four alleles is missing, making it impossible to assess the match.
#' }
#' @examples
#' # Example data demonstrating the following cases in each vector element:
#' # 1. no strand flips
#' # 2. effect allele switch
#' # 3. strand flip
#' # 4. effect allele switch AND strand flip
#' # 5. palindromic (ambiguous) alleles
#' # 6. unresolved mismatch
#' vcf.ref.allele <- c('A', 'A', 'A', 'A', 'A', 'A');
#' vcf.alt.allele <- c('G', 'G', 'G', 'G', 'T', 'G');
#' pgs.ref.allele <- c('A', 'G', 'T', 'C', 'T', 'A');
#' pgs.effect.allele <- c('G', 'A', 'C', 'T', 'A', 'C');
#' assess.pgs.vcf.allele.match(vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, pgs.effect.allele);
#' @export
assess.pgs.vcf.allele.match <- function(
    vcf.ref.allele,
    vcf.alt.allele,
    pgs.ref.allele,
    pgs.effect.allele,
    return.indels.as.missing = FALSE,
    return.ambiguous.as.missing = FALSE,
    max.strand.flips = 0
    ) {

    # check that all inputs are one dimensional character vectors
    if (!is.character(vcf.ref.allele) | !is.character(vcf.alt.allele) | !is.character(pgs.ref.allele) | !is.character(pgs.effect.allele)) {
        stop('alleles must be character vectors');
        }

    # check that all inputs are the same length
    if (length(vcf.ref.allele) != length(vcf.alt.allele)) {
        stop('vcf.ref.allele and vcf.alt.allele must be the same length.');
        }

    if (length(pgs.ref.allele) != length(pgs.effect.allele)) {
        stop('pgs.ref.allele and pgs.effect.allele must be the same length.');
        }

    if (length(vcf.ref.allele) != length(pgs.ref.allele)) {
        stop('vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, and pgs.effect.allele must be the same length.');
        }

    # check that all logical parameters are TRUE or FALSE
    if (!is.logical(return.indels.as.missing) || length(return.indels.as.missing) != 1) {
        stop('return.indels.as.missing must be a single logical value.');
        }

    if (!is.logical(return.ambiguous.as.missing) || length(return.ambiguous.as.missing) != 1) {
        stop('return.ambiguous.as.missing must be a single logical value.');
        }

    # check that max.strand.flips is a single integer
    if (!is.numeric(max.strand.flips) || length(max.strand.flips) != 1 || max.strand.flips < 0 || max.strand.flips != floor(max.strand.flips)) {
        stop('max.strand.flips must be a single non-negative integer.');
        }

    # verify valid alleles
    validate.allele.input(vcf.ref.allele, na.allowed = TRUE);
    validate.allele.input(pgs.ref.allele, na.allowed = TRUE);
    validate.allele.input(pgs.effect.allele, na.allowed = TRUE);

    # split multiallelic sites in VCF ALT field and validate
    vcf.alt.allele.split <- unlist(sapply(
        X = vcf.alt.allele,
        FUN = function(x) {
            if (is.na(x)) {
                return(NA);
                } else {
                    return(unlist(strsplit(x, ',')));
                }
            }
        ));
    validate.allele.input(vcf.alt.allele.split, na.allowed = TRUE);

    # initialize empty vectors to store flip info
    flipped.effect.allele <- rep(NA, length(vcf.ref.allele));
    flipped.other.allele <- rep(NA, length(vcf.ref.allele));
    # strand status
    flip.designation <- rep(NA, length(vcf.ref.allele));

    for (i in seq_along(vcf.ref.allele)) {

        current.pgs.ref.allele <- pgs.ref.allele[i];
        current.pgs.effect.allele <- pgs.effect.allele[i];
        current.vcf.ref.allele <- vcf.ref.allele[i];

        # VCF ALT field may contain multiallelic variants
        # Multiple alleles are separated by commas e.g. 'A,T'
        # Split the VCF ALT field into a vector of alleles
        current.vcf.alt.allele.split <- unlist(strsplit(vcf.alt.allele[i], ','));

        # NA handling.
        # If one of the four alleles is missing, it is not possible to assess the match.
        # Note that A VCF can have a missing ALT allele when in GVCF format and reporting
        # a site where all individuals are non-variant (homozygous reference)
        if (all(is.na(current.vcf.alt.allele.split)) || is.na(current.vcf.ref.allele) || is.na(current.pgs.ref.allele) || is.na(current.pgs.effect.allele)) {
            flip.designation[i] <- 'missing_allele';
                flipped.effect.allele[i] <- current.pgs.effect.allele;
                flipped.other.allele[i] <- current.pgs.ref.allele;
                next;
                }

        # check if default ref-ref / alt-effect alleles match
        default.ref.check <- current.pgs.ref.allele == current.vcf.ref.allele;
        default.alt.check <- any(current.vcf.alt.allele.split %in% current.pgs.effect.allele);

        if (default.ref.check && default.alt.check) {
            flip.designation[i] <- 'default_match';
            flipped.effect.allele[i] <- pgs.effect.allele[i];
            flipped.other.allele[i] <- pgs.ref.allele[i];
            next;
            }

        # check if PGS effect designation is on the VCF reference allele (effect switch)
        effect.switch.ref.check <- any(vcf.alt.allele.split %in% current.pgs.ref.allele);
        effect.switch.alt.check <- current.pgs.effect.allele == current.vcf.ref.allele;

        if (effect.switch.ref.check && effect.switch.alt.check) {
            # effect switch criteria is met, but could also be a palindromic strand flip
            effect.switch.candidate <- TRUE;
            } else {
                effect.switch.candidate <- FALSE;
            }


        # identify insertion/deletion alleles in PGS
        if (nchar(current.pgs.ref.allele) > 1 || nchar(current.pgs.effect.allele) > 1) {
            pgs.indel <- TRUE;
            if (effect.switch.candidate) {
                # if an INDEL effect switch is detected, do not continue to strand flip assessment
                # and return effect switch designation
                flip.designation[i] <- 'effect_switch';
                flipped.effect.allele[i] <- current.pgs.effect.allele;
                flipped.other.allele[i] <- current.pgs.ref.allele;
                next;
                }
            # no INDEL flipping supported, return either NA or the original allele
            warning('Mismatch detected in INDEL PGS allele. Skipping strand flip assessment.');
            flip.designation[i] <- 'indel_mismatch';
            if (return.indels.as.missing) {
                flipped.effect.allele[i] <- NA;
                flipped.other.allele[i] <- NA;
                next;
                } else {
                    flipped.effect.allele[i] <- pgs.effect.allele[i];
                    flipped.other.allele[i] <- pgs.ref.allele[i];
                    next;
                    }
            } else {
                pgs.indel <- FALSE;
            }

        # identify insertion/deletion alleles in VCF
        if (nchar(current.vcf.ref.allele) > 1 || all(nchar(current.vcf.alt.allele.split) > 1)) {
            if (effect.switch.candidate) {
                # if an INDEL effect switch is detected, do not continue to strand flip assessment
                # and return effect switch designation
                flip.designation[i] <- 'effect_switch';
                flipped.effect.allele[i] <- current.pgs.effect.allele;
                flipped.other.allele[i] <- current.pgs.ref.allele;
                next;
                }
            # no INDEL flipping supported, return either NA or the original allele
            warning('Mismatch detected in INDEL VCF allele. Skipping strand flip assessment.');
            flip.designation[i] <- 'indel_mismatch';
            if (return.indels.as.missing) {
                flipped.effect.allele[i] <- NA;
                flipped.other.allele[i] <- NA;
                next;
                } else {
                    flipped.effect.allele[i] <- current.pgs.effect.allele;
                    flipped.other.allele[i] <- current.pgs.ref.allele;
                    next;
                    }
            }


        # perform strand flip on SNP PGS alleles
        pgs.ref.flip <- flip.DNA.allele(pgs.ref.allele[i]);
        pgs.effect.flip <- flip.DNA.allele(pgs.effect.allele[i]);

        default.ref.flip.check <- current.vcf.ref.allele == pgs.ref.flip;
        default.alt.flip.check <- any(current.vcf.alt.allele.split %in% pgs.effect.flip);

        if (default.ref.flip.check && default.alt.flip.check) {
            # strand flip criteria is met, but could also be a palindromic effect switch
            strand.flip.candidate <- TRUE;
            } else {
                strand.flip.candidate <- FALSE;
            }

        # assess strand flip and effect switch candidates
        if (effect.switch.candidate && strand.flip.candidate) {
            # not possible to determine whether effect switch or strand flip has occured
            # This is an ambiguous case caused by palindromic SNPs
            # Default is to assume that this is an effect switch until strand flip threshold is exceeded
            flip.designation[i] <- 'ambiguous_flip';
            flipped.effect.allele[i] <- current.pgs.effect.allele;
            flipped.other.allele[i] <- current.pgs.ref.allele;
            next;
            } else if (effect.switch.candidate) {
            # if this is a clear-cut effect switch, return the default PGS alleles
            # apply.polygenic.score automatically handles effect switches during PGS application
            flip.designation[i] <- 'effect_switch';
            flipped.effect.allele[i] <- current.pgs.effect.allele;
            flipped.other.allele[i] <- current.pgs.ref.allele;
            next;
            } else if (strand.flip.candidate) {
            # if this is a clear-cut strand flip, return the flipped PGS alleles
            flip.designation[i] <- 'strand_flip';
            flipped.effect.allele[i] <- pgs.effect.flip;
            flipped.other.allele[i] <- pgs.ref.flip;
            next;
            }

        # The only remaining case is when a ref/ref and alt/alt mismatch was detected
        # but no strand flip or effect switch was detected.
        # One more resolvable case is when there is both an effect switch AND a strand flip

        # attempt strand flip on effect_switch ref-effect / alt-other alleles
        effect.switch.ref.flip.check <- current.vcf.ref.allele == pgs.effect.flip;
        effect.switch.alt.flip.check <- any(current.vcf.alt.allele.split %in% pgs.ref.flip);

        if (effect.switch.ref.flip.check && effect.switch.alt.flip.check) {
            flip.designation[i] <- 'effect_switch_with_strand_flip';
            flipped.effect.allele[i] <- pgs.effect.flip;
            flipped.other.allele[i] <- pgs.ref.flip;
            next;
            } else {
                # no solution found for ref/ref alt/alt mismatch
                # same return strategy as ambiguous case
                flip.designation[i] <- 'unresolved_mismatch';
                if (return.ambiguous.as.missing) {
                    flipped.effect.allele[i] <- NA;
                    flipped.other.allele[i] <- NA;
                    next;
                    } else {
                    flipped.effect.allele[i] <- current.pgs.effect.allele;
                    flipped.other.allele[i] <- current.pgs.ref.allele;
                    next;
                    }
                }

        }

    # implement max strand flip threshold for ambiguous cases
    # count number of unambiguous strand flips
    if (return.ambiguous.as.missing) {
        if (max.strand.flips > 0) {
            # if a max strand flips threshold is set, check if it is exceeded
            unambiguous.strand.flips <- sum(flip.designation == 'strand_flip' | flip.designation == 'effect_switch_with_strand_flip');
            if (unambiguous.strand.flips >= max.strand.flips) {
                # if the number of unambiguous strand flips equals or exceeds the threshold, return NA for all ambiguous cases
                flipped.effect.allele[flip.designation == 'ambiguous_flip'] <- NA;
                flipped.other.allele[flip.designation == 'ambiguous_flip'] <- NA;
                }
                # if max threshold is not exceeded, return all ambiguous cases with no flips (default)
            } else {
            # if max.strand.flips is 0 there is no tolerance for ambiguous cases, return all as NA
            flipped.effect.allele[flip.designation == 'ambiguous_flip'] <- NA;
            flipped.other.allele[flip.designation == 'ambiguous_flip'] <- NA;
            }
        # If return.ambiguous.as.missing is FALSE, keep all ambiguous cases as is (no flips)
        }

    return(
        list(
            match.status = flip.designation,
            new.pgs.effect.allele = flipped.effect.allele,
            new.pgs.other.allele = flipped.other.allele
            )
        );

    }
