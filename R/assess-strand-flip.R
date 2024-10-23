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
    snp.alleles <- alleles[nchar(alleles) == 1];
    accepted.alleles <- c('A', 'T', 'C', 'G', '*');
    allele.check <- snp.alleles[!is.na(snp.alleles)] %in% accepted.alleles;

    if (any(!allele.check, na.rm = TRUE)) {
        stop('Invalid allele: ', snp.alleles[!is.na(snp.alleles)][!allele.check]);
        }
    }

#' @title Flip DNA allele
#' @description Flip single base pair DNA alleles to their reverse complement. INDEL flipping is not supported.
#' @param alleles A character vector of DNA alleles.
#' @param return.indels.as.missing A logical value indicating whether to return NA for INDEL alleles. Default is FALSE.
#' @return A character vector of flipped DNA alleles. INDEL alleles are returned as is unless return.indels.as.missing is TRUE.
#' @examples
#' alleles <- c('A', 'T', 'C', 'G', ATG, NA);
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

assess.strand.flip <- function(
    vcf.ref.allele,
    vcf.alt.allele,
    pgs.ref.allele,
    pgs.effect.allele,
    return.indels.as.missing = FALSE,
    return.ambiguous.as.missing = FALSE
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

    # verify valid alleles
    validate.allele.input(vcf.ref.allele, na.allowed = FALSE);
    validate.allele.input(pgs.ref.allele, na.allowed = FALSE);
    validate.allele.input(pgs.effect.allele, na.allowed = FALSE);

    # split multiallelic sites in VCF ALT field and validate
    vcf.alt.allele.split <- sapply(
        X = vcf.alt.allele,
        FUN = function(x) {
            if (is.na(x)) {
                return(NA);
                } else {
                    return(unlist(strsplit(x, ',')));
                }
            }
        );
    validate.allele.input(vcf.alt.allele.split, na.allowed = FALSE);

    # initialize empty vectors to store flip info
    flipped.effect.allele <- rep(NA, length(vcf.ref.allele));
    flipped.other.allele <- rep(NA, length(vcf.ref.allele));
    # strand status
    flip.designation <- rep(NA, length(vcf.ref.allele));

    for (i in 1:length(vcf.ref.allele)) {

        current.pgs.ref.allele <- pgs.ref.allele[i];
        current.pgs.effect.allele <- pgs.effect.allele[i];
        current.vcf.ref.allele <- vcf.ref.allele[i];

        # VCF ALT field may contain multiallelic variants
        # Multiple alleles are separated by commas e.g. 'A,T'
        # Split the VCF ALT field into a vector of alleles
        current.vcf.alt.allele.split <- unlist(strsplit(vcf.alt.allele[i], ','));

        # check if default ref-ref / alt-effect alleles match
        default.ref.check <- current.pgs.ref.allele == current.vcf.ref.allele;
        default.alt.check <- any(current.vcf.alt.allele.split %in% current.pgs.effect.allele);

        if (default.ref.check && default.alt.check) {
            flip.designation[i] <- 'default_match';
            flipped.effect.allele[i] <- pgs.effect.allele[i];
            flipped.other.allele[i] <- pgs.ref.allele[i];
            break;
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
            # no INDEL flipping supported, return either NA or the original allele
            warning('Mismatch detected in INDEL PGS allele. Skipping strand flip assessment');
            flip.designation[i] <- 'indel_mismatch';
            if (return.indels.as.missing) {
                flipped.effect.allele[i] <- NA;
                flipped.other.allele[i] <- NA;
                break;
                } else {
                    flipped.effect.allele[i] <- pgs.effect.allele[i];
                    flipped.other.allele[i] <- pgs.ref.allele[i];
                    break;
                    }
            } else {
                pgs.indel <- FALSE;
            }

        # identify insertion/deletion alleles in VCF
        if (nchar(current.vcf.ref.allele) > 1 || all(nchar(current.vcf.alt.allele.split) > 1)) {
            vcf.indel <- TRUE;
            # no INDEL flipping supported, return either NA or the original allele
            warning('Mismatch detected in INDEL VCF allele. Skipping strand flip assessment');
            flip.designation[i] <- 'indel_mismatch';
            if (return.indels.as.missing) {
                flipped.effect.allele[i] <- NA;
                flipped.other.allele[i] <- NA;
                break;
                } else {
                    flipped.effect.allele[i] <- current.pgs.effect.allele;
                    flipped.other.allele[i] <- current.pgs.ref.allele;
                    break;
                    }
            } else {
                vcf.indel <- FALSE;
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
            flip.designation[i] <- 'ambiguous_flip';
            if (return.ambiguous.as.missing) {
                flipped.effect.allele[i] <- NA;
                flipped.other.allele[i] <- NA;
                break;
                } else {
                flipped.effect.allele[i] <- current.pgs.effect.allele;
                flipped.other.allele[i] <- current.pgs.ref.allele;
                break;
                }
            } else if (effect.switch.candidate) {
            # if this is a clear-cut effect switch, return the default PGS alleles
            # apply.polygenic.score automatically handles effect switches during PGS application
            flip.designation[i] <- 'effect_switch';
            flipped.effect.allele[i] <- current.pgs.effect.allele;
            flipped.other.allele[i] <- current.pgs.ref.allele;
            break;
            } else if (strand.flip.candidate) {
            # if this is a clear-cut strand flip, return the flipped PGS alleles
            flip.designation[i] <- 'strand_flip';
            flipped.effect.allele[i] <- pgs.effect.flip;
            flipped.other.allele[i] <- pgs.ref.flip;
            break;
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
            break;
            } else {
                # no solution found for ref/ref alt/alt mismatch
                # same return strategy as ambiguous case
                flip.designation[i] <- 'unresolved_mismatch';
                if (return.ambiguous.as.missing) {
                    flipped.effect.allele[i] <- NA;
                    flipped.other.allele[i] <- NA;
                    break;
                    } else {
                    flipped.effect.allele[i] <- current.pgs.effect.allele;
                    flipped.other.allele[i] <- current.pgs.ref.allele;
                    break;
                    }
                }

        }

    return(
        list(
            match.status = flip.designation,
            new.pgs.effect.allele = flipped.effect.allele,
            new.pgs.other.allele = flipped.other.allele
            )
        );

    }

# assess.strand.flip(
#     vcf.ref.allele = 'A',
#     vcf.alt.allele = 'G',
#     pgs.ref.allele = 'C',
#     pgs.effect.allele = 'TAA',
#     return.indels.as.missing = FALSE,
#     return.ambiguous.as.missing = FALSE
#     )
