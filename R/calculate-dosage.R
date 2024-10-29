#' @title Convert alleles to dosage
#' @description Convert genotype calls in the form of witten out alleles (e.g. 'A/T') to dosages (0, 1, 2) based on provided risk alleles from a PGS.
#' @param called.alleles A vector of genotypes in allelic notation separated by a slash or pipe.
#' @param risk.alleles A vector of risk alleles from a polygenic score corresponding to each genotype (by locus) in called.alleles.
#' @return A vector of dosages corresponding to each genotype in called.alleles.
#' @examples
#' called.alleles <- c('A/A', 'A/T', 'T/T');
#' risk.alleles <- c('T', 'T', 'T');
#' convert.alleles.to.pgs.dosage(called.alleles, risk.alleles);
#' @export
convert.alleles.to.pgs.dosage <- function(called.alleles, risk.alleles) {
    # check that risk.alleles is the same length as called.alleles
    if (length(called.alleles) != length(risk.alleles)) {
        stop('called.alleles and risk.alleles must be the same length.');
        }

    # check for missing risk alleles and warn
    if (any(is.na(risk.alleles))) {
        warning('risk.alleles contains missing values, NA will be returned for corresponding dosage.');
        }

    # check that risk.alleles is a vector of capitalized alphabetic characters
    if (!all(grepl('^[A-Z]+$', risk.alleles[!is.na(risk.alleles)]))) { # NA risk alleles allowed
        stop('unrecognized risk.allele format, must be capitalized letters.');
        }

    # handle totally missing genotypes
    # if the entire vector is NA or the entire vector is '.', return NA
    if (all(is.na(called.alleles)) | all(called.alleles == '.')) {
        split.alleles <- data.frame(called.alleles, called.alleles);
        } else {
            # check that called.alleles is a vector of genotypes in allelic notation or '.' separated by a slash or pipe
            # "*" characters represent overlapping deletions from an upstream indel and are accepted VCF format
            allowed.pattern <- '^((([A-Z]+|\\.|\\*)[/\\|]([A-Z]+|\\.|\\*))|\\.)$' # '|' are special chars in regular expressions
            passing.alleles <- grepl(allowed.pattern, called.alleles);
            passing.alleles[is.na(called.alleles)] <- TRUE; # NA allowed
            if (!all(passing.alleles)) {
                stop('unrecognized called.alleles format, must be capitalized letters, "." or "*" separated by a slash or pipe.');
                }
            split.alleles <- data.table::tstrsplit(called.alleles, split = c('/|\\|'), keep = c(1,2)); # '|' are special chars in regular expressions
            }
    names(split.alleles) <- c('called.allele.a', 'called.allele.b');

    # replace 'NA' with '.' for easier comparisons
    missing.label <- '.';
    split.alleles <- lapply(
        X = split.alleles,
        FUN = function(x) {
        x[is.na(x)] <- missing.label;
        return(x);
        }
        );

    dosage <- rep(NA, length(called.alleles));
    for (i in 1:length(called.alleles)) {
        if (is.na(risk.alleles[i])) {
            dosage[i] <- NA; # if the risk allele is missing, return NA, no dosage can be calculated
            } else if ((split.alleles$called.allele.a[i] == missing.label) & (split.alleles$called.allele.b[i] == missing.label)) {
            dosage[i] <- NA; # if both allelles are missing, no genotype was called, return NA
            } else if (split.alleles$called.allele.a[i] == missing.label | split.alleles$called.allele.b[i] == missing.label) {
                dosage[i] <- NA; # if one of the alleles is marked as missing but the other is not, this is an unrecognized format
                warning('one of two alleles is marked as missing at index ', i, ', this is an unrecognized format, returning NA for dosage.');
            } else if (split.alleles$called.allele.a[i] == risk.alleles[i] & split.alleles$called.allele.b[i] == risk.alleles[i]) {
                dosage[i] <- 2; # if both alleles are the risk allele, the genotype is homozygous for the effect allele and the dosage is 2.
            } else if (split.alleles$called.allele.a[i] == risk.alleles[i] | split.alleles$called.allele.b[i] == risk.alleles[i]) {
                dosage[i] <- 1; # if only one of the alleles is the risk allele, the genotype is heterozygous and the dosage is 1.
            } else {
                dosage[i] <- 0; # if neither allele is the risk allele, the genotype is homozygous for the non-effect allele and the dosage is 0.
            }
        }
    return(dosage);
    }

# The function for calculating a dosage value intended to replace missing genotypes.
calculate.missing.genotype.dosage <- function(dosage.matrix) {
    # calculate the mean dosage for each variant
    mean.dosage <- apply(
        X = dosage.matrix,
        MARGIN = 1,
        FUN = function(x) {
            # simple mean
            mean(x, na.rm = TRUE)
            }
        );
    return(mean.dosage);
    }

#' @title Convert allele frequency to mean dosage
#' @description Convert a population allele frequency to a mean dosage for that allele.
#' @param allele.frequency A numeric vector of allele frequencies.
#' @return A numeric vector of mean dosages for the allele frequencies.
#' @examples
#' allele.frequency <- seq(0.1, 0.9, 0.1);
#' convert.allele.frequency.to.dosage(allele.frequency);
#' @export
convert.allele.frequency.to.dosage <- function(allele.frequency) {
    # check that allele.frequency is a numeric vector
    if (!is.numeric(allele.frequency)) {
        stop('allele.frequency must be a numeric vector.');
        }
    # check for NA
    if (any(is.na(allele.frequency))) {
        stop('allele.frequency must not contain NA.');
        }
    # check that allele.frequency is between 0 and 1
    if (any(allele.frequency < 0) | any(allele.frequency > 1)) {
        stop('allele.frequency must be between 0 and 1.');
        }
    # calculate dosage
    dosage <- 2 * allele.frequency;
    return(dosage);
    }
