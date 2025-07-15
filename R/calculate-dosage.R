# Handling CRAN warnings for data.table syntax:
if (getRversion() >= '2.15.1') utils::globalVariables(c('dosage'));

#' @title Convert alleles to dosage
#' @description Convert genotype calls in the form of witten out alleles (e.g. 'A/T') to dosages (0, 1, 2) based on provided risk alleles from a PGS.
#' @param called.alleles A vector of genotypes in allelic notation separated by a slash or pipe.
#' @param risk.alleles A vector of risk alleles from a polygenic score corresponding to each genotype (by locus) in called.alleles.
#' @return A vector of dosages corresponding to each genotype in called.alleles. Hemizygous genotypes (one allele e.g. 'A') are counted as 1.
#' @examples
#' called.alleles <- c('A/A', 'A/T', 'T/T');
#' risk.alleles <- c('T', 'T', 'T');
#' convert.alleles.to.pgs.dosage(called.alleles, risk.alleles);
#' @export
convert.alleles.to.pgs.dosage <- function(called.alleles, risk.alleles) {

    # Check input class and convert to a matrix for consistent processing
    is.vector.input <- is.vector(called.alleles);
    if (is.vector.input) {
        # Fast-fail for all missing genotypes
        if (all(is.na(called.alleles)) | all(called.alleles == '.')) {
            return(rep(NA, length(called.alleles)));
        }
        called.alleles.matrix <- matrix(called.alleles, ncol = 1);
    } else if (is.matrix(called.alleles)) {
        # Fast-fail for all missing genotypes
        if (all(is.na(called.alleles)) | all(called.alleles == '.')) {
            return(matrix(NA, nrow = nrow(called.alleles), ncol = ncol(called.alleles), dimnames = dimnames(called.alleles)));
        }
        called.alleles.matrix <- called.alleles;
    } else {
        stop("Unrecognized 'called.alleles' format. Must be a vector or a matrix.");
    }

    # Check that called.alleles.matrix has rows corresponding to risk.alleles
    if (nrow(called.alleles.matrix) != length(risk.alleles)) {
        stop('Number of rows in called.alleles must equal length of risk.alleles.');
    }

    # # check that risk.alleles is the same length as called.alleles
    # if (length(called.alleles) != length(risk.alleles)) {
    #     stop('called.alleles and risk.alleles must be the same length.');
    #     }

    # check for missing risk alleles and warn
    if (any(is.na(risk.alleles))) {
        warning('risk.alleles contains missing values, NA will be returned for corresponding dosage.');
        }

    # check that risk.alleles is a vector of capitalized alphabetic characters
    if (!all(grepl('^[A-Z]+$', risk.alleles[!is.na(risk.alleles)]))) { # NA risk alleles allowed
        stop('unrecognized risk.allele format, must be capitalized letters.');
        }

    # # handle totally missing genotypes
    # # if the entire vector is NA or the entire vector is '.', return NA
    # if (all(is.na(called.alleles)) | all(called.alleles == '.')) {
    #     split.alleles <- data.frame(called.alleles, called.alleles);
    #     } else {
    #         # check that called.alleles is a vector of genotypes in allelic notation or '.' separated by a slash or pipe
    #         # "*" characters represent overlapping deletions from an upstream indel and are accepted VCF format
    #         allowed.pattern <- '^((([A-Z]+|\\.|\\*)[/\\|]([A-Z]+|\\.|\\*))|\\.|[A-Z]+)$' # '|' are special chars in regular expressions
    #         passing.alleles <- grepl(allowed.pattern, called.alleles);
    #         passing.alleles[is.na(called.alleles)] <- TRUE; # NA allowed
    #         if (!all(passing.alleles)) {
    #             stop('unrecognized called.alleles format, must be capitalized letters, "." or "*" separated by a slash or pipe.');
    #             }
    #         # replace hemizygous genotypes with a placeholder for easier splitting
    #         # index for non-NA alleles that are missing allele separators:
    #         no.sep.index <- (!grepl('/|\\|', called.alleles) & !is.na(called.alleles) & called.alleles != '.');
    #         called.alleles[no.sep.index] <- paste0(called.alleles[no.sep.index], '/-');
    #         split.alleles <- data.table::tstrsplit(called.alleles, split = c('/|\\|'), keep = c(1,2)); # '|' are special chars in regular expressions
    #         }
    # names(split.alleles) <- c('called.allele.a', 'called.allele.b');

    # Vectorized validation and handling of called alleles
    # "*" characters represent overlapping deletions from an upstream indel and are accepted VCF format
    allowed.pattern <- '^((([A-Z]+|\\.|\\*)[/\\|]([A-Z]+|\\.|\\*))|\\.|[A-Z]+)$';
    passing.alleles <- grepl(allowed.pattern, called.alleles);
    passing.alleles[is.na(called.alleles)] <- TRUE;
    if (!all(passing.alleles)) {
        stop('unrecognized called.alleles format, must be capitalized letters, "." or "*" separated by a slash or pipe.');
        }

    # Replace hemizygous genotypes with a placeholder for consistent splitting
    no.sep.index <- (!grepl('/|\\|', called.alleles) & !is.na(called.alleles) & called.alleles != '.');
    called.alleles[no.sep.index] <- paste0(called.alleles[no.sep.index], '/-');
    called.alleles.matrix <- matrix(called.alleles, nrow = nrow(called.alleles.matrix), ncol = ncol(called.alleles.matrix));

    # Split the entire matrix of alleles into two matrices, one for each allele
    alleles.split <- data.table::tstrsplit(as.vector(called.alleles.matrix), split = "/|\\|", fixed = FALSE);
    alleles.a <- matrix(alleles.split[[1]], nrow = nrow(called.alleles.matrix), ncol = ncol(called.alleles.matrix));
    alleles.b <- matrix(alleles.split[[2]], nrow = nrow(called.alleles.matrix), ncol = ncol(called.alleles.matrix));

    # Replicate risk.alleles across columns for vectorized comparison
    risk.alleles.matrix <- matrix(
        rep(risk.alleles, times = ncol(called.alleles.matrix)),
        nrow = nrow(called.alleles.matrix)
        );

    # Compute dosage (0, 1, 2)
    # Initialize dosage matrix with zeros
    dosage.matrix <- matrix(0L, nrow = nrow(called.alleles.matrix), ncol = ncol(called.alleles.matrix));
    # Add 1 to dosage for each instance of the risk allele
    dosage.matrix <- dosage.matrix + (alleles.a == risk.alleles.matrix);
    dosage.matrix <- dosage.matrix + (alleles.b == risk.alleles.matrix);

    # Handle special cases
    # 1. Check for missing alleles ('NA' or '.') for both NA assignment and warning
    is.missing.a <- is.na(alleles.a) | (alleles.a == '.');
    is.missing.b <- is.na(alleles.b) | (alleles.b == '.');

    # 2. Case where one allele is marked as missing and the other is not (e.g. `./A`)
    # This should return NA and issue a warning
    is.one.missing <- (is.missing.a & !is.missing.b) | (!is.missing.a & is.missing.b);
    if (any(is.one.missing)) {
        warning('some genotypes contain a missing allele, returning NA for corresponding dosage.');
    }

    # 3. Apply the final NA mask
    na.mask <- is.missing.a | is.missing.b | is.na(called.alleles.matrix) | is.na(risk.alleles.matrix);
    dosage.matrix[na.mask] <- NA;

    # Restore the matrix dimensions and dimnames
    dimnames(dosage.matrix) <- dimnames(called.alleles.matrix);

    # If the original input was a vector, convert the output back to a vector
    if (is.vector.input) {
        return(as.vector(dosage.matrix));
    } else {
        return(dosage.matrix);
    }

    return(dosage);
    }

# The function for calculating a dosage value intended to replace missing genotypes.
calculate.missing.genotype.dosage <- function(dosage.matrix) {
    # calculate the mean dosage for each variant
    mean.dosage <- rowMeans(x = dosage.matrix, na.rm = TRUE);

    # replace NaN (from all NA rows) with NA
    mean.dosage[is.nan(mean.dosage)] <- NA;

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
