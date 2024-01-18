#' @title Format chromosome names
#' @description Format chromosome names according to user specifications.
#' @param chromosome A character vector of chromosome names.
#' @param chr.prefix A logical indicating whether the 'chr' prefix should be used when formatting chromosome name.
#' @param numeric.sex.chr A logical indicating whether the sex chromosomes should be formatted numerically, as opposed to alphabetically.
#' @return A character vector of chromosome names formatted according to user specifications.
#' @export
format.chromosome.notation <- function(chromosome, chr.prefix, numeric.sex.chr) {
    chromosome <- as.charater(chromosome);

    # convert chromosome to default format (no 'chr' prefix, alphabetic sex chromosomes)
    chromosome <- gsub('^chr', '', chromosome);
    chromosome <- gsub('23$', 'X', chromosome);
    chromosome <- gsub('24$', 'Y', chromosome);

    # apply requested chromosome formatting
    if (chr.prefix) {
        chromosome <- paste0('chr', chromosome);
        }
    if (numeric.sex.chr) {
        chromosome <- gsub('X$', '23', chromosome);
        chromosome <- gsub('Y$', '24', chromosome);
        }

    return(chromosome);
    }
