#' @title Convert alleles to dosage
#' @description Convert genotype calls in the form of witten out alleles (e.g. 'A/T') to dosages based on provided risk alleles from a PGS (0, 1, or 2).
#' @param called.alleles A vector of genotypes in allelic notation separated by a slash or pipe.
#' @param risk.alleles A vector of risk alleles from a polygenic score corresponding to each genotype (by locus) in called.alleles.
#' @return A vector of dosages corresponding to each genotype in called.alleles.
#' @export
convert.alleles.to.pgs.dosage <- function(called.alleles, risk.alleles) {
    # check that risk.alleles is the same length as called.alleles
    if (length(called.alleles) != length(risk.alleles)) {
        stop('called.alleles and risk.alleles must be the same length.');
        }
    
    # check that risk.alleles is a vector of capitalized alphabetic characters
    if (!all(grepl('^[A-Z]+$', risk.alleles))) {
        stop('unrecognized risk.allele format, must be capitalized letters.');
        }

    # handle totally missing genotypes
    # if the entire vector is NA or the entire vector is '.', return NA
    if (all(is.na(called.alleles)) | all(called.alleles == '.')) {
        split.alleles <- data.frame(called.alleles, called.alleles);
        } else {
            # check that called.alleles is a vector of genotypes in allelic notation or '.' separated by a slash or pipe
            allowed.pattern <- '^((([A-Z]+|\\.)[/\\|]([A-Z]+|\\.))|\\.)$' # '|' are special chars in regular expressions
            if (!all(grepl(allowed.pattern, called.alleles))) {
                stop('unrecognized called.alleles format, must be capitalized letters or "." separated by a slash or pipe.');
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
        if ((split.alleles$called.allele.a[i] == missing.label) & (split.alleles$called.allele.b[i] == missing.label)) {
            dosage[i] <- NA; # if no genotype was called, return NA
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
