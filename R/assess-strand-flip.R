# default indel handling is to just return as is (no flip).
flip.DNA.allele <- function(allele, return.indels.as.missing = FALSE) {
    if (all(is.na(allele))) {
        return(NA);
        }
    # check that allele is a character vector
    if (!is.character(allele)) {
        stop('allele must be a character vector');
        }

    # Verify acceptable SNP alleles
    snp.allele <- allele[nchar(allele) == 1];
    accepted.alleles <- c('A', 'T', 'C', 'G');
    allele.check <- snp.allele[!is.na(snp.allele)] %in% accepted.alleles;

    if (any(!allele.check, na.rm = TRUE)) {
        stop('Invalid allele: ', snp.allele[!is.na(snp.allele)][!allele.check]);
        }

    flipped.allele <- sapply(
        X = allele,
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
                        `NA` = NA
                        );
                }
            }
        );

    names(flipped.allele) <- NULL;
    return(flipped.allele);
    }

assess.strand.flip <- function(
    vcf.ref.allele,
    vcf.alt.allele,
    pgs.ref.allele,
    pgs.effect.allele
    ) {
 
    # check that all inputs are one dimensional character vectors
    if (!is.character(vcf.ref.allele) | !is.character(vcf.alt.allele) | !is.character(pgs.ref.allele) | !is.character(pgs.effect.allele)) {
        stop('alleles must be character vectors');
        }

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

        all.alleles <- c(current.pgs.ref.allele, current.pgs.effect.allele, current.vcf.ref.allele, current.vcf.alt.allele.split);

        # verify that all alleles are character vectors and contain only A, T, C, G
        if (any(!is.character(all.alleles))) {
            stop('alleles must be character vectors');
            }

        accepted.alleles <- c('A', 'T', 'C', 'G', '*');
        allele.check <- unlist(lapply(X = all.alleles, FUN = function(x) all(strsplit(x, '')[[1]] %in% accepted.alleles)));
        if (any(!allele.check)) {
            stop('Invalid allele: ', all.alleles[!allele.check]);
            }

        # check for invalid indels (both Other and Effect cannot be more than one base long)
        if (nchar(current.pgs.ref.allele) > 1 && nchar(current.pgs.effect.allele) > 1) {
            stop('Invalid PGS alleles: ', current.pgs.ref.allele, ' / ', current.pgs.effect.allele);
            }

        if (nchar(current.vcf.ref.allele) > 1 && any(nchar(current.vcf.alt.allele.split) > 1)) {
            stop('Invalid VCF alleles: ', current.vcf.ref.allele, ' / ', current.vcf.alt.allele.split);
            }

        # check if default ref-ref / alt-effect alleles match
        default.ref.check <- current.pgs.ref.allele == current.vcf.ref.allele;
        default.alt.check <- any(current.vcf.alt.allele.split) %in% current.pgs.effect.allele;

        if (default.ref.check & default.alt.check) {
            flip.designation[i] <- 'default';
            flipped.effect.allele[i] <- pgs.effect.allele[i];
            flipped.other.allele[i] <- pgs.ref.allele[i];
            break;
            }

        # check if PGS effect designation is on the VCF reference allele (effect switch)
        effect.switch.ref.check <- current.pgs.ref.allele == current.vcf.effect.allele;
        effect.switch.alt.check <- any(vcf.alt.allele.split) %in% current.pgs.ref.allele;

        if (effect.switch.ref.check && effect.switch.alt.check) {
            effect.switch.candidate <- TRUE;
            } else {
                effect.switch.candidate <- FALSE;
            }

        # attempt strand flip on default ref-ref / alt-effect alleles
        # identify insertion/deletion alleles
        insertion <- FALSE;
        deletion <- FALSE;
        if (nchar(pgs.ref.allele[i] > 1 && nchar(pgs.effect.allele[i] == 1))) {
            insertion <- TRUE;
            } else if (nchar(pgs.ref.allele[i] == 1 && nchar(pgs.effect.allele[i] > 1))) {
                deletion <- TRUE;
        } else if (nchar(pgs.ref.allele[i] == 1 && nchar(pgs.effect.allele[i] == 1))) {
            insertion <- FALSE;
            deletion <- FALSE;
            } else {
                stop(paste0('Invalid PGS alleles: ', pgs.ref.allele[i], ' / ', pgs.effect.allele[i]));
            }

        pgs.ref.flip <- flip.DNA.allele(pgs.ref.allele[i]);
        pgs.effect.flip <- flip.DNA.allele(pgs.effect.allele[i]);

        default.ref.flip.check <- vcf.ref.allele[i] == pgs.ref.flip;
        default.alt.flip.check <- vcf.alt.allele[i] == pgs.effect.flip;

        if (default.ref.flip.check & default.alt.flip.check) {
            strand.flip.candidate <- TRUE;
            } else {
                strand.flip.candidate <- FALSE;
            }

        if (effect.switch.candidate & strand.flip.candidate) {
            # not possible to determine whether effect switch or strand flip has occured
            # This is an ambiguous case caused by palindromic SNPs
            flip.designation[i] <- 'ambiguous';
            flipped.effect.allele[i] <- NA;
            flipped.other.allele[i] <- NA;
            break;
            } else if (effect.switch.candidate) {
                flip.designation[i] <- 'effect_switch';
                flipped.effect.allele[i] <- pgs.effect.allele[i];
                flipped.other.allele[i] <- pgs.ref.allele[i];
                break;
                } else if (strand.flip.candidate) {
                    flip.designation[i] <- 'strand_flip';
                    flipped.effect.allele[i] <- pgs.effect.flip;
                    flipped.other.allele[i] <- pgs.ref.flip;
                    break;
                    }

        # attempt strand flip on effect_switch ref-effect / alt-other alleles
        effect.switch.ref.flip.check <- vcf.ref.allele[i] == pgs.effect.flip;
        effect.switch.alt.flip.check <- vcf.alt.allele[i] == pgs.ref.flip;

        if(effect.switch.ref.flip.check & effect.switch.alt.flip.check) {
            flip.designation[i] <- 'effect_switch_with_strand_flip';
            flipped.effect.allele[i] <- pgs.effect.flip;
            flipped.other.allele[i] <- pgs.ref.flip;
            break;
            }

        }

    return(
        list(
            strand.status = flip.designation,
            new.pgs.effect.allele = flipped.effect.allele,
            new.pgs.other.allele = flipped.other.allele
            )
        );

    }
