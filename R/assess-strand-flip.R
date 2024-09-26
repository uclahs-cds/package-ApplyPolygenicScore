flip.DNA.allele <- function(allele) {
    if (all(is.na(allele))) {
        return(NA);
        }
    # check that allele is a character vector
    if (!is.character(allele)) {
        stop('allele must be a character vector.');
        }

    accepted.alleles <- c('A', 'T', 'C', 'G');
    allele.check <- allele[!is.na(allele)] %in% accepted.alleles;
    if (any(!allele.check, na.rm = TRUE)) {
        stop('Invalid allele: ', allele[!is.na(allele)][!allele.check]);
        }

    flipped.allele <- sapply(
        X = allele,
        FUN = switch,
        'A' = 'T',
        'T' = 'A',
        'C' = 'G',
        'G' = 'C',
        `NA` = NA
            );
    return(flipped.allele);
    }

assess.strand.flip <- function(
    vcf.ref.allele,
    vcf.alt.allele,
    pgs.ref.allele,
    pgs.effect.allele
    ) {

    # initialize empty vectors to store flip info
    flipped.effect.allele <- rep(NA, length(vcf.ref.allele));
    flipped.other.allele <- rep(NA, length(vcf.ref.allele));
    # strand status
    flip.designation <- rep(NA, length(vcf.ref.allele));

    for (i in 1:length(vcf.ref.allele)) {
        # check if default ref-ref / alt-effect alleles match
        default.ref.check <- vcf.ref.allele[i] == pgs.ref.allele[i];
        default.alt.check <- vcf.alt.allele[i] == pgs.effect.allele[i];

        if (default.ref.check & default.alt.check) {
            flip.designation[i] <- 'default';
            flipped.effect.allele[i] <- pgs.effect.allele[i];
            flipped.other.allele[i] <- pgs.ref.allele[i];
            break;
            }

        # check if PGS effect designation is on the VCF reference allele (effect switch)
        effect.switch.ref.check <- vcf.ref.allele[i] == pgs.effect.allele[i];
        effect.switch.alt.check <- vcf.alt.allele[i] == pgs.ref.allele[i];

        if (effect.switch.ref.check & effect.switch.alt.check) {
            effect.switch.candidate <- TRUE;
            } else {
                effect.switch.candidate <- FALSE;
            }

        # attempt strand flip on default ref-ref / alt-effect alleles
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
