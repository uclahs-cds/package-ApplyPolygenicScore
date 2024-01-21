test_that(
    'format.chromosome.notation correctly formats choromome names', {
        # check no prefix to prefix conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', 'X', 'Y'),
                chr.prefix = TRUE,
                numeric.sex.chr = FALSE
                ),
            c('chr1', 'chr2', 'chrX', 'chrY')
            );
        # check no prefix to no prefix conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', 'X', 'Y'),
                chr.prefix = FALSE,
                numeric.sex.chr = FALSE
                ),
            c('1', '2', 'X', 'Y')
            );
        # check prefix to no prefix conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('chr1', 'chr2', 'chrX', 'chrY'),
                chr.prefix = FALSE,
                numeric.sex.chr = FALSE
                ),
            c('1', '2', 'X', 'Y')
            );
        # check prefix to prefix conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('chr1', 'chr2', 'chrX', 'chrY'),
                chr.prefix = TRUE,
                numeric.sex.chr = FALSE
                ),
            c('chr1', 'chr2', 'chrX', 'chrY')
            );

        # check alphabetic sex to numeric sex conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', 'X', 'Y'),
                chr.prefix = FALSE,
                numeric.sex.chr = TRUE
                ),
            c('1', '2', '23', '24')
            );
        # check alphabetic sex to alphabetic sex conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', 'X', 'Y'),
                chr.prefix = FALSE,
                numeric.sex.chr = FALSE
                ),
            c('1', '2', 'X', 'Y')
            );
        # check numeric sex to alphabetic sex conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', '23', '24'),
                chr.prefix = FALSE,
                numeric.sex.chr = FALSE
                ),
            c('1', '2', 'X', 'Y')
            );
        # check numeric sex to numeric sex conversion
        expect_equal(
            format.chromosome.notation(
                chromosome = c('1', '2', '23', '24'),
                chr.prefix = FALSE,
                numeric.sex.chr = TRUE
                ),
            c('1', '2', '23', '24')
            );
        }
    );
