
test_that(
    'convert.alleles.to.pgs.dosage correctly checks input format', {
        # check for correct input lengths
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T', 'T/T', 'A/A', 'A/T', 'T/T'),
                risk.alleles = c('A', 'T', 'A', 'T', 'A', 'T', 'A')
                ),
            'called.alleles and risk.alleles must be the same length.'
            );

        # check for correct risk.alleles format
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T'),
                risk.alleles = c('A', 't')
                ),
            'unrecognized risk.allele format, must be capitalized letters.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T'),
                risk.alleles = c('A', '1')
                ),
            'unrecognized risk.allele format, must be capitalized letters.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T'),
                risk.alleles = c('A', NA)
                ),
            'unrecognized risk.allele format, must be capitalized letters.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T'),
                risk.alleles = c('A', '.')
                ),
            'unrecognized risk.allele format, must be capitalized letters.'
            );

        # check for correct called.alleles format
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A'),
                risk.alleles = c('A', 'T')
                ),
            'unrecognized called.alleles format, must be capitalized letters or "." separated by a slash or pipe.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A,'),
                risk.alleles = c('A', 'T')
                ),
            'unrecognized called.alleles format, must be capitalized letters or "." separated by a slash or pipe.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A-A'),
                risk.alleles = c('A', 'T')
                ),
            'unrecognized called.alleles format, must be capitalized letters or "." separated by a slash or pipe.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'a/A'),
                risk.alleles = c('A', 'T')
                ),
            'unrecognized called.alleles format, must be capitalized letters or "." separated by a slash or pipe.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/'),
                risk.alleles = c('A', 'T', 'A')
                ),
            'called.alleles and risk.alleles must be the same length.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', '/A'),
                risk.alleles = c('A', 'T', 'A')
                ),
            'called.alleles and risk.alleles must be the same length.'
            );
        expect_error(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', NA),
                risk.alleles = c('A', 'T', 'A')
                ),
            'called.alleles and risk.alleles must be the same length.'
            );
        expect_warning(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', './A'),
                risk.alleles = c('A', 'A')
                ),
            'one of two alleles is marked as missing at index 2, this is an unrecognized format, returning NA for dosage.'
            );
        expect_warning(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/.'),
                risk.alleles = c('A', 'A')
                ),
            'one of two alleles is marked as missing at index 2, this is an unrecognized format, returning NA for dosage.'
            );

        # check that correct input is accepted
        expect_silent(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A|T', 'TA/T', 'A/ATTTT', './.', '.', '*/T', 'T/*', '*/*'),
                risk.alleles = c('A', 'T', 'A', 'T', 'A', 'T', 'A', 'T', 'A')
                )
            );
    }
)

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from unphased genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T', 'T/T', 'A/A', 'A/T', 'T/T'),
                risk.alleles = c('A', 'T', 'A', 'T', 'A', 'T')
                ),
            c(2, 1, 0, 0, 1, 2)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from phased genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A|A', 'A|T', 'T|T', 'A|A', 'A|T', 'T|T'),
                risk.alleles = c('A', 'T', 'A', 'T', 'A', 'T')
                ),
            c(2, 1, 0, 0, 1, 2)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from short indel alleles', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/AT', 'AT/A', 'AT/AT', 'A/AT', 'AT/A', 'AT/AT','A/ATCG', 'ATCG/A', 'AT/ATCG'),
                risk.alleles = c('A', 'A', 'A', 'AT', 'AT', 'AT', 'A', 'ATCG', 'A')
                ),
            c(1, 1, 0, 1, 1, 2, 1, 1, 0)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from missing genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('./.', './A', 'A/.', '.', NA),
                risk.alleles = c('A', 'A', 'T', 'T', 'T')
                ),
            c(NA, NA, NA, NA, NA)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly for overlapping deletion alleles', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/*', '*/A', '*/*', 'A/*', '*/A'),
                risk.alleles = c('A', 'A', 'A', 'T', 'T')
                ),
            c(1, 1, 0, 0, 0)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage correctly handles a scenario where no genotypes are provided', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c(NA, NA, NA),
                risk.alleles = c('A', 'A', 'T')
                ),
            c(NA, NA, NA)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage works on real data', {
        test.vcf <- ApplyPolygenicScore::import.vcf('data/HG001_GIAB.vcf.gz')

        expect_silent(
            convert.alleles.to.pgs.dosage(
                called.alleles = test.vcf$dat$gt_GT_alleles,
                risk.alleles = test.vcf$dat$ALT
                )
            );
        }
    );

test_that(
    'calculate.missing.genotype.dosage correctly calculates mean dosage', {
        expect_equal(
            calculate.missing.genotype.dosage(
                dosage.matrix = data.frame(
                    sample1 = c(0, 1, NA),
                    sample2 = c(0, 0, NA),
                    sample3 = c(0, 1, NA),
                    sample4 = c(0, NA, NA),
                    row.names = c('variant1', 'variant2', 'variant3')
                    )
                ),
            c(variant1 = 0, variant2 = 2 / 3, variant3 = NA)
            );
        }
    );

test_that(
    'convert.allele.frequency.to.dosage correctly checks input format', {
        # check for correct input format
        expect_error(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1, 1.1)
                ),
            'allele.frequency must be between 0 and 1.'
            );
        expect_error(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1, -0.1)
                ),
            'allele.frequency must be between 0 and 1.'
            );
        expect_error(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1, 'A')
                ),
            'allele.frequency must be a numeric vector.'
            );
        expect_error(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1, NA)
                ),
            'allele.frequency must not contain NA'
            );

        # check that correct input is accepted
        expect_silent(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1)
                )
            );
        }
    );

test_that(
    'convert.allele.frequency.to.dosage correctly calculates dosage', {
        expect_equal(
            convert.allele.frequency.to.dosage(
                allele.frequency = c(0, 0.5, 1)
                ),
            c(0, 1, 2)
            );
        }
    );
