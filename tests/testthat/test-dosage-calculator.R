test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from unphased genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A/T', 'T/T', 'A/A', 'A/T', 'T/T'),
                risk.alleles = c('A', 'A', 'T', 'A', 'A', 'T')
                ),
            c(2, 1, 0, 2, 1, 0)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from phased genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A|A', 'A|T', 'T|T', 'A|A', 'A|T', 'T|T'),
                risk.alleles = c('A', 'A', 'T', 'A', 'A', 'T')
                ),
            c(2, 1, 0, 2, 1, 0)
            );
        }
    );

test_that(
    'convert.alleles.to.pgs.dosage calculates dosage correctly from missing genotypes', {
        expect_equal(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('./.', './A', 'A/.'),
                risk.alleles = c('A', 'A', 'T')
                ),
            c(NA, 1, 1)
            );
        }
    );
