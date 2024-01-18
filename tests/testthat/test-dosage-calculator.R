
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

        # check that correct input is accepted
        expect_silent(
            convert.alleles.to.pgs.dosage(
                called.alleles = c('A/A', 'A|T', 'TA/T', 'A/ATTTT', './.', '.'),
                risk.alleles = c('A', 'T', 'A', 'T', 'A', 'T')
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
                called.alleles = c('./.', './A', 'A/.', '.'),
                risk.alleles = c('A', 'A', 'T', 'T')
                ),
            c(NA, 1, 0, NA)
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
        test.vcf <- ApplyPolygenicScore::import.vcf('data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10_duplicated-sample.vcf.gz')

        expect_silent(
            convert.alleles.to.pgs.dosage(
                called.alleles = test.vcf$dat$gt_GT_alleles,
                risk.alleles = test.vcf$dat$ALT
                )
            );
        }
    );
