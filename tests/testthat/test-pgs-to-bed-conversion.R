test_that(
    'convert.pgs.to.bed correctly checks input', {
        # check that data is a data.frame
        expect_error(
            convert.pgs.to.bed(pgs.weight.data = 'not.a.data.frame'),
            'data must be a data.frame'
            );

        # check that data has CHROM and POS columns
        expect_error(
            convert.pgs.to.bed(pgs.weight.data = data.frame(
                chr = c('1', '2', '3'),
                POS = c(1, 2, 3)
                )),
            'data must have CHROM and POS columns'
            );

        # check that slop is a non-negative integer
        expect_error(
            convert.pgs.to.bed(pgs.weight.data = data.frame(
                CHROM = c('1', '2', '3'),
                POS = c(1, 2, 3)
                ),
                slop = 'not.a.non-negative.integer'
                ),
            'slop must be a non-negative integer'
            );

        # check that slop is a non-negative integer
        expect_error(
            convert.pgs.to.bed(pgs.weight.data = data.frame(
                CHROM = c('1', '2', '3'),
                POS = c(1, 2, 3)
                ),
                slop = -1
                ),
            'slop must be a non-negative integer'
            );

        # check that negative start coordinates are caught
        expect_error(
            convert.pgs.to.bed(pgs.weight.data = data.frame(
                CHROM = c('1', '2', '3'),
                POS = c(0, 2, 3)
                ),
                slop = 0
                ),
            '0-indexing caused negative start coordinates.'
            );

        # check that correct input is accepted
        expect_silent(
            convert.pgs.to.bed(pgs.weight.data = data.frame(
                CHROM = c('1', '2', '3'),
                POS = c(1, 2, 3)
                ))
            );
        }
    );

test_that(
    'convert.pgs.to.bed returns a BED formatted data frame', {
        load('data/tiny.pgs.test.data.Rda');

        # convert tiny.pgs.test.data to BED format
        tiny.pgs.test.data.bed <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data);

        # check that the output is a data.frame
        expect_equal(
            class(tiny.pgs.test.data.bed),
            'data.frame'
            );

        # check that the output has the correct number of rows and columns
        expect_equal(
            dim(tiny.pgs.test.data.bed),
            c(8, 6)
            );

        # check that the output has the correct column names
        expect_equal(
            colnames(tiny.pgs.test.data.bed),
            c('chr', 'start', 'end', 'effect_allele', 'effect_weight', 'beta')
            );

        }
    );


test_that(
    'convert.pgs.to.bed correctly formats output chr column', {
        load('data/tiny.pgs.test.data.Rda');
        test.no.chr <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, chr.prefix = FALSE);
        test.yes.chr <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, chr.prefix = TRUE);
        test.numeric.sex.chr <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, numeric.sex.chr = TRUE);
        test.no.numeric.sex.chr <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, numeric.sex.chr = FALSE);

        # check that CHROM prefix is correctly added
        expect_equal(
            test.no.chr$chr,
            c('1', '2', '3', '4', 'X', 'Y', 'X', 'Y')
            );

        # check that CHROM prefix is correctly removed
        expect_equal(
            test.yes.chr$chr,
            c('chr1', 'chr2', 'chr3', 'chr4', 'chrX', 'chrY', 'chrX', 'chrY')
            );

        # check that CHROM sex chromosomes are correctly converted to numeric
        expect_equal(
            test.numeric.sex.chr$chr,
            c('chr1', 'chr2', 'chr3', 'chr4', 'chr23', 'chr24', 'chr23', 'chr24')
            );

        # check that CHROM sex chromosomes are correctly converted to alphabetic
        expect_equal(
            test.no.numeric.sex.chr$chr,
            c('chr1', 'chr2', 'chr3', 'chr4', 'chrX', 'chrY', 'chrX', 'chrY')
            );

        }

    );

test_that(
    'convert.pgs.to.bed correctly formats output coordinates', {
        load('data/tiny.pgs.test.data.Rda');
        test.no.slop <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 0);

        # check that no negative start coordinate warnings are issued
        expect_silent(
            convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 0)
            );

        # check that coordinates are correctly converted to 0-index
        expect_equal(
            test.no.slop$start,
            c(0, 9, 99, 999, 9999, 99999, 999999, 9999999)
            );
        expect_equal(
            test.no.slop$end,
            c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
            );

        }
    );

test_that(
    'convert.pgs.to.bed correctly adds slop', {
        load('data/tiny.pgs.test.data.Rda');
        test.small.slop <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 1);

        # check that warning is issued for negative slop coordinates
        expect_warning(
            convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 1)
            );
        # check that slop is correctly added
        expect_equal(
            test.small.slop$start,
            c(0, 8, 98, 998, 9998, 99998, 999998, 9999998)
            );
        expect_equal(
            test.small.slop$end,
            c(2, 11, 101, 1001, 10001, 100001, 1000001, 10000001)
            );

        test.large.slop <- convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 100);

        # check that a warning is issued for negative start coordinates
        expect_warning(
            convert.pgs.to.bed(pgs.weight.data = tiny.pgs.test.data, slop = 100),
            'Slop caused negative start coordinates; replacing with 0.'
            );

        # check that slop is correctly added
        expect_equal(
            test.large.slop$start,
            c(0, 0, 0, 899, 9899, 99899, 999899, 9999899)
            );
        expect_equal(
            test.large.slop$end,
            c(101, 110, 200, 1100, 10100, 100100, 1000100, 10000100)
            );

        }
    );

test_that(
    'convert.pgs.to.bed works on real data', {
        pgs <- import.pgs.weight.file('data/PGS000662_hmPOS_GRCh38.txt')
        expect_silent(
            convert.pgs.to.bed(pgs.weight.data = pgs$pgs.weight.data)
            );
        }
    );
