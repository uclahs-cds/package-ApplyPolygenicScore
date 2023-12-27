test_that(
    'check.pgs.weight.columns detects missing generic columns', {
        # check that required columns are present
        expect_error(
            check.pgs.weight.columns(x = c('chr_name', 'chr_position', 'effect_allele', 'foo', 'bar'), harmonized = FALSE),
            'The following required columns are missing from the PGS weight file: effect_weight'
            );
        }
    );

test_that(
    'check.pgs.weight.columns detects missing harmonized columns', {
        # check that required columns are present
        expect_error(
            check.pgs.weight.columns(x = c('chr_name', 'chr_position', 'effect_allele', 'effect_weight', 'hm_chr'), harmonized = TRUE),
            'The following required columns are missing from the PGS weight file: hm_pos'
            );
        }
    );

test_that(
    'check.pgs.weight.columns accepts complete column set', {
        # check that required columns are present
        expect_silent(
            check.pgs.weight.columns(x = c('chr_name', 'chr_position', 'effect_allele', 'effect_weight', 'hm_chr', 'hm_pos'), harmonized = TRUE)
            );
        }
    );

test_that(
    'parse.pgs.input.header works correctly on unzipped input', {
        load('data/import.test.data.Rda');
        expect_equal(
            parse.pgs.input.header(input = 'data/PGS000662_hmPOS_GRCh38.txt'),
            import.test.data$PGS00662.metadata
            );
        }
    );

test_that(
    'parse.pgs.input.header works correctly on zipped input', {
        load('data/import.test.data.Rda');
        expect_equal(
            parse.pgs.input.header(input = 'data/PGS000662_hmPOS_GRCh38.txt.gz'),
            import.test.data$PGS00662.metadata
            );
        }
    );
