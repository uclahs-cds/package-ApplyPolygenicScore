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