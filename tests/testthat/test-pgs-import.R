TEST.FILE.COLNAMES <- c(
    'rsID',
    'chr_name',
    'chr_position',
    'effect_allele',
    'other_allele',
    'effect_weight',
    'allelefrequency_effect_European',
    'allelefrequency_effect_African',
    'allelefrequency_effect_Asian',
    'allelefrequency_effect_Hispanic',
    'hm_source',
    'hm_rsID',
    'hm_chr',
    'hm_pos',
    'hm_inferOtherAllele'
    );

TEST.FILE.FIRST.ROW <- c(
    'rs7542260',
    '1',
    '5743196',
    'T',
    'C',
    '0.102298257',
    '0.067',
    '0.439',
    '0.113',
    '0.157',
    'ENSEMBL',
    'rs7542260',
    '1',
    '5683136',
    NA
    );

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
        expect_true(
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

test_that(
    'import.pgs.weight.file correctly outputs a list of two data frames', {
        # import beta weights
        beta.weights <- import.pgs.weight.file(input = 'data/PGS000662_hmPOS_GRCh38.txt', use.harmonized.data = TRUE);

        # check that the output is a list of two data frames
        expect_equal(
            length(beta.weights),
            2
            );
        expect_equal(
            class(beta.weights),
            'list'
            );
        expect_equal(
            class(beta.weights[[1]]),
            'data.frame'
            );
        expect_equal(
            class(beta.weights[[2]]),
            'data.frame'
            );
        expect_equal(
            names(beta.weights),
            c('file.metadata', 'pgs.weight.data')
            );

        # check that pgs.weight.data is correctly formatted
        expect_equal(
            dim(beta.weights$pgs.weight.data),
            c(269, 18)
            );

        expected.colnames <- c(TEST.FILE.COLNAMES, 'CHROM', 'POS', 'beta');

        expect_equal(
            colnames(beta.weights$pgs.weight.data),
            expected.colnames
            );

        }
    );

test_that(
    'import.pgs.weight.file correctly formats columns for harmonized data', {
        # import beta weights
        beta.weights <- import.pgs.weight.file(input = 'data/PGS000662_hmPOS_GRCh38.txt', use.harmonized.data = TRUE);

        # check that the harmonized data columns are correctly formatted
        expect_equal(
            beta.weights$pgs.weight.data$CHROM,
            beta.weights$pgs.weight.data$hm_chr
            );
        expect_equal(
            beta.weights$pgs.weight.data$POS,
            beta.weights$pgs.weight.data$hm_pos
            );

        # check first row of data frame

        expected.first.row <- c(TEST.FILE.FIRST.ROW, '1', '5683136', '0.102298257');

        expect_equal(
            as.vector(unlist(beta.weights$pgs.weight.data[1,])),
            expected.first.row
            );
        }
    );

test_that(
    'import.pgs.weight.file correctly formats columns for non-harmonized data', {
        # import beta weights
        beta.weights <- import.pgs.weight.file(input = 'data/PGS000662_hmPOS_GRCh38.txt', use.harmonized.data = FALSE);

        # check that the harmonized data columns are correctly formatted
        expect_equal(
            beta.weights$pgs.weight.data$CHROM,
            beta.weights$pgs.weight.data$chr_name
            );
        expect_equal(
            beta.weights$pgs.weight.data$POS,
            beta.weights$pgs.weight.data$chr_position
            );

        # check first row of data frame
        expected.first.row <- c(TEST.FILE.FIRST.ROW, '1', '5743196', '0.102298257');

        expect_equal(
            as.vector(unlist(beta.weights$pgs.weight.data[1,])),
            expected.first.row
            );
        }
    );


test_that(
    'import.pgs.weight.file correctly formats columns for non-harmonized data', {
        # import beta weights
        beta.weights <- import.pgs.weight.file(input = 'data/PGS000662_hmPOS_GRCh38.txt', use.harmonized.data = FALSE);

        # check that the harmonized data columns are correctly formatted
        expect_equal(
            beta.weights$pgs.weight.data$CHROM,
            beta.weights$pgs.weight.data$chr_name
            );
        expect_equal(
            beta.weights$pgs.weight.data$POS,
            beta.weights$pgs.weight.data$chr_position
            );
        }
    );

test_that(
    'import.pgs.weight.file recognizes unreported weight format', {
      # check that a warning is issued
        expect_warning(
            import.pgs.weight.file(input = 'data/PGS003378_hmPOS_GRCh38_weight-NR.txt', use.harmonized.data = FALSE),
            'Weight format was not reported in the PGS file header. Assuming beta weights.'
            );

        # check that the weight format is assumed to be beta
        beta.weights <- import.pgs.weight.file(input = 'data/PGS003378_hmPOS_GRCh38_weight-NR.txt', use.harmonized.data = FALSE);
        expect_equal(
            beta.weights$pgs.weight.data$beta,
            as.numeric(beta.weights$pgs.weight.data$effect_weight)
            );

        }
    );

test_that(
    'import.pgs.weight.file correctly formats OR/HR weights', {
        # import beta weights
        beta.weights <- import.pgs.weight.file(input = 'data/PGS003378_hmPOS_GRCh38_weight-OR.txt', use.harmonized.data = FALSE);

        # check that warning is issued
        expect_warning(
            import.pgs.weight.file(input = 'data/PGS003378_hmPOS_GRCh38_weight-OR.txt', use.harmonized.data = FALSE),
            'OR/HR weights were converted to beta weights.'
            );

        # check that OR weights are converted to beta
        expect_equal(
            beta.weights$pgs.weight.data$beta,
            log(as.numeric(beta.weights$pgs.weight.data$effect_weight))
            );
        }
    );
