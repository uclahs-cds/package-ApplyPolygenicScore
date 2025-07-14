test_that(
    'get.variant.by.sample.matrix correctly transforms long data to SNP by sample matrix', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample2', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.variant.by.sample.matrix <- get.variant.by.sample.matrix(
            long.data = test.long.data,
            row.id.cols = 'variant.id',
            value.col = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(2, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );

test_that(
    'get.variant.by.sample.matrix utility correctly handles duplicates', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        expect_error(
            get.variant.by.sample.matrix(
                long.data = test.long.data,
                variant.id = test.long.data$variant.id,
                value.var = 'dosage'
                )
            );
        }
    );

test_that(
    'get.combined.multiallelic.variant.by.sample.matrix correctly transforms long data to SNP by sample matrix', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample2', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.variant.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data,
            variant.id = test.long.data$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(2, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );

test_that(
    'get.combined.multiallelic.variant.by.sample.matrix utility correctly handles duplicates', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.long.data.one.missing.duplicate <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, NA, 3, 4, NA)
            );
        test.long.data.two.missing.duplicates <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(NA, NA, 3, 4, NA)
            );

        test.variant.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data,
            variant.id = test.long.data$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(3 / 2, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );

        test.variant.by.sample.matrix.one.missing.duplicate <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data.one.missing.duplicate,
            variant.id = test.long.data.one.missing.duplicate$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix.one.missing.duplicate,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );

        test.variant.by.sample.matrix.two.missing.duplicates <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data.two.missing.duplicates,
            variant.id = test.long.data.two.missing.duplicates$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix.two.missing.duplicates,
            as.matrix(data.frame(
                sample1 = c(NA, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );
