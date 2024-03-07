test_that(
    'get.variant.by.sample.matrix correctly transforms long data to SNP by sample matrix', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample2', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.variant.by.sample.matrix <- get.variant.by.sample.matrix(
            long.data = test.long.data,
            variant.id = test.long.data$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(2, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                )
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
