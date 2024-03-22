test_that(
    'get.pgs.percentiles correctly formats output', {
        pgs <- c(5, 2, 3, 4, 1);
        percentile.data <- get.pgs.percentiles(pgs);

        # check that the output is a data frame
        expect_equal(class(percentile.data), 'data.frame');

        # check that the data frame has the correct number of rows
        expect_equal(nrow(percentile.data), length(pgs));

        # check that the data frame has the correct number of columns
        expect_equal(ncol(percentile.data), 3);

        # check that the data frame has the correct column names
        expect_equal(
            colnames(percentile.data),
            c('percentile', 'decile', 'quartile')
            );
        
        percentile.data.custom <- get.pgs.percentiles(pgs, n.percentiles = 5);

        # check that the data frame has the correct number of columns
        expect_equal(ncol(percentile.data.custom), 4);

        # check that the data frame has the correct column names
        expect_equal(
            colnames(percentile.data.custom),
            c('percentile', 'decile', 'quartile', 'percentile.5')
            );
        }
    );

test_that(
    'get.pgs.percentiles correctly calculates default percentiles', {
        pgs <- c(5, 2, 3, 4, 1);
        percentile.data <- get.pgs.percentiles(pgs);

        # check percentiles
        expect_equal(
            percentile.data$percentile,
            c(1, 0.4, 0.6, 0.8, 0.2)
            );

        # check deciles
        expect_equal(
            percentile.data$decile,
            c(10, 4, 6, 8, 2)
            );

        # check quantiles
        expect_equal(
            percentile.data$quartile,
            c(4, 2, 3, 4, 1)
            );

        }
    );

test_that(
    'get.pgs.percentiles correctly calculates custom percentiles', {
        pgs <- c(5, 2, 3, 4, 1);
        percentile.data <- get.pgs.percentiles(pgs, n.percentiles = 5);

        # check percentiles
        expect_equal(
            percentile.data$percentile.5,
            c(5, 2, 3, 4, 1)
            );
        }
    );