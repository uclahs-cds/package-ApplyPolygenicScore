test_that(
    'get.pgs.percentiles correctly checks inputs', {
        # check that pgs is numeric
        expect_error(get.pgs.percentiles('a'), 'pgs must be a numeric vector');

        # check that n.percentiles is an integer
        expect_error(get.pgs.percentiles(1, n.percentiles = 1.5), 'n.percentiles must be an integer');

        }
    );

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

test_that(
    'run.pgs.regression correctly identifies continuous and binary phenotypes', {
        # load test data
        load('data/phenotype.test.data.Rda');
        pgs <- seq(0, 1, length.out = nrow(phenotype.test.data$phenotype.data));
        phenotype.columns <- c('continuous.phenotype', 'binary.phenotype');
        phenotype.data <- phenotype.test.data$phenotype.data[ , phenotype.columns];
        phenotype.data$categorical.phenotype <- rep(c('a', 'b', 'c', 'd', 'e'), 2); # should not be included in the regression

        # run function
        regression.data <- run.pgs.regression(pgs = pgs, phenotype.data = phenotype.data);

        # check that the output is a data frame
        expect_equal(class(regression.data), 'data.frame');

        # check that the data frame has the correct number of rows
        expect_equal(nrow(regression.data), 2);

        # check that the data frame has the correct number of columns
        expect_equal(ncol(regression.data), 6);

        # check that the data frame has the correct column names
        regression.data.expected.colnames <- c('phenotype', 'model', 'beta', 'se', 'p.value', 'r.squared')
        expect_equal(
            colnames(regression.data),
            regression.data.expected.colnames
            );

        # check that the data frame has the correct values
        expect_equal(
            regression.data$phenotype,
            phenotype.columns
            );

        expect_equal(
            regression.data$model,
            c('linear.regression', 'logistic.regression')
            );

        linear.model.expected.results <- summary(lm(continuous.phenotype ~ pgs, data = phenotype.data));
        logistic.model.expected.results <- summary(glm(binary.phenotype ~ pgs, data = phenotype.data, family = binomial));

        expect_equal(
            regression.data$beta,
            c(linear.model.expected.results$coefficients['pgs', 'Estimate'], logistic.model.expected.results$coefficients['pgs', 'Estimate'])
            );
        expect_equal(
            regression.data$se,
            c(linear.model.expected.results$coefficients['pgs', 'Std. Error'], logistic.model.expected.results$coefficients['pgs', 'Std. Error'])
            );
        expect_equal(
            regression.data$p.value,
            c(linear.model.expected.results$coefficients['pgs', 'Pr(>|t|)'], logistic.model.expected.results$coefficients['pgs', 'Pr(>|z|)'])
            );
        expect_equal(
            regression.data$r.squared,
            c(linear.model.expected.results$r.squared, NA)
            );

        }
    );
