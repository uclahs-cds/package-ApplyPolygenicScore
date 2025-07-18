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

        # edge case when one of the pgss is NA
        pgs.with.na <- c(5, 2, 3, 4, 1, NA);
        percentile.data.with.na <- get.pgs.percentiles(pgs.with.na, n.percentiles = 5);
        expect_equal(
            percentile.data.with.na$percentile.5,
            c(5, 2, 3, 4, 1, NA)
            );
        }
    );

test_that(
    'classify.variable.type correctly identifies continuous and binary variables', {
        # create test data
        continuous.var <- c(1, 2, 3, 4, 5);
        binary.var <- c(1, 2, 1, 2, 1);
        binary.factor <- c('a', 'b', 'a', 'b', 'a');
        categorical.var <- c('a', 'b', 'c', 'd', 'e');
        data <- data.frame(continuous.var, binary.var, binary.factor, categorical.var);

        # run function
        variable.type <- classify.variable.type(data);

        # check that the output is a list
        expect_equal(class(variable.type), 'list');

        # check that the list has the correct number of elements
        expect_equal(length(variable.type), 3);

        # check that the list has the correct element names
        expect_equal(
            names(variable.type),
            c('continuous', 'binary', 'other')
            );

        # check that the list has the correct element values
        expect_equal(
            as.vector(variable.type$continuous),
            c(TRUE, FALSE, FALSE, FALSE)
            );
        expect_equal(
            as.vector(variable.type$binary),
            c(FALSE, TRUE, TRUE, FALSE)
            );
        expect_equal(
            as.vector(variable.type$other),
            c(FALSE, FALSE, FALSE, TRUE)
            );

        }
    );


test_that(
    'run.pgs.regression correctly formats outputs', {
        # load test data
        load('data/phenotype.test.data.Rda');
        pgs <- seq(0, 1, length.out = nrow(phenotype.test.data$phenotype.data));
        phenotype.columns <- c('continuous.phenotype', 'binary.phenotype');
        phenotype.data <- phenotype.test.data$phenotype.data[ , phenotype.columns];
        phenotype.data$binary.factor <- c(rep(c('a', 'b'), 5), 'a'); # should be included in the regression
        phenotype.data$categorical.phenotype <- c(rep(c('a', 'b', 'c', 'd', 'e'), 2), 'f'); # should not be included in the regression
        phenotype.data$binary.factor <- c(rep(c('a', 'b'), 5), 'a');

        # run function
        regression.data <- run.pgs.regression(pgs = pgs, phenotype.data = phenotype.data);

        # check that the output is a data frame
        expect_equal(class(regression.data), 'data.frame');

        # check that the data frame has the correct number of rows
        expect_equal(nrow(regression.data), 3);

        # check that the data frame has the correct number of columns
        expect_equal(ncol(regression.data), 7);

        # check that the data frame has the correct column names
        regression.data.expected.colnames <- c('phenotype', 'model', 'beta', 'se', 'p.value', 'r.squared', 'AUC');
        expect_equal(
            colnames(regression.data),
            regression.data.expected.colnames
            );

        # check that only continuous and binary phenotypes are outputted
        expect_equal(
            regression.data$phenotype,
            c(phenotype.columns, 'binary.factor')
            );

        }
    );

test_that(
    'run.pgs.regression correctly runs regressions', {
        # load test data
        load('data/phenotype.test.data.Rda');
        pgs <- seq(0, 1, length.out = nrow(phenotype.test.data$phenotype.data));
        phenotype.columns <- c('continuous.phenotype', 'binary.phenotype');
        phenotype.data <- phenotype.test.data$phenotype.data[ , phenotype.columns];
        phenotype.data$binary.factor <- c(rep(c('a', 'b'), 5), 'a'); # should be included in the regression
        phenotype.data$categorical.phenotype <- c(rep(c('a', 'b', 'c', 'd', 'e'), 2), 'f'); # should not be included in the regression

        # run function
        regression.data <- run.pgs.regression(pgs = pgs, phenotype.data = phenotype.data);

        # check correct models for continuous and binary phenotypes
        expect_equal(
            regression.data$model,
            c('linear.regression', 'logistic.regression', 'logistic.regression')
            );

        linear.model.expected.results <- summary(lm(continuous.phenotype ~ pgs, data = phenotype.data));
        logistic.model.expected.results <- summary(glm(binary.phenotype ~ pgs, data = phenotype.data, family = binomial));
        logistic.model.factors.expected.results <- summary(glm(factor(binary.factor) ~ pgs, data = phenotype.data, family = binomial));

        expect_equal(
            regression.data$beta,
            c(
                linear.model.expected.results$coefficients['pgs', 'Estimate'],
                logistic.model.expected.results$coefficients['pgs', 'Estimate'],
                logistic.model.factors.expected.results$coefficients['pgs', 'Estimate']
                )
            );
        expect_equal(
            regression.data$se,
            c(
                linear.model.expected.results$coefficients['pgs', 'Std. Error'],
                logistic.model.expected.results$coefficients['pgs', 'Std. Error'],
                logistic.model.factors.expected.results$coefficients['pgs', 'Std. Error']
                )
            );
        expect_equal(
            regression.data$p.value,
            c(
                linear.model.expected.results$coefficients['pgs', 'Pr(>|t|)'],
                logistic.model.expected.results$coefficients['pgs', 'Pr(>|z|)'],
                logistic.model.factors.expected.results$coefficients['pgs', 'Pr(>|z|)']
                )
            );
        expect_equal(
            regression.data$r.squared,
            c(linear.model.expected.results$r.squared, NA, NA)
            );

        }
    );

## Tests for analyze.pgs.binary.predictiveness ##

set.seed(123); # For reproducibility of test data

# Dummy data for testing
test.data <- data.frame( # Changed to dot case
    ID = 1:100,
    PGS.A = rnorm(100), # Changed to dot case
    PGS.B = rnorm(100), # Changed to dot case
    Pheno.Binary.01 = factor(sample(c(0, 1), 100, replace = TRUE), levels = c(0, 1)), # Changed to dot case
    Pheno.Binary.TF = factor(sample(c(FALSE, TRUE), 100, replace = TRUE), levels = c(FALSE, TRUE)), # Changed to dot case
    Pheno.Binary.Numeric.01 = sample(c(0, 1), 100, replace = TRUE), # Changed to dot case
    Pheno.Binary.Char.YesNo = sample(c("Yes", "No"), 100, replace = TRUE), # Changed to dot case
    Pheno.Binary.3Levels = factor(sample(c(0, 1, 2), 100, replace = TRUE), levels = c(0,1,2)), # Changed to dot case
    Pheno.Continuous.Num = rnorm(100, mean = 50, sd = 10), # Changed to dot case
    Cov.Age = rnorm(100, mean = 40, sd = 10), # Changed to dot case
    Cov.Sex = factor(sample(c("M", "F"), 100, replace = TRUE)), # Changed to dot case
    NonNumericPGS = as.character(rnorm(100)), # Changed to dot case
    stringsAsFactors = FALSE
    );

# Introduce NAs randomly into some columns for robustness testing
# Introduce NAs into PGS columns
test.data$PGS.A[sample(1:nrow(test.data), 5)] <- NA;
test.data$PGS.B[sample(1:nrow(test.data), 3)] <- NA;

# Introduce NAs into phenotype columns
test.data$Pheno.Binary.01[sample(1:nrow(test.data), 4)] <- NA;
test.data$Pheno.Binary.Numeric.01[sample(1:nrow(test.data), 2)] <- NA;
test.data$Pheno.Continuous.Num[sample(1:nrow(test.data), 6)] <- NA;

# Introduce NAs into covariate columns
test.data$Cov.Age[sample(1:nrow(test.data), 3)] <- NA;

# Create a temporary directory for plotting output for tests that might try to write files
temp.output.dir <- tempdir();

test_that(
    'analyze.pgs.binary.predictiveness correctly validates core inputs', {
        # Test 1: `data` must be a data frame
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = list(),
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary'
            ),
            '`data` must be a data frame.'
        );

        # Test 2: Not all specified `pgs.columns` found in the data frame
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = c('PGS.A', 'NonExistentPGS'),
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary'
            ),
            'Not all specified `pgs.columns` found in the data frame.'
        );

        # Test 3: All `pgs.columns` must be numeric
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'NonNumericPGS',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary'
            ),
            'All `pgs.columns` must be numeric.'
        );

        # Test 4: Not all specified `phenotype.columns` found in the data frame
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = c('Pheno.Binary.01', 'NonExistentPheno'),
                phenotype.type = 'binary'
            ),
            'Not all specified `phenotype.columns` found in the data frame.'
        );

        # Test 5: `covariate.columns` not found
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                covariate.columns = 'NonExistentCovariate',
                phenotype.type = 'binary'
            ),
            'Not all specified `covariate.columns` found in the data frame.'
        );

        # Test 6: `phenotype.type` is invalid (using match.arg's error)
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'invalid_type'
            ),
            '`phenotype.type` must be either "binary" or "continuous"'
        );
    }
);

test_that(
    'analyze.pgs.binary.predictiveness validates binary phenotype columns', {

        # Test 7a: Binary phenotype column is a factor but not 2 levels
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.3Levels',
                phenotype.type = 'binary'
            ),
            'must have exactly two levels.'
        );

        # Test 7b: Binary phenotype column not a factor and not 0/1 numeric
        # Pheno.Binary.Char.YesNo is character, not 0/1 numeric
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.Char.YesNo',
                phenotype.type = 'binary'
            ),
            'is specified as binary but is not a factor or 0/1 numeric. Convert to factor'
        );

        # Test 7c: Numeric binary phenotype is automatically converted to factor.
        expect_warning( # Expect no error or warning for valid conversion
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.Numeric.01',
                phenotype.type = 'binary'
            ),
            'is not a factor. Attempting to convert to factor.'
        );

        # Test 7d: Binary phenotype that converts successfully (no error expected, but good for coverage)
        expect_silent( # Expect no error or warning for valid conversion
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary'
            )
        );
    }
);

test_that(
    'analyze.pgs.binary.predictiveness validates continuous phenotype columns and cutoff.threshold', {

        # Test 8: Continuous phenotype column not numeric
        test.data.bad.continuous <- test.data;
        test.data.bad.continuous$Pheno.Continuous.Num <- as.character(test.data.bad.continuous$Pheno.Continuous.Num);
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data.bad.continuous,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Continuous.Num',
                phenotype.type = 'continuous',
                cutoff.threshold = 50
            ),
            'is specified as continuous but is not numeric.'
        );

        # Test 9: `cutoff.threshold` missing for continuous phenotype
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Continuous.Num',
                phenotype.type = 'continuous',
                cutoff.threshold = NULL # Explicitly set to NULL
            ),
            '`cutoff.threshold` must be provided for `continuous` phenotype.type.'
        );

        # Test 10: `cutoff.threshold` is a list but missing entry for specific phenotype
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Continuous.Num',
                phenotype.type = 'continuous',
                cutoff.threshold = list(AnotherPheno = 60) # Missing entry for Pheno.Continuous.Num
            ),
            'must contain an entry for \'Pheno.Continuous.Num\' when it is a named list.'
        );

        # Test 11: `cutoff.threshold` is not numeric or a named list
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Continuous.Num',
                phenotype.type = 'continuous',
                cutoff.threshold = 'not_numeric' # Invalid type
            ),
            '`cutoff.threshold` must be numeric or a named list for `continuous` phenotype.type.'
        );

        # Test 12: `cutoff.threshold` is a list but value for specific phenotype is NULL/not numeric
        expect_error(
            analyze.pgs.binary.predictiveness(
                data = test.data,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Continuous.Num',
                phenotype.type = 'continuous',
                cutoff.threshold = list(Pheno.Continuous.Num = NULL) # Value is NULL inside list
            ),
            'Cutoff threshold for phenotype \'Pheno.Continuous.Num\' is missing or not numeric.'
        );
    }
);

test_that(
    'analyze.pgs.binary.predictiveness handles edge cases and warnings', {

        # Test 13: Too many PGS columns for default.colours (if 'scales' is NOT installed)
        # This test is conditional on 'scales' package availability.
        # If 'scales' is installed, this test will be skipped.
        if (isFALSE(requireNamespace("scales", quietly = TRUE))) {
            many.pgs.data <- data.frame(
                ID = 1:10,
                Pheno.Binary.01 = factor(sample(c(0,1), 10, replace = TRUE), levels = c(0,1))
            );
            for (i in 1:13) { # 13 PGS columns > max.colors (12)
                many.pgs.data[[paste0('PGS', i)]] <- rnorm(10);
            }
            expect_error(
                analyze.pgs.binary.predictiveness(
                    data = many.pgs.data,
                    pgs.columns = paste0('PGS', 1:13),
                    phenotype.columns = 'Pheno.Binary.01',
                    phenotype.type = 'binary',
                    output.dir = temp.output.dir,
                    filename.prefix = "test"
                ),
                'Too many PGS columns selected, default.colours supports 12 or fewer.'
            );
        } else {
            message("Skipping test for >12 PGS columns as 'scales' package is available and will be used.");
        }


        # Test 14: No complete cases after filtering (expects warning and NA results)
        test.data.no.complete.cases <- data.frame(
            ID = 1:10,
            PGS.A = c(rnorm(2), rep(NA, 8)), # Mostly NA
            Pheno.Binary.01 = factor(c(rep(NA, 8), sample(c(0,1), 2, replace=TRUE)), levels=c(0,1))
        );
        expect_warning(
            results <- analyze.pgs.binary.predictiveness(
                data = test.data.no.complete.cases,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary',
                output.dir = temp.output.dir
            ),
            'No complete cases for PGS \'PGS.A\' and phenotype \'Pheno.Binary.01\' after NA removal. Skipping.'
        );
        expect_true(is.na(results$results.df$AUC));


        # Test 15: Logistic regression fails (e.g., due to no variance in PGS)
        test.data.bad.glm <- data.frame(
            ID = 1:10,
            PGS.A = rep(1, 10), # No variance in PGS
            Pheno.Binary.01 = factor(sample(c(0,1), 10, replace=TRUE), levels=c(0,1))
        );
        expect_warning(
            results <- analyze.pgs.binary.predictiveness(
                data = test.data.bad.glm,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary',
                output.dir = temp.output.dir
            ),
            'PGS term \'PGS.A\' in model for phenotype \'Pheno.Binary.01\' is NA or Inf'
        );
        expect_true(is.na(results$results.df$AUC));
        expect_true(is.na(results$results.df$OR));
        expect_true(is.na(results$results.df$p.value));

        # Test 16: ROC calculation fails (e.g., phenotype is constant after complete.cases)
        test.data.bad.roc <- data.frame(
            ID = 1:10,
            PGS.A = rnorm(10),
            Pheno.Binary.01 = factor(rep(0, 10), levels=c(0,1)) # All phenotype values are 0
            );
        expect_warning(
            results <- analyze.pgs.binary.predictiveness(
                data = test.data.bad.roc,
                pgs.columns = 'PGS.A',
                phenotype.columns = 'Pheno.Binary.01',
                phenotype.type = 'binary',
                output.dir = temp.output.dir
            ),
            'ROC curve calculation failed for PGS \'PGS.A\' and phenotype \'Pheno.Binary.01\''
        );
        expect_true(is.na(results$results.df$AUC));


        # Test 17: PGS term not found in model coefficients
        # Catastrophic error, don't know how to trigger.

    }
);

## Tests for analyze.pgs.binary.predictiveness outputs ##

test_that(
    'analyze.pgs.binary.predictiveness correctly formats outputs', {

        # Scenario 1: Binary phenotypes, with covariates, return plot object
        results.binary <- analyze.pgs.binary.predictiveness(
            data = test.data, # Using the global test.data
            pgs.columns = c('PGS.A', 'PGS.B'),
            phenotype.columns = c('Pheno.Binary.01', 'Pheno.Binary.TF'),
            covariate.columns = c('Cov.Age', 'Cov.Sex'),
            phenotype.type = 'binary',
            output.dir = NULL # Return plot object
            );

        # Check overall return type
        expect_type(results.binary, 'list');
        expect_named(results.binary, c('results.df', 'roc.plot'));

        # Check results.df
        expect_s3_class(results.binary$results.df, 'data.frame');
        # Expected rows: (number of PGS) * (number of phenotypes) = 2 * 2 = 4
        expect_equal(nrow(results.binary$results.df), 4);
        # Expected columns: Phenotype, PGS, AUC, OR, OR.Lower.CI, OR.Upper.CI, P.Value
        expect_equal(ncol(results.binary$results.df), 7);

        expected.colnames <- c('phenotype', 'PGS', 'AUC', 'OR', 'OR.Lower.CI', 'OR.Upper.CI', 'p.value');
        expect_equal(colnames(results.binary$results.df), expected.colnames);

        # Check column types (accounting for possible NAs in input leading to NA results)
        expect_type(results.binary$results.df$phenotype, 'character');
        expect_type(results.binary$results.df$PGS, 'character');
        # AUC, OR, P.Value should be numeric even if values are NA
        expect_type(results.binary$results.df$AUC, 'double');
        expect_type(results.binary$results.df$OR, 'double');
        expect_type(results.binary$results.df$OR.Lower.CI, 'double');
        expect_type(results.binary$results.df$OR.Upper.CI, 'double');
        expect_type(results.binary$results.df$p.value, 'double');

        # Check plausible values for AUC (between 0 and 1, or NA if calculation failed)
        # We expect some NA values due to NAs in test.data, so cannot check all non-NA values
        expect_true(all(is.na(results.binary$results.df$AUC) | (results.binary$results.df$AUC >= 0 & results.binary$results.df$AUC <= 1)));

        # Check roc.plot object (should be a BoutrosLab.plotting.general multipanel object)
        # Only check if plotting is enabled or if the plotting package is available.
        if (requireNamespace('BoutrosLab.plotting.general', quietly = TRUE)) {
            expect_s3_class(results.binary$roc.plot, 'multipanel'); # BoutrosLab specific class
            }


        # Scenario 2: Continuous phenotype, no covariates, save to file (roc.plot should be NULL in return)
        # Use existing test.data
        results.continuous <- analyze.pgs.binary.predictiveness(
            data = test.data, # Using the global test.data
            pgs.columns = 'PGS.A',
            phenotype.columns = 'Pheno.Continuous.Num',
            covariate.columns = NULL,
            phenotype.type = 'continuous',
            cutoff.threshold = 50, # Use a common threshold
            output.dir = temp.output.dir, # Save to file
            filename.prefix = 'TEST-cont-ROC'
            );

        # Check overall return type and elements
        expect_type(results.continuous, 'list');
        expect_named(results.continuous, c('results.df', 'roc.plot'));

        # Check results.df format (should be consistent)
        expect_s3_class(results.continuous$results.df, 'data.frame');
        expect_equal(nrow(results.continuous$results.df), 1 * 1); # 1 PGS * 1 Pheno
        expect_equal(ncol(results.continuous$results.df), 7);
        expect_equal(colnames(results.continuous$results.df), expected.colnames);
        expect_true(all(is.na(results.continuous$results.df$AUC) | (results.continuous$results.df$AUC >= 0 & results.continuous$results.df$AUC <= 1)));

        # Check roc.plot object (should be NULL because output.dir was provided)
        expect_null(results.continuous$roc.plot);

        # Check that the output file was created
        output.file <- generate.filename(
            project.stem = 'TEST-cont-ROC',
            file.core = 'pgs-roc-curves',
            extension = 'png'
            );

        expect_true(file.exists(file.path(temp.output.dir, output.file)));

        # Clean up the generated plot file
        file.remove(list.files(temp.output.dir, pattern = '*\\.png$', full.names = TRUE));
    }
);
