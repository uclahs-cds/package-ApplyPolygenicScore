# plotting functions take a long time to run, this var toggles off plotting tests for faster testing
SKIP.PLOTS <- FALSE;
SKIP.COMPREHENSIVE.CASES <- FALSE;
skip.plotting.tests <- function(skip.plots = FALSE) {
    if (skip.plots) {
        skip('Plotting tests disabled');
        }
    }

set.seed(123);

# get test dataset
load('data/phenotype.test.data.Rda');

# run pgs statistics
pgs.test <- apply.polygenic.score(
    vcf.data = phenotype.test.data$vcf.data,
    pgs.weight.data = phenotype.test.data$pgs.weight.data,
    phenotype.data = phenotype.test.data$phenotype.data,
    missing.genotype.method = c('mean.dosage', 'normalize', 'none'),
    n.percentiles = 2
    )$pgs.output;
# add missing genotpye counts
pgs.test$n.missing.genotypes <- sample(1:10, nrow(pgs.test), replace = TRUE);
pgs.test$percent.missing.genotypes <- round(sample(1:100, nrow(pgs.test), replace = TRUE) / 100, 2);
# add some missing data
pgs.test$continuous.phenotype[1:2] <- NA;
pgs.test$binary.phenotype[2:3] <- NA;
pgs.test$categorical.phenotype[3:4] <- NA;
pgs.test$decile[4:5] <- NA;
pgs.test$percentile[4:5] <- NA;

test_that(
    'create.pgs.density.plot correctly validates inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check that input data is a data frame
        expect_error(
            create.pgs.density.plot(
                pgs.data = list()
                ),
            'pgs.data must be a data.frame'
            );

        # check that the required columns are present
        expect_error(
            create.pgs.density.plot(
                pgs.data = data.frame(not.recognized.PGS.colum = 1:10)
                ),
            'No recognized PGS columns found in pgs.data'
            );

        # check that phenotype.columns is a character vector
        expect_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                phenotype.columns = c(1,2,3)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.columns are present in pgs.data
        expect_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('missing.phenotype')
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that phenotype.columns do not contain recognized PGS columns
        expect_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('PGS.with.replaced.missing')
                ),
            'phenotype.columns cannot contain recognized PGS column names'
            );
        # check that output.dir is a real directory
        expect_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                output.dir = 'not/a/real/directory'
                ),
            'not/a/real/directory does not exist'
            );
        }
    );

test_that(
    'create.pgs.density.plot runs with no error with basic inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs density
        expect_no_error(
            create.pgs.density.plot(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing, PGS.with.normalized.missing)),
                phenotype.columns = NULL,
                output.dir = temp.dir,
                filename.prefix = 'TEST'
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'pgs-density',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        # check returned object
        test.plot.object <- create.pgs.density.plot(
            pgs.data = pgs.test,
            filename.prefix = 'TEST',
            output.dir = NULL
            );
        expect_equal(
            class(test.plot.object),
            'multipanel'
            );

        }
    );

test_that(
    'create.pgs.densityplot runs correctly with large number of phenotype categories', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # extend test data data frame by repeating itself a bunch of times
        very.categorical.pgs.test <- pgs.test[rep(1:nrow(pgs.test), 73 * 4), ];
        # simulate a PGS column
        very.categorical.pgs.test$PGS <- rnorm(nrow(very.categorical.pgs.test));
        # add more categories to categorical phenotype
        very.categorical.pgs.test$very.categorical.phenotype <- as.character(rep(1:73, each = nrow(pgs.test)));

        expect_warning(
            create.pgs.density.plot(
                pgs.data = very.categorical.pgs.test,
                phenotype.columns = 'very.categorical.phenotype',
                output.dir = temp.dir,
                filename.prefix = 'TEST-many-categories'
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-many-categories',
            file.core = 'pgs-density',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.density.plot runs correctly with tidy titles enabled', {
        skip.plotting.tests(skip.plots = SKIP.COMPREHENSIVE.CASES || SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs density
        expect_no_error(
            create.pgs.density.plot(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing, PGS.with.normalized.missing)),
                phenotype.columns = NULL,
                output.dir = temp.dir,
                filename.prefix = 'TEST-tidy-titles',
                tidy.titles = TRUE
                )
            );

        }
    );


test_that(
    'create.pgs.density.plot runs correctily with user provided phenotypes', {
        skip.plotting.tests(skip.plots = SKIP.COMPREHENSIVE.CASES || SKIP.PLOTS);

        temp.dir <- tempdir();

        # check handling of many phenotypes
        expect_no_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype', 'categorical.phenotype'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-all-phenotypes'
                )
            );

        # check handling of only continuous phenotype (not supposed to be plotted)
        expect_no_error(
            create.pgs.density.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-continuous-phenotype'
                )
            );

        test.filename.all.phenotypes <- generate.filename(
            project.stem = 'TEST-all-phenotypes',
            file.core = 'pgs-density',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename.all.phenotypes))
            );

        test.filename.continuous.phenotype <- generate.filename(
            project.stem = 'TEST-continuous-phenotype',
            file.core = 'pgs-density',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename.continuous.phenotype))
            );


        # check handling of categorical phenotype with fewer than 2 samples in all categories
        low.n.test.data <- pgs.test;
        low.n.test.data$categorical.phenotype <- paste0('cat', 1:nrow(low.n.test.data));
        expect_warning(
            create.pgs.density.plot(
                pgs.data = low.n.test.data,
                phenotype.columns = c('categorical.phenotype'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-all-small-categories'
                )
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot correctly validates inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check that input data is a data frame
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = list(),
                phenotype.columns = 'continuous.phenotype'
                ),
            'pgs.data must be a data.frame'
            );

        # check that the required columns are present
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = data.frame(not.recognized.PGS.column = 1:10, continuous.phenotype = 1:10),
                phenotype.columns = 'continuous.phenotype'
                ),
            'No recognized PGS columns found in pgs.data'
            );

        # check that phenotype.column is a character vector
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = c(1,2,3)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.column is present in pgs.data
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'missing.phenotype'
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that phenotype.column does not contain recognized PGS columns
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'PGS.with.replaced.missing'
                ),
            'phenotype.columns cannot contain recognized PGS column names'
            );
        # check that at least one of the phenotype columns provided is a continuous variable
        expect_warning(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'binary.phenotype'
                ),
            'No continuous phenotype variables detected; returning NULL'
            );
        expect_equal(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'binary.phenotype'
                ),
            NULL
            );
        # handle NULL phenotype.columns
        expect_warning(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = NULL
                ),
            'No continuous phenotype variables detected; returning NULL'
            );
        expect_equal(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = NULL
                ),
            NULL
            );
        # check that output.dir is a real directory
        expect_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'continuous.phenotype',
                output.dir = 'not/a/real/directory'
                ),
            'not/a/real/directory does not exist'
            );
        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot runs with no error with basic inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing, PGS.with.normalized.missing)),
                phenotype.columns = 'continuous.phenotype',
                output.dir = temp.dir,
                filename.prefix = 'TEST'
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'pgs-scatter',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        # check returned object
        test.plot.object <- create.pgs.with.continuous.phenotype.plot(
            pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing, PGS.with.normalized.missing)),
            phenotype.columns = 'continuous.phenotype',
            filename.prefix = 'TEST',
            output.dir = NULL
            );
        expect_equal(
            class(test.plot.object),
            'multipanel'
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot runs correctly with correlation disabled', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'continuous.phenotype',
                output.dir = temp.dir,
                filename.prefix = 'TEST-no-correlation',
                compute.correlation = FALSE
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-no-correlation',
            file.core = 'pgs-scatter',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot runs correctly with tidy titles enabled', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'continuous.phenotype',
                output.dir = temp.dir,
                filename.prefix = 'TEST-tidy-titles',
                tidy.titles = TRUE
                )
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot runs correctly with include.origin enabled', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # add phenotype with non-zero minimum
        pgs.test$continuous.phenotype.shift <- pgs.test$continuous.phenotype + 10;

        # plot pgs with continuous phenotype
        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'continuous.phenotype.shift',
                output.dir = temp.dir,
                filename.prefix = 'TEST-include-origin',
                include.origin = TRUE
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-include-origin',
            file.core = 'pgs-scatter',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotype.plot runs correctly with multiple phenotypes', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # add another continuous phenotype
        pgs.test$continuous.phenotype2 <- rnorm(nrow(pgs.test));

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'continuous.phenotype2'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-two-continuous-phenotypes'
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-two-continuous-phenotypes',
            file.core = 'pgs-scatter',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.with.continuous.phenotypes.plot correctly switches to hexbin plot with default parameters', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # generate large test dataset
        large.pgs.test <- pgs.test[rep(1:nrow(pgs.test), 100), ];
        large.pgs.test$hexbin.phenotype <- rnorm(nrow(large.pgs.test));

        expect_no_error(
            create.pgs.with.continuous.phenotype.plot(
                pgs.data = large.pgs.test,
                phenotype.columns = c('hexbin.phenotype'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-hexbin',
                width = 15,
                height = 10
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-hexbin',
            file.core = 'pgs-scatter',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.rank.plot correctly validates inputs', {
        # check that input data is a data frame
        expect_error(
            create.pgs.rank.plot(
                pgs.data = list()
                ),
            'pgs.data must be a data frame'
            );

        # check that required columns are present
        expect_error(
            create.pgs.rank.plot(
                pgs.data = data.frame(sample = 1:10)
                ),
            'pgs.data must contain columns for Indiv, percentile, decile, quartile, n.missing.genotypes, percent.missing.genotypes'
            );
        # check that missing genotype style is correct
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                missing.genotype.style = 'not.a.style'
                ),
            'missing.genotype.style must be either "count" or "percent"'
            );
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                missing.genotype.style = c('count', 'percent')
                ),
            'missing.genotype.style must be either "count" or "percent"'
            );

        # check that phenotype.columns is a character vector
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = c(4,5,6)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.columns are present in data
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'missing.phenotype'
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that categorical.palette has valid colors
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                categorical.palette = c('red', 'blue', 'green', 'not.a.color')
                ),
            'categorical.palette must be a vector of valid colors'
            );
        # check that there are sufficient colors for the number of categories
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = 'categorical.phenotype',
                categorical.palette = c('red', 'blue')
                ),
            'Number of unique categories in a phenotype covariate exceeds the number of colors in the color palette. Please provide a larger color palette.'
            );
        # check that binary.palette has valid colors
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                binary.palette = c('red', 'blue', 'green', 'not.a.color')
                ),
            'binary.palette must be a vector of valid colors'
            );
        # check that there are sufficient colors for the number of phenotypes
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('binary.phenotype', 'continuous.phenotype'),
                binary.palette = c('red')
                ),
            'Number of binary and continuous phenotype covariates exceeds the number of binary color palettes. Please provide a larger color palette.'
            );
        # check that output.dir is a directory
        expect_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                output.dir = 'not/a/directory'
                ),
            'not/a/directory does not exist'
            );
        }
    );

test_that(
    'create.pgs.rank.plot runs with no error with basic inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check file writing
        temp.dir <- tempdir();

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                filename.prefix = 'TEST',
                output.dir = temp.dir
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        # check returned object
        test.plot.object <- create.pgs.rank.plot(
            pgs.data = pgs.test,
            filename.prefix = 'TEST',
            output.dir = NULL
            );
        expect_equal(
            class(test.plot.object),
            'multipanel'
            );

        }
    );

test_that(
    'create.pgs.rank.plot correctly switches between missing genotype barplot styles', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                filename.prefix = 'TEST-missing-genotype-percent',
                output.dir = temp.dir,
                missing.genotype.style = 'percent' # default is count
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-missing-genotype-percent',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );
        }
    );

test_that(
    'create.pgs.rank.plot correctly switches between large and small sample size modes', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

    # create test file with more than 50 individuals
    large.pgs.test <- pgs.test[rep(1:nrow(pgs.test), 6), ];
    large.pgs.test$Indiv <- paste0(large.pgs.test$Indiv, rep(1:6, each = nrow(pgs.test)));

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = large.pgs.test,
                filename.prefix = 'TEST-large-sample-size',
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype'),
                output.dir = temp.dir
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-large-sample-size',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.rank.plot correctly handles user-provided color palette', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('categorical.phenotype', 'binary.phenotype', 'continuous.phenotype'),
                filename.prefix = 'TEST-color-palette',
                output.dir = temp.dir,
                categorical.palette = c('red', 'orange','yellow', 'green', 'blue', 'purple'),
                binary.palette = c('red2', 'orange2')
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-color-palette',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.rank.plot runs correctly with user provided phenotypes',{
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype', 'categorical.phenotype'),
                filename.prefix = 'TEST-all-phenotypes',
                output.dir = temp.dir
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-all-phenotypes',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );

test_that(
    'create.pgs.rank.plot runs correctly when no missing genotypes are present', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS || SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        # remove missing genotypes
        pgs.no.missing.test <- pgs.test;
        pgs.no.missing.test$n.missing.genotypes <- 0;

        expect_no_error(
            create.pgs.rank.plot(
                pgs.data = pgs.no.missing.test,
                filename.prefix = 'TEST-no-missing',
                output.dir = temp.dir
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST-no-missing',
            file.core = 'pgs-rank-plot',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );
