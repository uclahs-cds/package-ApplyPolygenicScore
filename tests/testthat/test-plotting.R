# plotting functions take a long time to run, this var toggles off plotting tests for faster testing
SKIP.PLOTS <- TRUE;
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
    missing.genotype.method = c('mean.dosage', 'normalize'),
    n.percentiles = 2
    )$pgs.output;
# add missing genotpye counts
pgs.test$n.missing.genotypes <- sample(1:10, nrow(pgs.test), replace = TRUE);
# add some missing data
pgs.test$continuous.phenotype[1:2] <- NA;
pgs.test$binary.phenotype[2:3] <- NA;
pgs.test$categorical.phenotype[3:4] <- NA;
pgs.test$decile[4:5] <- NA;
pgs.test$percentile[4:5] <- NA;

test_that(
    'plot.pgs.density correctly validates inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check that input data is a data frame
        expect_error(
            plot.pgs.density(
                pgs.data = list()
                ),
            'pgs.data must be a data.frame'
            );

        # check that the required columns are present
        expect_error(
            plot.pgs.density(
                pgs.data = data.frame(not.recognized.PGS.colum = 1:10)
                ),
            'No recognized PGS columns found in pgs.data'
            );

        # check that phenotype.columns is a character vector
        expect_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                phenotype.columns = c(1,2,3)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.columns are present in pgs.data
        expect_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                phenotype.columns = c('missing.phenotype')
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that phenotype.columns do not contain recognized PGS columns
        expect_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                phenotype.columns = c('PGS.with.replaced.missing')
                ),
            'phenotype.columns cannot contain recognized PGS column names'
            );
        # check that output.dir is a real directory
        expect_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                output.dir = 'not/a/real/directory'
                ),
            'not/a/real/directory does not exist'
            );
        }
    );

test_that(
    'plot.pgs.density runs with no error with basic inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs density
        expect_no_error(
            plot.pgs.density(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing)),
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
        test.plot.object <- plot.pgs.density(
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
    'plot.pgs.density runs correctly with tidy titles enabled', {
        skip.plotting.tests(skip.plots = SKIP.COMPREHENSIVE.CASES || SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs density
        expect_no_error(
            plot.pgs.density(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing)),
                phenotype.columns = NULL,
                output.dir = getwd(),#temp.dir,
                filename.prefix = 'TEST-tidy-titles',
                tidy.titles = TRUE
                )
            );

        }
    );


test_that(
    'plot.pgs.density runs correctily with user provided phenotypes', {
        skip.plotting.tests(skip.plots = SKIP.COMPREHENSIVE.CASES || SKIP.PLOTS);

        temp.dir <- tempdir();

        # check handling of many phenotypes
        expect_no_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype', 'categorical.phenotype'),
                output.dir = temp.dir,
                filename.prefix = 'TEST-all-phenotypes'
                )
            );

        # check handling of only continuous phenotype (not supposed to be plotted, how can I test for multipanel dimensions over a multipanel object?)
        expect_no_error(
            plot.pgs.density(
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

        }
    );

test_that(
    'plot.pgs.with.continuous.phenotype correctly validates inputs', {
        #skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check that input data is a data frame
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = list(),
                phenotype.columns = 'continuous.phenotype'
                ),
            'pgs.data must be a data.frame'
            );

        # check that the required columns are present
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = data.frame(not.recognized.PGS.column = 1:10, continuous.phenotype = 1:10),
                phenotype.columns = 'continuous.phenotype'
                ),
            'No recognized PGS columns found in pgs.data'
            );

        # check that phenotype.column is a character vector
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = c(1,2,3)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.column is present in pgs.data
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = 'missing.phenotype'
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that phenotype.column does not contain recognized PGS columns
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = 'PGS.with.replaced.missing'
                ),
            'phenotype.columns cannot contain recognized PGS column names'
            );
        # check that at least one of the phenotype columns provided is a continuous variable
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = 'binary.phenotype'
                ),
            'No continuous phenotype variables detected'
            );
        # check that output.dir is a real directory
        expect_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = 'continuous.phenotype',
                output.dir = 'not/a/real/directory'
                ),
            'not/a/real/directory does not exist'
            );
        }
    );

test_that(
    'plot.pgs.with.continuous.phenotype runs with no error with basic inputs', {
        #skip.plotting.tests(skip.plots = SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = subset(pgs.test, select = -c(PGS.with.replaced.missing)),
                phenotype.columns = 'continuous.phenotype',
                output.dir = getwd(),#temp.dir,
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
        test.plot.object <- plot.pgs.with.continuous.phenotype(
            pgs.data = pgs.test,
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
    'plot.pgs.with.continuous.phenotype runs correctly with multiple phenotypes', {
        #skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # add another continuous phenotype
        pgs.test$continuous.phenotype2 <- rnorm(nrow(pgs.test));

        temp.dir <- tempdir();

        # plot pgs with continuous phenotype
        expect_no_error(
            plot.pgs.with.continuous.phenotype(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'continuous.phenotype2'),
                output.dir = getwd(),#temp.dir,
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