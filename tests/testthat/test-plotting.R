# plotting functions take a long time to run, this var toggles off plotting tests for faster testing
SKIP.PLOTS <- FALSE;
SKIP.COMPREHENSIVE.CASES <- FALSE;
skip.plotting.tests <- function(skip.plots = FALSE) {
    # skip plotting tests if 
    if (skip.plots) {
        skip('Plotting tests disabled');
        }
    }


# get test dataset
load('data/phenotype.test.data.Rda');

# run pgs statistics
pgs.test <- apply.polygenic.score(
    vcf.data = phenotype.test.data$vcf.data,
    pgs.weight.data = phenotype.test.data$pgs.weight.data,
    phenotype.data = phenotype.test.data$phenotype.data,
    missing.genotype.method = c('mean.dosage', 'normalize'),
    n.percentiles = 2
    );
# add missing genotpye counts
pgs.test$n.missing.genotypes <- sample(1:10, nrow(pgs.test), replace = TRUE);
# add some missing data
pgs.test$continuous.phenotype[1:2] <- NA;
pgs.test$binary.phenotype[2:3] <- NA;
pgs.test$categorical.phenotype[3:4] <- NA;
pgs.test$decile[4:5] <- NA;
pgs.test$percentile[5:6] <- NA;


test_that(
    'plot.pgs.density runs with no error', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);
        skip();

        temp.dir <- tempdir();

        # plot pgs density
        expect_no_error(
            plot.pgs.density(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype'),
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

        }
    );

test_that(
    'plot.pgs.rank correctly validates inputs', {
        # check that input data is a data frame
        expect_error(
            plot.pgs.rank(
                pgs.data = list()
                ),
            'pgs.data must be a data frame'
            );

        # check that required columns are present
        expect_error(
            plot.pgs.rank(
                pgs.data = data.frame(sample = 1:10)
                ),
            'pgs.data must contain columns for Indiv, percentile, decile, quartile, and n.missing.genotypes'
            );
        # check that phenotype.columns is a character vector
        expect_error(
            plot.pgs.rank(
                pgs.data = pgs.test,
                phenotype.columns = c(4,5,6)
                ),
            'phenotype.columns must be a character vector'
            );
        # check that phenotype.columns are present in data
        expect_error(
            plot.pgs.rank(
                pgs.data = pgs.test,
                phenotype.columns = 'missing.phenotype'
                ),
            'phenotype.columns must be a subset of the column names in pgs.data'
            );
        # check that output.dir is a directory
        expect_error(
            plot.pgs.rank(
                pgs.data = pgs.test,
                output.dir = 'not/a/directory'
                ),
            'not/a/directory does not exist'
            );
        }
    );

test_that(
    'plot.pgs.rank runs with no error with basic inputs', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        # check file writing
        temp.dir <- tempdir();

        expect_no_error(
            plot.pgs.rank(
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
        test.plot.object <- plot.pgs.rank(
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
    'plot.pgs.rank runs correctly with user provided phenotypes',{
        skip.plotting.tests(skip.plots = SKIP.COMPREHENSIVE.CASES);

        temp.dir <- tempdir();

        expect_no_error(
            plot.pgs.rank(
                pgs.data = pgs.test,
                phenotype.columns = c('continuous.phenotype', 'binary.phenotype', 'categorical.phenotype'),
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

        }
    );
