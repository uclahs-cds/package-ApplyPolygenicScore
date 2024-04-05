# plotting functions take a long time to run, this var toggles off plotting tests for faster testing
SKIP.PLOTS <- TRUE#FALSE;
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
        }
    );

test_that(
    'plot.pgs.rank runs with no error', {
        skip.plotting.tests(skip.plots = SKIP.PLOTS);

        temp.dir <- tempdir();

        # plot pgs rank
        expect_no_error(
            plot.pgs.rank(
                pgs.data = pgs.test,
                phenotype.columns = 'categorical.phenotype',#c('continuous.phenotype','binary.phenotype', 'binary.factor.phenotype', 'categorical.phenotype'),
                #output.dir = temp.dir,
                filename.prefix = 'TEST'
                )
            );

        test.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'pgs-rank',
            extension = 'png'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        }
    );
