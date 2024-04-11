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
