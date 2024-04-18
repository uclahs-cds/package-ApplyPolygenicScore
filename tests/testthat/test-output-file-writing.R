test_that(
    'write.per.sample.pgs.table correctly checks inputs', {
        temp.dir <- tempdir();
        # check that input is a data.frame
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = 'not a data.frame',
                output.dir = temp.dir
                ),
            'input data must be a data.frame'
            );

        # check that output.path is a character string
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.dir = 1
                )
            );

        # check that output.path is a directory
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.dir = 'not/a/real/dir/'
                )
            );

        # check that file.prefix is a character string
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.dir = temp.dir,
                file.prefix = 1
                ),
            'file.prefix must be a character string'
            );
        }
    );

test_that(
    'write.apply.polygenic.score.output.to.file correctly writes file', {
        temp.dir <- tempdir();
        test.pgs.output <- list(
            pgs.output = data.frame(
                sample = c('sample1', 'sample2'),
                PGS = c(1, 3)
                ),
            regression.output = data.frame(
                sample = c('sample1', 'sample2'),
                beta = c(0.1, 0.2),
                se = c(0.01, 0.02),
                p.value = c(0.001, 0.002)
                )
            );

        expect_no_error(
            write.apply.polygenic.score.output.to.file(
                apply.polygenic.score.output = test.pgs.output,
                output.dir = temp.dir,
                file.prefix = 'TEST'
                )
            );

        test.pgs.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'per-sample-pgs-summary',
            extension = 'txt'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.pgs.filename))
            );

        test.regression.filename <- generate.filename(
            project.stem = 'TEST',
            file.core = 'pgs-regression-output',
            extension = 'txt'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.regression.filename))
            );

        # read file back in and check that it is correct
        read.in.pgs.per.sample <- read.table(
            file = file.path(temp.dir, test.pgs.filename),
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
            );
        expect_equal(
            read.in.pgs.per.sample,
            test.pgs.output$pgs.output
            );

        read.in.pgs.regression <- read.table(
            file = file.path(temp.dir, test.regression.filename),
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
            );
        expect_equal(
            read.in.pgs.regression,
            test.pgs.output$regression.output
            );

        # check that when no regression data is present, regression file is not written
        test.pgs.output <- list(
            pgs.output = data.frame(
                sample = c('sample1', 'sample2'),
                PGS = c(1, 3)
                ),
            regression.output = NULL
            );

        write.apply.polygenic.score.output.to.file(
            apply.polygenic.score.output = test.pgs.output,
            output.dir = temp.dir,
            file.prefix = 'TEST-no-regression'
            );

        test.no.regression.filename <- generate.filename(
            project.stem = 'TEST-no-regression',
            file.core = 'pgs-regression-output',
            extension = 'txt'
            );

        test.no.regression.pgs.filename <- generate.filename(
            project.stem = 'TEST-no-regression',
            file.core = 'per-sample-pgs-summary',
            extension = 'txt'
            );

        expect_false(
            file.exists(file.path(temp.dir, test.no.regression.filename))
            );
        expect_true(
            file.exists(file.path(temp.dir, test.no.regression.pgs.filename))
            );

        }
    );
