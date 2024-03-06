test_that(
    'write.per.sample.pgs.table correctly checks inputs', {
        temp.dir <- tempdir();
        # check that input is a data.frame
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = 'not a data.frame',
                output.path = temp.dir
                ),
            'pgs.data must be a data.frame'
            );

        # check that output.path is a character string
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.path = 1
                )
            );

        # check that output.path is a directory
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.path = 'not/a/real/dir/'
                )
            );

        # check that file.prefix is a character string
        expect_error(
            write.per.sample.pgs.table(
                per.sample.pgs.summary.data = data.frame(),
                output.path = temp.dir,
                file.prefix = 1
                ),
            'file.prefix must be a character string'
            );
        }
    );

test_that(
    'write.per.sample.pgs.table correctly writes file', {
        temp.dir <- tempdir();
        test.pgs.per.sample <- data.frame(
            sample = c('sample1', 'sample2'),
            PGS = c(1, 3)
            );
        write.per.sample.pgs.table(
            per.sample.pgs.summary.data = test.pgs.per.sample,
            output.path = temp.dir
            );
        test.filename <- generate.filename(
            project.stem = NULL,
            file.core = 'per-sample-pgs-summary',
            extension = 'txt'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.filename))
            );

        # read file back in and check that it is correct
        read.in.pgs.per.sample <- read.table(
            file = file.path(temp.dir, test.filename),
            header = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
            );
        expect_equal(
            read.in.pgs.per.sample,
            test.pgs.per.sample
            );
        }
    );
