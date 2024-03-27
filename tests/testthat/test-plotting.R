test_that(
    'plot.pgs.density runs with no error', {
        temp.dir <- tempdir();
        # load test data
        load('data/phenotype.test.data.Rda');

        # run pgs statistics
        pgs.test <- apply.polygenic.score(
            vcf.data = phenotype.test.data$vcf.data,
            pgs.weight.data = phenotype.test.data$pgs.weight.data,
            phenotype.data = phenotype.test.data$phenotype.data
            );

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
