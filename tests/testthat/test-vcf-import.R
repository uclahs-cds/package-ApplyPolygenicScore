test_that(
    'import.vcf successfully imports VCF', {
        # check that correct input is accepted
        expect_no_error(
            import.vcf(
                input = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
                info.fields = NULL,
                format.fields = NULL
                )
            );
        }
    );

test_that(
    'import.vcf outputs a vcfR tidy object', {
        test.vcf <- import.vcf(
            input = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = NULL,
            format.fields = NULL
            );

        expect_equal(
            class(test.vcf),
            'list'
            );

        expect_equal(
            names(test.vcf),
            c('dat', 'meta')
            );

        expect_s3_class(
            test.vcf$dat,
            'data.frame'
            );

        expect_s3_class(
            test.vcf$meta,
            'data.frame'
            );

        }
    );

test_that(
    'import.vcf accurately imports VCF fields', {
        test.vcf <- import.vcf(
            input = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = NULL,
            format.fields = NULL
            );

        # check that the correct number of variants were imported
        n.variants <- 63;
        n.samples <- 1;
        expect_equal(
            nrow(test.vcf$dat),
            n.variants * n.samples
            );

        # check that standard VCF columns were imported
        expect_true(
            all(c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER') %in% colnames(test.vcf$dat))
            );

        # check that required vcfR columns were added (sample names, genotype, and genotype alleles)
        expect_true(
            all(c('Indiv', 'gt_GT', 'gt_GT_alleles') %in% colnames(test.vcf$dat))
            );

        # check that the sample name is correct
        expect_equal(
            unique(test.vcf$dat$Indiv),
            'HG001'
            );
        
        # check the first three variants for correct coordiinates
        expect_equal(
            test.vcf$dat$CHROM[1:3],
            c('chr1', 'chr1', 'chr1')
            );
        expect_equal(
            test.vcf$dat$POS[1:3],
            c(87734095, 111721652, 154895579)
            );

        # check the first three variants for correct alleles
        expect_equal(
            test.vcf$dat$gt_GT[1:3],
            c('1/1', '0/1', '0/1')
            );
        expect_equal(
            test.vcf$dat$gt_GT_alleles[1:3],
            c('A/A', 'T/C', 'G/T')
            );

        }
    );

