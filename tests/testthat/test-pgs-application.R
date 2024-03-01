test_that(
    'apply.polygenic.score correctly checks inputs', {
        test.vcf.data <- import.vcf('data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10_duplicated-sample.vcf.gz')
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');

        # check that only data frame inputs are accepted
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'vcf.data must be a data.frame'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data
                ),
            'pgs.weight.data must be a data.frame'
            );

        # check that required columns are present
        test.vcf.data.missing.columns <- test.vcf.data$dat;
        test.vcf.data.missing.columns$gt_GT_alleles <- NULL;
        test.pgs.weight.data.missing.columns <- test.pgs.weight.data$pgs.weight.data;
        test.pgs.weight.data.missing.columns$beta <- NULL;

        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data.missing.columns,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'vcf.data must contain columns named CHROM, POS, REF, ALT, Indiv, and gt_GT_alleles'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data.missing.columns
                ),
            'pgs.weight.data must contain columns named CHROM, POS, effect_allele, and beta'
            );

        # check for duplicate variants in PGS data
        test.pgs.weight.data.duplicated.variants <- rbind(test.pgs.weight.data$pgs.weight.data, test.pgs.weight.data$pgs.weight.data[1, ]);
        expect_warning(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data.duplicated.variants
                ),
            'Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.'
            );


        # Verify correct number of variants and samples
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat[1:3, ],
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.'
            );
        }
    );

test_that(
    'apply.polygenic.score correctly formats outputs', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        # check that output is a data.frame
        expect_s3_class(
            test.pgs.per.sample,
            'data.frame'
            );

        # check that output has correct number of rows and columns
        expect_equal(
            nrow(test.pgs.per.sample),
            2
            );
        expect_equal(
            ncol(test.pgs.per.sample),
            2
            );
        }
    );

test_that(
    'apply.polygenic.score correctly calculates pgs', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        # check that output is correct
        expect_equal(
            test.pgs.per.sample$sample,
            c('sample1', 'sample2')
            );
        expect_equal(
            test.pgs.per.sample$PGS,
            c(1, 3)
            );
        }
    );

test_that(
    'apply.polygenic.score correctly handles multiallic sites', {
        load('data/merged.multiallelic.site.test.data.Rda');
        # test that 
        ref.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.as.single.risk.allele.multiallelic.pgs.weight.data
            );
        expect_equal(
            ref.as.single.risk.allele.test$PGS,
            c(4, 3, 0)
            );

        alt.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data
            );
        expect_equal(
            alt.as.single.risk.allele.test$PGS,
            c(2, 1, 2)
            );

        alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data
            );
        expect_equal(
            alt.as.two.risk.alleles.test$PGS,
            c(2, 1.5, 3)
            );
    
        ref.and.alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data
            );
        expect_equal(
            ref.and.alt.as.two.risk.alleles.test$PGS,
            c(2, 1.5, 2)
            );
        }
    );

test_that(
    'apply.polygenic.score works correctly on real data', {
        test.vcf.data <- import.vcf('data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10_duplicated-sample.vcf.gz')
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');

        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            )

        }
    );
