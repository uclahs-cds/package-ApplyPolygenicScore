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

        # Verify correct number of variants and samples
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat[1:3, ],
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.'
            )
    }
)

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
