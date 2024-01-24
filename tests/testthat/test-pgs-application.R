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
