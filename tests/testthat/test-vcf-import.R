test_that(
    'import.vcf successfully imports VCF', {
        # check that correct input is accepted
        expect_no_error(
            import.vcf(
                vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
                info.fields = NULL,
                format.fields = NULL
                )
            );
        expect_no_warning(
            import.vcf(
                vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
                info.fields = NULL,
                format.fields = NULL
                )
            );
        # turn off verbose output
        expect_silent(
            import.vcf(
                vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
                info.fields = NULL,
                format.fields = NULL,
                verbose = FALSE
                )
            );
        }
    );

test_that(
    'import.vcf outputs a vcfR tidy object', {
        test.vcf <- import.vcf(
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
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
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
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

test_that(
    'import.vcf accurately imports INFO fields', {
        all.info.test.vcf <- import.vcf(
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = NULL,
            format.fields = NULL
            );

        select.info.test.vcf <- import.vcf(
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = c('DPSum', 'platforms', 'arbitrated'),
            format.fields = NULL
            );

        # check that all requested INFO columns are present
        expect_true(
            all(c('DPSum', 'platforms', 'arbitrated') %in% colnames(select.info.test.vcf$dat))
            );
        expect_true(
            all(c('DPSum', 'platforms', 'platformnames', 'platformbias', 'datasets','datasetnames', 'datasetsmissingcall', 'callsets', 'callsetnames', 'varType','filt', 'callable', 'difficultregion', 'arbitrated', 'callsetwiththisuniqgenopassing', 'callsetwithotheruniqgenopassing') %in% colnames(all.info.test.vcf$dat))
            );

        # check that the total number of columns is correct
        generic.vcf.ncol <- 7;
        generic.vcfR.ncol <- 2;
        info.ncol <- 16;
        select.info.ncol <- 3;
        format.ncol <- 6;
        ncol.all <- generic.vcf.ncol + generic.vcfR.ncol + info.ncol + format.ncol;
        ncol.select <- generic.vcf.ncol + generic.vcfR.ncol + select.info.ncol + format.ncol;

        expect_equal(
            ncol(select.info.test.vcf$dat),
            ncol.select
            );
        expect_equal(
            ncol(all.info.test.vcf$dat),
            ncol.all
            );

        # check that one of the selected INFO columns has the correct value
        expect_equal(
            select.info.test.vcf$dat$platforms[1:3],
            c(5, 5, 5)
            );
        expect_equal(
            all.info.test.vcf$dat$platforms[1:3],
            c(5, 5, 5)
            );

        }
    );

test_that(
    'import.vcf accurately imports FORMAT fields', {
        all.format.test.vcf <- import.vcf(
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = NULL,
            format.fields = NULL
            );

        select.format.test.vcf <- import.vcf(
            vcf.path = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            info.fields = NULL,
            format.fields = c('DP', 'GQ', 'GT')
            );

        # check that all requested FORMAT columns are present
        expect_true(
            all(c('gt_DP', 'gt_GQ', 'gt_GT') %in% colnames(select.format.test.vcf$dat))
            );
        expect_true(
            all(c('gt_DP', 'gt_GQ', 'gt_ADALL', 'gt_AD', 'gt_GT', 'gt_PS') %in% colnames(all.format.test.vcf$dat))
            );

        # check that one of the selected FORMAT columns has the correct value
        expect_equal(
            select.format.test.vcf$dat$gt_DP[1:3],
            c(762, 658, 966)
            );
        expect_equal(
            all.format.test.vcf$dat$gt_DP[1:3],
            c(762, 658, 966)
            );

        }
    );

test_that(
    'check.for.no.info.fields utility correctly catches empty INFO bug', {
        test.vcf.yes.info <- vcfR::read.vcfR(
            file = 'data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10.vcf.gz',
            convertNA = TRUE,
            verbose = FALSE
            );

        test.vcf.no.info <- test.vcf.yes.info;
        test.vcf.no.info@meta <- test.vcf.no.info@meta[-grep('^##INFO', test.vcf.no.info@meta)];

        expect_silent(
            check.for.no.info.fields(test.vcf.yes.info)
            );

        expect_warning(
            check.for.no.info.fields(test.vcf.no.info),
            'No INFO fields present in VCF file. Adding dummy INFO field "DUMMY" to avoid vcfR import errors.'
            );
        }
    );
