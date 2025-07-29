## get.variant.by.sample.matrix tests ##
test_that(
    'get.variant.by.sample.matrix correctly transforms long data to SNP by sample matrix', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample2', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.variant.by.sample.matrix <- get.variant.by.sample.matrix(
            long.data = test.long.data,
            row.id.cols = 'variant.id',
            value.col = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(2, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );

test_that(
    'get.variant.by.sample.matrix utility correctly handles duplicates', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        expect_error(
            get.variant.by.sample.matrix(
                long.data = test.long.data,
                row.id.cols = 'variant.id',
                value.col = 'dosage'
                )
            );
        }
    );

## get.combined.multiallelic.variant.by.sample.matrix tests ##
test_that(
    'get.combined.multiallelic.variant.by.sample.matrix correctly transforms long data to SNP by sample matrix', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample2', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.variant.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data,
            variant.id = test.long.data$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(2, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );

test_that(
    'get.combined.multiallelic.variant.by.sample.matrix utility correctly handles duplicates', {
        test.long.data <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, 2, 3, 4, NA)
            );
        test.long.data.one.missing.duplicate <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(1, NA, 3, 4, NA)
            );
        test.long.data.two.missing.duplicates <- data.frame(
            variant.id = c('variant1', 'variant1', 'variant2', 'variant2', 'variant3'),
            Indiv = c('sample1', 'sample1', 'sample1', 'sample2', NA),
            dosage = c(NA, NA, 3, 4, NA)
            );

        test.variant.by.sample.matrix <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data,
            variant.id = test.long.data$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix,
            as.matrix(data.frame(
                sample1 = c(3 / 2, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );

        test.variant.by.sample.matrix.one.missing.duplicate <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data.one.missing.duplicate,
            variant.id = test.long.data.one.missing.duplicate$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix.one.missing.duplicate,
            as.matrix(data.frame(
                sample1 = c(1, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );

        test.variant.by.sample.matrix.two.missing.duplicates <- get.combined.multiallelic.variant.by.sample.matrix(
            long.data = test.long.data.two.missing.duplicates,
            variant.id = test.long.data.two.missing.duplicates$variant.id,
            value.var = 'dosage'
            );
        expect_equal(
            test.variant.by.sample.matrix.two.missing.duplicates,
            as.matrix(data.frame(
                sample1 = c(NA, 3, NA),
                sample2 = c(NA, 4, NA),
                row.names = c('variant1', 'variant2', 'variant3')
                ))
            );
        }
    );

## convert long vcf to wide vcf tests ##
test_that(
    'convert.long.vcf.to.wide.vcf correctly transforms long VCF data to wide format', {
        # load simple test VCF data
        load('data/simple.pgs.application.test.data.Rda');

        expect_silent(
            converted.vcf <- convert.long.vcf.to.wide.vcf(long.vcf = simple.pgs.application.test.data$vcf.data)
            );

        expect_equal(
            length(converted.vcf),
            2
            );

        expect_equal(
            names(converted.vcf),
            c('genotyped.alleles', 'vcf.fixed.fields')
            );

        # check data types of outputed objects
        expect_true(
            is.matrix(converted.vcf$genotyped.alleles)
            );
        expect_true(
            is.data.frame(converted.vcf$vcf.fixed.fields)
            );

        # check dimensions of genotyped alleles matrix
        expect_equal(
            dim(converted.vcf$genotyped.alleles),
            c(2, 2) # 2 variants, 2 samples
            );

        # check dimensions of fixed fields data frame
        expect_equal(
            dim(converted.vcf$vcf.fixed.fields),
            c(2, 7) # 2 variants, 7 fixed fields
            );

        # check row names of genotyped alleles matrix
        expect_equal(
            rownames(converted.vcf$genotyped.alleles),
            unique(paste0(
                converted.vcf$vcf.fixed.fields$CHROM, ':',
                converted.vcf$vcf.fixed.fields$POS
                ))
            );

        # check column names of genotyped alleles matrix
        expect_equal(
            colnames(converted.vcf$genotyped.alleles),
            unique(simple.pgs.application.test.data$vcf.data$Indiv)
            );

        # check fixed fields data frame
        expect_equal(
            colnames(converted.vcf$vcf.fixed.fields),
            c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'variant.id', 'allele.matrix.row.index')
            );

        expected.fixed.fields <- data.frame(
            CHROM = c('chr1', 'chr2'),
            POS = c(1, 2),
            ID = c('rs1', 'rs2'),
            REF = c('A', 'T'),
            ALT = c('T', 'A'),
            variant.id = c('chr1:1', 'chr2:2'),
            allele.matrix.row.index = c(1, 2)
            );
        expected.fixed.fields <- data.table::as.data.table(expected.fixed.fields);

        expected.genotype.matrix <- matrix(
            data = c('A/A', 'T/T', 'A/T', 'T/A'),
            nrow = 2,
            ncol = 2,
            byrow = TRUE,
            dimnames = list(
                c('chr1:1', 'chr2:2'),
                c('sample1', 'sample2')
                )
            );

        expect_equal(
            converted.vcf$vcf.fixed.fields,
            expected.fixed.fields
            );

        expect_equal(
            converted.vcf$genotyped.alleles,
            expected.genotype.matrix
            );

        }
    );

test_that(
    'convert.long.vcf.to.wide.vcf correctly handles VCFs with multiallelic sites', {
        # load test VCF data with multiallelic sites
        load('data/merged.multiallelic.site.test.data.Rda');

        expect_silent(
            converted.vcf <- convert.long.vcf.to.wide.vcf(long.vcf = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data)
            );

        expected.fixed.fields <- data.frame(
            CHROM = c('chr1', 'chr2', 'chr3'),
            POS = c(1, 2, 3),
            ID = c('rs1', 'rs2', 'rs3'),
            REF = c('T', 'T', 'T'),
            ALT = c('A', 'A,C', 'A,G'),
            variant.id = c('chr1:1', 'chr2:2', 'chr3:3'),
            allele.matrix.row.index = c(1, 2, 3)
            );

        expected.matrix <- matrix(
            data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data$gt_GT_alleles,
            nrow = 3,
            ncol = 4,
            byrow = TRUE,
            dimnames = list(
                c('chr1:1', 'chr2:2', 'chr3:3'),
                c('sample1', 'sample2', 'sample3', 'sample4')
                )
            );

        expect_equal(
            converted.vcf$vcf.fixed.fields,
            data.table::as.data.table(expected.fixed.fields)
            );

        expect_equal(
            converted.vcf$genotyped.alleles,
            expected.matrix
            );

        }
    );

test_that(
    'convert.long.vcf.to.wide.vcf correctly handles split multiallelic sites', {
        load('data/split.multiallelic.site.test.data.Rda');
        split.multiallelic.site.test.data$split.multiallelic.vcf.data$ID <- c(rep('rs1', 2), rep('rs2', 4));

        expect_error(
            convert.long.vcf.to.wide.vcf(long.vcf = split.multiallelic.site.test.data$split.multiallelic.vcf.data),
            regexp = 'chr2_2_sample1.*chr2_2_sample2'
            );

        }
    );
