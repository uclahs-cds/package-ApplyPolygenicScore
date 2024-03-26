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

        # check that missing genotype method input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = 'not a valid method'
                ),
            'missing.genotype.method must be either "mean.dosage", "normalize", or "none"'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('mean.dosage', 'normalize', 'not a valid method')
                ),
            'missing.genotype.method must be either "mean.dosage", "normalize", or "none"'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('mean.dosage', 'none')
                ),
            'If "none" is included in missing.genotype.method, it must be the only method included'
            );

        # check that percentile.source input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                percentile.source = c('mean.dosage', 'normalize')
                ),
            'percentile.source must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                percentile.source = 'not a valid method'
                ),
            'percentile.source must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = 'normalize',
                percentile.source = 'mean.dosage'
                ),
            'percentile.source must be one of the chosen missing genotype methods'
            );
        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                percentile.source = 'mean.dosage'
                )
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

        # check for effect allele frequency column when use.external.effect.allele.frequency option is selected
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                use.external.effect.allele.frequency = TRUE
                ),
            'pgs.weight.data must contain a column named allelefrequency_effect if use.external.effect.allele.frequency is TRUE'
            );

        # check for duplicate coordinates in PGS data
        duplicate.row <- test.pgs.weight.data$pgs.weight.data[1, ];
        duplicate.row.as.multiallelic <- duplicate.row;
        duplicate.row.as.multiallelic$effect_allele <- 'C';
        test.pgs.weight.data.duplicated.coordinates <- rbind(test.pgs.weight.data$pgs.weight.data, duplicate.row.as.multiallelic);
        expect_warning(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data.duplicated.coordinates
                ),
            'Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.'
            );

        # check for duplicate variants in PGS data
        test.pgs.weight.data.duplicated.variants <- rbind(test.pgs.weight.data$pgs.weight.data, duplicate.row);
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data.duplicated.variants
                ),
            'Duplicate variants detected in the PGS weight data. Please ensure only unique coordinate:effect allele combinations are present.'
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
            5
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
        # test case with VCF multiallelic site with no extra beta, REF allele is the risk allele
        ref.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.as.single.risk.allele.multiallelic.pgs.weight.data
            );
        expect_equal(
            ref.as.single.risk.allele.test$PGS,
            c(4, 3, 0)
            );

        # test case with VCF multiallelic sites with no extra beta, ALT allele is the risk allele
        alt.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data
            );
        expect_equal(
            alt.as.single.risk.allele.test$PGS,
            c(2, 1, 2)
            );
        expect_no_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data
                )
            );

        # test case with a VCF multiallelic site that also has an extra beta, both ALT alleles are risk alleles
        alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data
            );
        expect_equal(
            alt.as.two.risk.alleles.test$PGS,
            c(2, 1.5, 3)
            );
        expect_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data
                ),
            'Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.'
            );

        # test case with a VCF multiallelic site that also has an extra beta, one REF and one ALT allele are risk alleles
        ref.and.alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data
            );
        expect_equal(
            ref.and.alt.as.two.risk.alleles.test$PGS,
            c(2, 1.5, 2)
            );
        expect_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data
                ),
            'Multiple effect alleles found in sample1 genotype, choosing effect allele with highest beta for dosage calculation. Check coordinates chr2:2'
            );
        }
    );

test_that(
    'apply.polygenic.score correctly handles missing genotypes', {
        load('data/missing.genotype.test.data.Rda');
        test.missing.genotype.mean.dosage <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage'
            );
        test.missing.genotype.normalize <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'normalize'
            );
        test.missing.genotype.both <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize')
            );
        test.missing.genotype.none <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'none'
            );

        # check column names
        percentile.colnames <- c('percentile', 'decile', 'quartile');
        expect_equal(
            colnames(test.missing.genotype.mean.dosage),
            c('sample', 'PGS.with.replaced.missing', percentile.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.normalize),
            c('sample', 'PGS.with.normalized.missing', percentile.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.both),
            c('sample', 'PGS.with.normalized.missing', 'PGS.with.replaced.missing', percentile.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.none),
            c('sample', 'PGS', percentile.colnames)
            );

        expect_equal(
            test.missing.genotype.mean.dosage$PGS.with.replaced.missing,
            c(1, 4, 2.5, 2.5)
            );
        expect_equal(
            test.missing.genotype.normalize$PGS.with.normalized.missing,
            c(1 / 6, 4 / 6, NA, 0.5)
            );
        expect_equal(
            test.missing.genotype.both$PGS.with.normalized.missing,
            c(1 / 6, 4 / 6, NA, 0.5)
            );
        expect_equal(
            test.missing.genotype.both$PGS.with.replaced.missing,
            c(1, 4, 2.5, 2.5)
            );
        expect_equal(
            test.missing.genotype.none$PGS,
            c(1, 4, 0, 2)
            );
        }
    );

test_that(
    'apply.polygenic.score correctly handles external effect allele frequency', {
        load('data/missing.genotype.test.data.Rda');
        # add effect allele frequency column to PGS weight data
        missing.genotype.test.data$missing.genotype.pgs.weight.data$allelefrequency_effect <- c(0.5, 0.25, 0.75, 0.5);
        test.missing.genotype.mean.dosage <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage',
            use.external.effect.allele.frequency = TRUE
            );

        expect_equal(
            test.missing.genotype.mean.dosage$PGS.with.replaced.missing,
            c(1, 4, 3, 3.5)
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

        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                n.percentiles = 5
                )
            )

        }
    );
