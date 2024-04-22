PGS.OUTPUT.INDEX <- 1;
REGRESSION.OUTPUT.INDEX <- 2;
test_that(
    'apply.polygenic.score correctly checks inputs', {
        test.vcf.data <- import.vcf('data/HG001_GIAB.vcf.gz')
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');
        test.phenotype.data <- data.frame(Indiv = c('HG001', '2:HG001'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1));

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

        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('none', 'normalize', 'mean.dosage')
                )
            );

        # check that analysis.source.pgs input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = c('mean.dosage', 'normalize')
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = 'not a valid method'
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = 'normalize',
                analysis.source.pgs = 'mean.dosage'
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = 'mean.dosage'
                )
            );

        # check that phenotype data input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = 'not a data frame'
                ),
            'phenotype.data must be a data.frame'
            );

        # check for matching samples between phenotype and vcf data
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = data.frame(Indiv = c('sample11'))
                ),
            'No matching Indiv between phenotype.data and vcf.data'
            );

        # check required columns in phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = subset(test.phenotype.data, select = -Indiv)
                ),
            'phenotype.data must contain columns named Indiv'
            );

        # check for correct phenotype analysis columns
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = test.phenotype.data,
                phenotype.analysis.columns = c('not a valid column')
                ),
            'phenotype.analysis.columns must be columns in phenotype.data'
            );

        # check for missing phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.analysis.columns = 'continuous.phenotype'
                ),
            'phenotype.analysis.columns provided but no phenotype data detected'
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
    'apply.polygenic.score correctly formats general output', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        # check that output is a list
        expect_equal(
            class(test.pgs.per.sample),
            'list'
            );

        # check that output has correct number of elements
        expect_equal(
            length(test.pgs.per.sample),
            2
            );

        # check that output has correct names
        output.names <- c('pgs.output', 'regression.output');
        expect_equal(
            names(test.pgs.per.sample),
            output.names
            );

        }
    );

test_that(
    'apply.polygenic.score correctly formats pgs output', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        # check that output is a data.frame
        expect_s3_class(
            test.pgs.per.sample[[PGS.OUTPUT.INDEX]],
            'data.frame'
            );

        # check that output has correct number of rows and columns
        expect_equal(
            nrow(test.pgs.per.sample[[PGS.OUTPUT.INDEX]]),
            2
            );
        expect_equal(
            ncol(test.pgs.per.sample[[PGS.OUTPUT.INDEX]]),
            7
            );
        }
    );

test_that(
    'apply.polygenic.score correctly formats regression output', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample.with.phenotype <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
            phenotype.data = data.frame(Indiv = c('sample1', 'sample2'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1)),
            phenotype.analysis.columns = c('continuous.phenotype', 'binary.phenotype')
            );

        # check that output is a data.frame
        expect_s3_class(
            test.pgs.per.sample.with.phenotype[[REGRESSION.OUTPUT.INDEX]],
            'data.frame'
            );

        # check that output has correct number of rows and columns
        expect_equal(
            nrow(test.pgs.per.sample.with.phenotype[[REGRESSION.OUTPUT.INDEX]]),
            2
            );
        expect_equal(
            ncol(test.pgs.per.sample.with.phenotype[[REGRESSION.OUTPUT.INDEX]]),
            7
            );

        # check NULL output when no phenotype data is provided
        test.pgs.per.sample.no.phenotype <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        expect_equal(
            test.pgs.per.sample.no.phenotype[[REGRESSION.OUTPUT.INDEX]],
            NULL
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
            test.pgs.per.sample[[PGS.OUTPUT.INDEX]]$Indiv,
            c('sample1', 'sample2')
            );
        expect_equal(
            test.pgs.per.sample[[PGS.OUTPUT.INDEX]]$PGS,
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
            ref.as.single.risk.allele.test[[PGS.OUTPUT.INDEX]]$PGS,
            c(4, 3, 0)
            );

        # test case with VCF multiallelic sites with no extra beta, ALT allele is the risk allele
        alt.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data
            );
        expect_equal(
            alt.as.single.risk.allele.test[[PGS.OUTPUT.INDEX]]$PGS,
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
            alt.as.two.risk.alleles.test[[PGS.OUTPUT.INDEX]]$PGS,
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
            ref.and.alt.as.two.risk.alleles.test[[PGS.OUTPUT.INDEX]]$PGS,
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
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'mean.dosage'
            );
        test.missing.genotype.both.percentile.check <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'normalize'
            );
        test.missing.genotype.none <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'none'
            );

        # check column names
        percentile.colnames <- c('percentile', 'decile', 'quartile');
        extra.colnames <- c('n.missing.genotypes', 'percent.missing.genotypes');
        expect_equal(
            colnames(test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]),
            c('Indiv', 'PGS.with.replaced.missing', percentile.colnames, extra.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]),
            c('Indiv', 'PGS.with.normalized.missing', percentile.colnames, extra.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.both[[PGS.OUTPUT.INDEX]]),
            c('Indiv', 'PGS.with.normalized.missing', 'PGS.with.replaced.missing', percentile.colnames, extra.colnames)
            );
        expect_equal(
            colnames(test.missing.genotype.none[[PGS.OUTPUT.INDEX]]),
            c('Indiv', 'PGS', percentile.colnames, extra.colnames)
            );

        # check that PGS values are calculated correctly
        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$PGS.with.replaced.missing,
            c(1, 4, 2.5, 2.5)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$PGS.with.normalized.missing,
            c(1 / 6, 4 / 6, NA, 0.5)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$PGS.with.normalized.missing,
            c(1 / 6, 4 / 6, NA, 0.5)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$PGS.with.replaced.missing,
            c(1, 4, 2.5, 2.5)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$PGS,
            c(1, 4, 0, 2)
            );

        # check that percentiles are calculated correctly
        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$percentile,
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$percentile
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$percentile,
            test.missing.genotype.both.percentile.check[[PGS.OUTPUT.INDEX]]$percentile
            );


        # check missing genotype counts
        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 4, 2)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 4, 2)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 4, 2)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 4, 2)
            );

        # check missing genotype percentages
        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.25, 0.25, 1, 0.5)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.25, 0.25, 1, 0.5)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.25, 0.25, 1, 0.5)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.25, 0.25, 1, 0.5)
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
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$PGS.with.replaced.missing,
            c(1, 4, 3, 3.5)
            );
        }
    );

test_that(
    'apply.polygenic.score correctly validates phenotype data', {

        load('data/phenotype.test.data.Rda');

        # check that phenotype data is a data frame
        expect_error(
            apply.polygenic.score(
                vcf.data = phenotype.test.data$vcf.data,
                pgs.weight.data = phenotype.test.data$pgs.weight.data,
                phenotype.data = 'not a data frame'
                ),
            'phenotype.data must be a data.frame'
            );

        # check that phenotype data has correct columns
        phenotype.test.data.missing.columns <- phenotype.test.data$phenotype.data;
        phenotype.test.data.missing.columns$Indiv <- NULL;
        expect_error(
            apply.polygenic.score(
                vcf.data = phenotype.test.data$vcf.data,
                pgs.weight.data = phenotype.test.data$pgs.weight.data,
                phenotype.data = phenotype.test.data.missing.columns
                ),
            'phenotype.data must contain columns named Indiv'
            );

        # check for at least one matching sample in phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = phenotype.test.data$vcf.data,
                pgs.weight.data = phenotype.test.data$pgs.weight.data,
                phenotype.data = data.frame(Indiv = c('sample11'))
                ),
            'No matching Indiv between phenotype.data and vcf.data'
            );

        }
    );

test_that(
    'apply.polygenic.score correctly aggregates phenotype data', {
        load('data/phenotype.test.data.Rda');
        test.phenotype.output <- apply.polygenic.score(
            vcf.data = phenotype.test.data$vcf.data,
            pgs.weight.data = phenotype.test.data$pgs.weight.data,
            phenotype.data = phenotype.test.data$phenotype.data,
            missing.genotype.method = 'none'
            );

        # check that output is a data.frame
        expect_s3_class(
            test.phenotype.output[[PGS.OUTPUT.INDEX]],
            'data.frame'
            );

        # check for the correct column names
        expect_true(
            all(c(colnames(phenotype.test.data), 'PGS') %in% colnames(test.phenotype.output[[PGS.OUTPUT.INDEX]]))
            );

        # check that phenotype values are correctly matched to samples
        expect_equal(
            test.phenotype.output[[PGS.OUTPUT.INDEX]]$continuous.phenotype[match(phenotype.test.data$phenotype.data$Indiv, test.phenotype.output[[PGS.OUTPUT.INDEX]]$Indiv)],
            phenotype.test.data$phenotype.data$continuous.phenotype
            );

        expect_equal(
            test.phenotype.output[[PGS.OUTPUT.INDEX]]$binary.phenotype[match(phenotype.test.data$phenotype.data$Indiv, test.phenotype.output[[PGS.OUTPUT.INDEX]]$Indiv)],
            phenotype.test.data$phenotype.data$binary.phenotype
            );

        }
    );

test_that(
    'apply.polygenic.score works correctly on real data', {
        test.vcf.data <- import.vcf('data/HG001_GIAB.vcf.gz')
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');
        test.phenotype.data <- data.frame(Indiv = c('HG001', '2:HG001'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1));

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

        expect_no_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = test.phenotype.data,
                phenotype.analysis.columns = c('continuous.phenotype')
                )
            )

        }
    );
