PGS.OUTPUT.INDEX <- 1;
REGRESSION.OUTPUT.INDEX <- 2;

vcf.import.split.name <- 'split.wide.vcf.matrices';
vcf.import.long.name <- 'combined.long.vcf.df';

test_that(
    'apply.polygenic.score correctly checks inputs', {
        test.vcf.data <- import.vcf('data/HG001_GIAB.vcf.gz', long.format = TRUE);
        vcf.data.wide <- test.vcf.data[[vcf.import.split.name]];
        vcf.data.long <- test.vcf.data[[vcf.import.long.name]]$dat;

        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');
        test.phenotype.data <- data.frame(Indiv = c('HG001', '2:HG001'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1));

        # check that only data frame inputs are accepted
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'vcf.data must be a data.frame'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data,
                vcf.long.format = FALSE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'vcf.data must contain named elements: genotyped.alleles, vcf.fixed.fields'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data
                ),
            'pgs.weight.data must be a data.frame'
            );

        # check that missing genotype method input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = 'not a valid method'
                ),
            'missing.genotype.method must be either "mean.dosage", "normalize", or "none"'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('mean.dosage', 'normalize', 'not a valid method')
                ),
            'missing.genotype.method must be either "mean.dosage", "normalize", or "none"'
            );

        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('none', 'normalize', 'mean.dosage')
                )
            );

        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.long,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = c('none', 'normalize', 'mean.dosage')
                )
            );

        # check that analysis.source.pgs input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = c('mean.dosage', 'normalize')
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = 'not a valid method'
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                missing.genotype.method = 'normalize',
                analysis.source.pgs = 'mean.dosage'
                ),
            'analysis.source.pgs must be one of the chosen missing genotype methods'
            );
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                analysis.source.pgs = 'mean.dosage'
                )
            );

        # check that phenotype data input is correct
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = 'not a data frame'
                ),
            'phenotype.data must be a data.frame'
            );

        # check for matching samples between phenotype and vcf data
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = data.frame(Indiv = c('sample11'))
                ),
            'No matching Indiv between phenotype.data and vcf.data'
            );

        # check required columns in phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = subset(test.phenotype.data, select = -Indiv)
                ),
            'phenotype.data must contain columns named Indiv'
            );

        # check for correct phenotype analysis columns
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = test.phenotype.data,
                phenotype.analysis.columns = c('not a valid column')
                ),
            'phenotype.analysis.columns must be columns in phenotype.data'
            );

        # check for missing phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.analysis.columns = 'continuous.phenotype'
                ),
            'phenotype.analysis.columns provided but no phenotype data detected'
            );

        # check that required columns are present
        test.vcf.data.missing.columns <- vcf.data.long;
        test.vcf.data.missing.columns$gt_GT_alleles <- NULL;
        test.pgs.weight.data.missing.columns <- test.pgs.weight.data$pgs.weight.data;
        test.pgs.weight.data.missing.columns$beta <- NULL;
        test.pgs.weight.data.missing.columns$other_allele <- NULL;

        test.vcf.data.wide.missing.columns <- vcf.data.wide;
        test.vcf.data.wide.missing.columns$vcf.fixed.fields$CHROM <- NULL;

        # test vcf data for missing columns (long format)
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data.missing.columns,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                correct.strand.flips = FALSE
                ),
            'vcf.data must contain columns named CHROM, POS, REF, ALT, Indiv, gt_GT_alleles'
            );

        # test vcf data for missing columns (wide format)
        expect_error(
            apply.polygenic.score(
                vcf.data = test.vcf.data.wide.missing.columns,
                vcf.long.format = FALSE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                correct.strand.flips = FALSE
                ),
            'vcf.data\\$vcf.fixed.fields must contain columns named CHROM, POS, REF, ALT, allele.matrix.row.index'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.missing.columns,
                correct.strand.flips = FALSE
                ),
            'pgs.weight.data must contain columns named CHROM, POS, effect_allele, and beta'
            );

        # check for effect allele frequency column when use.external.effect.allele.frequency option is selected
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                use.external.effect.allele.frequency = TRUE,
                correct.strand.flips = FALSE
                ),
            'pgs.weight.data must contain a column named allelefrequency_effect if use.external.effect.allele.frequency is TRUE'
            );

        # check for other_allele column when remove.ambiguous.allele.matches, or remove.mismatched.indels are selected
        test.pgs.weight.data.missing.other.allele <- test.pgs.weight.data$pgs.weight.data;
        test.pgs.weight.data.missing.other.allele$other_allele <- NULL;
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.missing.other.allele,
                correct.strand.flips = TRUE
                ),
            'pgs.weight.data must contain a column named other_allele if correct.strand.flips, remove.ambiguous.allele.matches, or remove.mismatched.indels is TRUE'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.missing.other.allele,
                remove.ambiguous.allele.matches = TRUE,
                correct.strand.flips = FALSE
                ),
            'pgs.weight.data must contain a column named other_allele if correct.strand.flips, remove.ambiguous.allele.matches, or remove.mismatched.indels is TRUE'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.missing.other.allele,
                remove.mismatched.indels = TRUE,
                correct.strand.flips = FALSE
                ),
            'pgs.weight.data must contain a column named other_allele if correct.strand.flips, remove.ambiguous.allele.matches, or remove.mismatched.indels is TRUE'
            );

        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.missing.other.allele,
                remove.ambiguous.allele.matches = TRUE,
                remove.mismatched.indels = TRUE,
                correct.strand.flips = TRUE
                ),
            'pgs.weight.data must contain a column named other_allele if correct.strand.flips, remove.ambiguous.allele.matches, or remove.mismatched.indels is TRUE'
            );


        # check for duplicate coordinates in PGS data
        duplicate.row <- test.pgs.weight.data$pgs.weight.data[1, ];
        duplicate.row.as.multiallelic <- duplicate.row;
        duplicate.row.as.multiallelic$effect_allele <- 'C';
        test.pgs.weight.data.duplicated.coordinates <- rbind(test.pgs.weight.data$pgs.weight.data, duplicate.row.as.multiallelic);
        expect_warning(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.duplicated.coordinates
                ),
            'Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.'
            );

        # check for duplicate variants in PGS data
        test.pgs.weight.data.duplicated.variants <- rbind(test.pgs.weight.data$pgs.weight.data, duplicate.row);
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                pgs.weight.data = test.pgs.weight.data.duplicated.variants
                ),
            'Duplicate variants detected in the PGS weight data. Please ensure only unique coordinate:effect allele combinations are present.'
            );

        # Verify correct number of variants and samples
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.long[1:3, ],
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'Number of vcf data rows is not equivalent to number of samples times number of variants. Please ensure that all samples have variant data represented for all variants.'
            );

        vcf.data.wide.incorrect.matrix <- vcf.data.wide;
        vcf.data.wide.incorrect.matrix$genotyped.alleles <- vcf.data.wide.incorrect.matrix$genotyped.alleles[1:3, ];
        expect_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide.incorrect.matrix,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                ),
            'vcf.data\\$genotyped.alleles and vcf.data\\$vcf.fixed.fields must have the same number of rows'
            );

        }
    );

test_that(
    'apply.polygenic.score correctly outputs validation only option', {
        load('data/simple.pgs.application.test.data.Rda');
        simple.pgs.application.test.data$wide.vcf.data <- convert.long.vcf.to.wide.vcf(simple.pgs.application.test.data$vcf.data);
        expect_equal(
            apply.polygenic.score(
                vcf.data = simple.pgs.application.test.data$vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
                validate.inputs.only = TRUE
                ),
            TRUE
            );
        expect_message(
            apply.polygenic.score(
                vcf.data = simple.pgs.application.test.data$vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
                validate.inputs.only = TRUE
                ),
            'Input data passed validation'
            );

        expect_equal(
            apply.polygenic.score(
                vcf.data = simple.pgs.application.test.data$wide.vcf.data,
                vcf.long.format = FALSE,
                pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
                validate.inputs.only = TRUE
                ),
            TRUE
            );
        expect_message(
            apply.polygenic.score(
                vcf.data = simple.pgs.application.test.data$wide.vcf.data,
                vcf.long.format = FALSE,
                pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
                validate.inputs.only = TRUE
                ),
            'Input data passed validation'
            );
        }
    );

test_that(
    'apply.polygenic.score correctly formats general output', {
        load('data/simple.pgs.application.test.data.Rda');
        simple.pgs.application.test.data$wide.vcf.data <- convert.long.vcf.to.wide.vcf(simple.pgs.application.test.data$vcf.data);
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            vcf.long.format = TRUE,
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

        # wide format
        test.pgs.per.sample.wide <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );

        # check that wide format output is idenical to long format output
        expect_equal(
            test.pgs.per.sample,
            test.pgs.per.sample.wide
            );


        }
    );

test_that(
    'apply.polygenic.score correctly formats pgs output', {
        load('data/simple.pgs.application.test.data.Rda')
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            vcf.long.format = TRUE,
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
            9
            );
        }
    );

test_that(
    'apply.polygenic.score correctly formats regression output', {
        load('data/simple.pgs.application.test.data.Rda')
        simple.pgs.application.test.data$wide.vcf.data <- convert.long.vcf.to.wide.vcf(simple.pgs.application.test.data$vcf.data);
        test.pgs.per.sample.with.phenotype <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            vcf.long.format = TRUE,
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

        # wide format equivalence check
        test.pgs.per.sample.with.phenotype.wide <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data,
            phenotype.data = data.frame(Indiv = c('sample1', 'sample2'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1)),
            phenotype.analysis.columns = c('continuous.phenotype', 'binary.phenotype')
            );

        expect_equal(
            test.pgs.per.sample.with.phenotype,
            test.pgs.per.sample.with.phenotype.wide
            );

        # check NULL output when no phenotype data is provided
        test.pgs.per.sample.no.phenotype <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            vcf.long.format = TRUE,
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
        simple.pgs.application.test.data$wide.vcf.data <- convert.long.vcf.to.wide.vcf(simple.pgs.application.test.data$vcf.data);
        test.pgs.per.sample <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$vcf.data,
            vcf.long.format = TRUE,
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

        # wide format equivalence check
        test.pgs.per.sample.wide.format <- apply.polygenic.score(
            vcf.data = simple.pgs.application.test.data$wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = simple.pgs.application.test.data$pgs.weight.data
            );
        expect_equal(
            test.pgs.per.sample,
            test.pgs.per.sample.wide.format
            );
        }
    );

test_that(
    'apply.polygenic.score correctly handles multiallic sites', {
        load('data/merged.multiallelic.site.test.data.Rda');
        merged.multiallelic.site.test.data$merged.multiallelic.wide.vcf.data <- convert.long.vcf.to.wide.vcf(merged.multiallelic.site.test.data$merged.multiallelic.vcf.data);
        # test case with VCF multiallelic site with no extra beta, REF allele is the risk allele
        ref.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.as.single.risk.allele.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            ref.as.single.risk.allele.test[[PGS.OUTPUT.INDEX]]$PGS,
            c(4, 3, 0, 0)
            );

        ref.as.single.risk.allele.wide.format.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.as.single.risk.allele.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            ref.as.single.risk.allele.test,
            ref.as.single.risk.allele.wide.format.test
            );

        # test case with VCF multiallelic sites with no extra beta, ALT allele is the risk allele
        alt.as.single.risk.allele.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            alt.as.single.risk.allele.test[[PGS.OUTPUT.INDEX]]$PGS,
            c(2, 1, 2, 0)
            );
        expect_no_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data
                )
            );

        # wide format equivalence check
        alt.as.single.risk.allele.wide.format.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.single.risk.allele.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            alt.as.single.risk.allele.test,
            alt.as.single.risk.allele.wide.format.test
            );

        # test case with a VCF multiallelic site that also has an extra beta, both ALT alleles are risk alleles
        alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            alt.as.two.risk.alleles.test[[PGS.OUTPUT.INDEX]]$PGS,
            c(2, 1.5, 3, 0)
            );
        expect_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data
                ),
            'Duplicate variants detected in the PGS weight data. These will be treated as multiallelic sites.'
            );

        # wide format equivalence check
        alt.as.two.risk.alleles.wide.format.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = merged.multiallelic.site.test.data$alt.as.two.risk.alleles.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            alt.as.two.risk.alleles.test,
            alt.as.two.risk.alleles.wide.format.test
            );

        # test case with a VCF multiallelic site that also has an extra beta, one REF and one ALT allele are risk alleles
        ref.and.alt.as.two.risk.alleles.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            ref.and.alt.as.two.risk.alleles.test[[PGS.OUTPUT.INDEX]]$PGS,
            c(2, 1.5, 2, 0)
            );
        expect_warning(
            apply.polygenic.score(
                vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.vcf.data,
                vcf.long.format = TRUE,
                pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data
                ),
            'Multiple effect alleles found in a sample, choosing effect allele with highest beta for dosage calculation. Check coordinates: chr2:2'
            );

        # wide format equivalence check
        ref.and.alt.as.two.risk.alleles.wide.format.test <- apply.polygenic.score(
            vcf.data = merged.multiallelic.site.test.data$merged.multiallelic.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = merged.multiallelic.site.test.data$ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            ref.and.alt.as.two.risk.alleles.test,
            ref.and.alt.as.two.risk.alleles.wide.format.test
            );
        }
    );

test_that(
    'apply.polygenic.score correctly handles missing genotypes', {
        load('data/missing.genotype.test.data.Rda');
        missing.genotype.test.data$missing.genotype.wide.vcf.data <- convert.long.vcf.to.wide.vcf(missing.genotype.test.data$missing.genotype.vcf.data);
        test.missing.genotype.mean.dosage <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage'
            );
        test.missing.genotype.normalize <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'normalize'
            );
        test.missing.genotype.both <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'mean.dosage'
            );
        test.missing.genotype.both.percentile.check <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'normalize'
            );
        test.missing.genotype.none <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'none'
            );

        # check column names
        percentile.colnames <- c('percentile', 'decile', 'quartile');
        extra.colnames <- c('n.pgm.sites', 'n.missing.genotypes', 'percent.missing.genotypes', 'n.non.missing.alleles');
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
            c(1, 5, 19 / 6, 3.5)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$PGS.with.normalized.missing,
            c(1 / 8, 5 / 8, NA, 3 / 6)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$PGS.with.normalized.missing,
            c(1 / 8, 5 / 8, NA, 3 / 6)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$PGS.with.replaced.missing,
            c(1, 5, 19 / 6, 3.5)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$PGS,
            c(1, 5, 0, 3)
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
            c(1, 1, 5, 2)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 5, 2)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 5, 2)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$n.missing.genotypes,
            c(1, 1, 5, 2)
            );

        # check missing genotype percentages
        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.20, 0.20, 1, 0.40)
            );
        expect_equal(
            test.missing.genotype.normalize[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.20, 0.20, 1, 0.40)
            );
        expect_equal(
            test.missing.genotype.both[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.20, 0.20, 1, 0.40)
            );
        expect_equal(
            test.missing.genotype.none[[PGS.OUTPUT.INDEX]]$percent.missing.genotypes,
            c(0.20, 0.20, 1, 0.40)
            );

        # wide format equivalence checks
        test.missing.genotype.mean.dosage.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage'
            );
        expect_equal(
            test.missing.genotype.mean.dosage,
            test.missing.genotype.mean.dosage.wide
            );

        test.missing.genotype.normalize.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'normalize'
            );
        expect_equal(
            test.missing.genotype.normalize,
            test.missing.genotype.normalize.wide
            );

        test.missing.genotype.both.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'mean.dosage'
            );
        expect_equal(
            test.missing.genotype.both,
            test.missing.genotype.both.wide
            );

        test.missing.genotype.both.percentile.check.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = c('mean.dosage', 'normalize'),
            analysis.source.pgs = 'normalize'
            );
        expect_equal(
            test.missing.genotype.both.percentile.check,
            test.missing.genotype.both.percentile.check.wide
            );

        test.missing.genotype.none.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            test.missing.genotype.none,
            test.missing.genotype.none.wide
            );
    }
    );


test_that(
    'apply.polygenic.score correctly handles external effect allele frequency', {
        load('data/missing.genotype.test.data.Rda');
        missing.genotype.test.data$missing.genotype.wide.vcf.data <- convert.long.vcf.to.wide.vcf(missing.genotype.test.data$missing.genotype.vcf.data);
        # add effect allele frequency column to PGS weight data
        missing.genotype.test.data$missing.genotype.pgs.weight.data$allelefrequency_effect <- c(0.5, 0.25, 0.75, 0.5, 0.25);
        test.missing.genotype.mean.dosage <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage',
            use.external.effect.allele.frequency = TRUE
            );

        expect_equal(
            test.missing.genotype.mean.dosage[[PGS.OUTPUT.INDEX]]$PGS.with.replaced.missing,
            #c(1, 5, 3.5, 4.5)
            # c(1.5, 5.5, 4.5, 5)
            c(2, 6, 4.5, 5.5)
            );

        # wide format equivalence check
        test.missing.genotype.mean.dosage.wide <- apply.polygenic.score(
            vcf.data = missing.genotype.test.data$missing.genotype.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = missing.genotype.test.data$missing.genotype.pgs.weight.data,
            missing.genotype.method = 'mean.dosage',
            use.external.effect.allele.frequency = TRUE
            );
        expect_equal(
            test.missing.genotype.mean.dosage,
            test.missing.genotype.mean.dosage.wide
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
                vcf.long.format = TRUE,
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
                vcf.long.format = TRUE,
                pgs.weight.data = phenotype.test.data$pgs.weight.data,
                phenotype.data = phenotype.test.data.missing.columns
                ),
            'phenotype.data must contain columns named Indiv'
            );

        # check for at least one matching sample in phenotype data
        expect_error(
            apply.polygenic.score(
                vcf.data = phenotype.test.data$vcf.data,
                vcf.long.format = TRUE,
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
        phenotype.test.data$vcf.wide.data <- convert.long.vcf.to.wide.vcf(phenotype.test.data$vcf.data);
        test.phenotype.output <- apply.polygenic.score(
            vcf.data = phenotype.test.data$vcf.data,
            vcf.long.format = TRUE,
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

        # wide format equivalence check
        test.phenotype.output.wide <- apply.polygenic.score(
            vcf.data = phenotype.test.data$vcf.wide.data,
            vcf.long.format = FALSE,
            pgs.weight.data = phenotype.test.data$pgs.weight.data,
            phenotype.data = phenotype.test.data$phenotype.data,
            missing.genotype.method = 'none'
            );
        expect_equal(
            test.phenotype.output,
            test.phenotype.output.wide
            );

        }
    );

test_that(
    'apply.polygenic.score works correctly on real data', {
        test.vcf.data <- import.vcf('data/HG001_GIAB.vcf.gz', long.format = TRUE)
        vcf.data.wide <- test.vcf.data[[vcf.import.split.name]];
        vcf.data.long <- test.vcf.data[[vcf.import.long.name]]$dat;
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');
        test.phenotype.data <- data.frame(Indiv = c('HG001', '2:HG001'), continuous.phenotype = c(1, 2), binary.phenotype = c(0, 1));

        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.long,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            )
        # wide format
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                vcf.long.format = FALSE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            )

        # custom percentiles
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.long,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                n.percentiles = 5
                )
            )

        # wide format
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                vcf.long.format = FALSE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                n.percentiles = 5
                )
            )

        # with phenotype data
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.long,
                vcf.long.format = TRUE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = test.phenotype.data,
                phenotype.analysis.columns = c('continuous.phenotype')
                )
            )

        # wide format
        expect_no_error(
            apply.polygenic.score(
                vcf.data = vcf.data.wide,
                vcf.long.format = FALSE,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
                phenotype.data = test.phenotype.data,
                phenotype.analysis.columns = c('continuous.phenotype')
                )
            );

        # check file writing
        temp.dir <- tempdir();

        apply.polygenic.score(
            vcf.data = vcf.data.long,
            vcf.long.format = TRUE,
            pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
            phenotype.data = test.phenotype.data,
            phenotype.analysis.columns = c('continuous.phenotype'),
            output.dir = temp.dir,
            file.prefix = 'TEST-apply-pgs'
            );

        test.pgs.filename <- generate.filename(
            project.stem = 'TEST-apply-pgs',
            file.core = 'per-sample-pgs-summary',
            extension = 'txt'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.pgs.filename))
            );

        test.regression.filename <- generate.filename(
            project.stem = 'TEST-apply-pgs',
            file.core = 'pgs-regression-output',
            extension = 'txt'
            );
        expect_true(
            file.exists(file.path(temp.dir, test.regression.filename))
            );

        # wide format
        apply.polygenic.score(
            vcf.data = vcf.data.wide,
            vcf.long.format = FALSE,
            pgs.weight.data = test.pgs.weight.data$pgs.weight.data,
            phenotype.data = test.phenotype.data,
            phenotype.analysis.columns = c('continuous.phenotype'),
            output.dir = temp.dir,
            file.prefix = 'TEST-apply-pgs-wide'
            );

        test.pgs.wide.filename <- generate.filename(
            project.stem = 'TEST-apply-pgs-wide',
            file.core = 'per-sample-pgs-summary',
            extension = 'txt'
            );

        expect_true(
            file.exists(file.path(temp.dir, test.pgs.wide.filename))
            );

        test.regression.wide.filename <- generate.filename(
            project.stem = 'TEST-apply-pgs-wide',
            file.core = 'pgs-regression-output',
            extension = 'txt'
            );

        expect_true(
            file.exists(file.path(temp.dir, test.regression.wide.filename))
            );

        }
    );

test_that(
    'apply.polygenic.score correctly handles strand flipping', {
        load('data/strand.flip.test.data.Rda');
        strand.flip.test.data$strand.flip.wide.vcf.data <- convert.long.vcf.to.wide.vcf(strand.flip.test.data$strand.flip.vcf.data);
        # no strand flip correction
        strand.flip.raw.data <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            correct.strand.flips = FALSE
            );

        expect_equal(
            strand.flip.raw.data$pgs.output$PGS,
            c(4, 4, 4)
            );

        expect_equal(
            strand.flip.raw.data$pgs.output$n.missing.genotypes,
            c(0, 0, 0)
            );

        # wide format equivalence check
        strand.flip.raw.wide.data <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            correct.strand.flips = FALSE
            );
        expect_equal(
            strand.flip.raw.data,
            strand.flip.raw.wide.data
            );

        # strand flip correction
        strand.flip.correct.flips <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            correct.strand.flips = TRUE
            );

        expect_equal(
            strand.flip.correct.flips$pgs.output$PGS,
            c(6, 6, 6)
            );

        expect_equal(
            strand.flip.correct.flips$pgs.output$n.missing.genotypes,
            c(0, 0, 0)
            );

        # wide format equivalence check
        strand.flip.correct.flips.wide <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            correct.strand.flips = TRUE
            );
        expect_equal(
            strand.flip.correct.flips,
            strand.flip.correct.flips.wide
            );

        # ambiguous allele matches removed
        strand.flip.remove.ambiguous <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.vcf.data,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            vcf.long.format = TRUE,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE
            );

        expect_equal(
            strand.flip.remove.ambiguous$pgs.output$PGS,
            c(4, 5, 6)
            );

        expect_equal(
            strand.flip.remove.ambiguous$pgs.output$n.missing.genotypes,
            c(2, 2, 2)
            );

        # wide format equivalence check
        strand.flip.remove.ambiguous.wide <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.wide.vcf.data,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            vcf.long.format = FALSE,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE
            );
        expect_equal(
            strand.flip.remove.ambiguous,
            strand.flip.remove.ambiguous.wide
            );

        # ambiguous allele matches not removed due to threshold
        strand.flip.threshold.ambiguous <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE,
            max.strand.flips = 3 # there are two unambiguous strand flips in the test data, setting threshold above that
            );

        # expect equal to not removing ambiguous allele matches
        expect_equal(
            strand.flip.threshold.ambiguous$pgs.output$PGS, # dosage is unaffected because unresolved mismatch does not contain an effect allele so dosage is 0
            strand.flip.correct.flips$pgs.output$PGS
            );

        expect_equal(
            strand.flip.threshold.ambiguous$pgs.output$n.missing.genotypes,
            c(1, 1, 1) # the unresolved mismatch is marked missing
            );

        # wide format equivalence check
        strand.flip.threshold.ambiguous.wide <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE,
            max.strand.flips = 3
            );
        expect_equal(
            strand.flip.threshold.ambiguous,
            strand.flip.threshold.ambiguous.wide
            );

        # indels removed
        strand.flip.remove.indel.mismatches <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.vcf.data,
            vcf.long.format = TRUE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE,
            remove.mismatched.indels = TRUE
            );

        expect_equal(
            strand.flip.remove.indel.mismatches$pgs.output$PGS,
            c(4, 5, 6)
            );

        expect_equal(
            strand.flip.remove.indel.mismatches$pgs.output$n.missing.genotypes,
            c(3, 3, 3)
            );
        # wide format equivalence check
        strand.flip.remove.indel.mismatches.wide <- apply.polygenic.score(
            vcf.data = strand.flip.test.data$strand.flip.wide.vcf.data,
            vcf.long.format = FALSE,
            pgs.weight.data = strand.flip.test.data$strand.flip.pgs.weight.data,
            correct.strand.flips = TRUE,
            missing.genotype.method = c('none', 'mean.dosage', 'normalize'),
            remove.ambiguous.allele.matches = TRUE,
            remove.mismatched.indels = TRUE
            );
        expect_equal(
            strand.flip.remove.indel.mismatches,
            strand.flip.remove.indel.mismatches.wide
            );
        }
    );
