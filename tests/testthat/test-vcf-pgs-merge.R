test_that(
    'merge.vcf.with.pgs correctly checks inputs', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that inputs are data.frames
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = 'not a data.frame',
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must be a data.frame'
            );
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = 'not a data.frame'
                ),
            'pgs.weight.data must be a data.frame'
            );

        # check required columns
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = data.frame(
                    CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
                    foo = c(1, 2, 3, 4),
                    REF = c('A', 'T', 'C', 'G'),
                    ALT = c('T', 'A', 'G', 'C')
                    ),
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must contain columns named CHROM and POS'
            );
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = data.frame(
                    foo = c('chr1', 'chr1', 'chr1', 'chr1'),
                    POS = c(1, 2, 3, 4),
                    REF = c('A', 'T', 'C', 'G'),
                    ALT = c('T', 'A', 'G', 'C')
                    ),
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must contain columns named CHROM and POS'
            );
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = data.frame(
                    CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
                    foo = c(1, 2, 3, 4),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
                ),
            'pgs.weight.data must contain columns named CHROM and POS'
            );
        expect_error(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = data.frame(
                    foo = c('chr1', 'chr1', 'chr1', 'chr1'),
                    POS = c(1, 2, 3, 4),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
                ),
            'pgs.weight.data must contain columns named CHROM and POS'
            );
        }
    )

test_that(
    'merge.vcf.with.pgs correctly formats CHROM column', {

        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        chr.prefix.vcf.no.prefix.pgs <- merge.vcf.with.pgs(
            vcf.data = test.vcf.data,
            pgs.weight.data = data.frame(
                CHROM = c('1', '1', '1', '1'),
                POS = c(1, 2, 3, 4),
                effect_allele = c('A', 'T', 'C', 'G'),
                beta = c(1, 2, 3, 4)
                )
            )$merged.vcf.with.pgs.data;
        no.prefix.vcf.chr.prefix.pgs <- merge.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('1', '1', '1', '1'),
                POS = c(1, 2, 3, 4),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = test.pgs.weight.data
            )$merged.vcf.with.pgs.data;
        numeric.sex.vcf.alpha.sex.pgs <- merge.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('23', '23', '24', '24'),
                POS = c(1, 2, 3, 4),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = data.frame(
                    CHROM = c('X', 'X', 'Y', 'Y'),
                    POS = c(1, 2, 3, 4),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
            )$merged.vcf.with.pgs.data;
        alpha.sex.vcf.numeric.sex.pgs <- merge.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('X', 'X', 'Y', 'Y'),
                POS = c(1, 2, 3, 4),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = data.frame(
                    CHROM = c('23', '23', '24', '24'),
                    POS = c(1, 2, 3, 4),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
            )$merged.vcf.with.pgs.data;

        expect_equal(
            chr.prefix.vcf.no.prefix.pgs$CHROM,
            c('chr1', 'chr1', 'chr1', 'chr1')
            );
        expect_equal(
            no.prefix.vcf.chr.prefix.pgs$CHROM,
            c('1', '1', '1', '1')
            );
        expect_equal(
            numeric.sex.vcf.alpha.sex.pgs$CHROM,
            c('23', '23', '24', '24')
            );
        expect_equal(
            alpha.sex.vcf.numeric.sex.pgs$CHROM,
            c('X', 'X', 'Y', 'Y')
            );
        }
    );

test_that(
    'merge.vcf.with.pgs correctly formats output', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 5),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that no error is thrown when inputs are correct
        expect_no_error(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = test.pgs.weight.data
                )
            );

        test.merge.vcf.with.pgs <- merge.vcf.with.pgs(
            vcf.data = test.vcf.data,
            pgs.weight.data = test.pgs.weight.data
            );

        # check that merge.vcf.with.pgs returns a list
        expect_type(
            test.merge.vcf.with.pgs,
            'list'
            );

        # check that merge.vcf.with.pgs returns a data.frame
        expect_s3_class(
            test.merge.vcf.with.pgs$merged.vcf.with.pgs.data,
            'data.frame'
            );

        # check that merge.vcf.with.pgs returns a data.frame
        expect_s3_class(
            test.merge.vcf.with.pgs$missing.snp.data,
            'data.frame'
            );

        }
    );

test_that(
    'merge.vcf.with.pgs correctly performs merge', {
        test.vcf.data.no.missing <- data.frame(
            CHROM = c('chr1', 'chr2', 'chr3', 'chr4'),
            POS = c(1, 2, 3, 4),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.vcf.data.missing <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 2, 5),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 2, 4),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that merge.vcf.with.pgs throws a warning for missing SNPs
        expect_warning(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data.missing,
                pgs.weight.data = test.pgs.weight.data
                ),
            'PGS is missing 1 SNPs from VCF'
            );

        test.merge.vcf.with.pgs.no.missing <- merge.vcf.with.pgs(
            vcf.data = test.vcf.data.no.missing,
            pgs.weight.data = test.pgs.weight.data
            );

        test.merge.vcf.with.pgs.missing <- merge.vcf.with.pgs(
            vcf.data = test.vcf.data.missing,
            pgs.weight.data = test.pgs.weight.data
            );

        # check that merge.vcf.with.pgs returns the correct number of rows
        expect_equal(
            nrow(test.merge.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.merge.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.merge.vcf.with.pgs.missing$missing.snp.data),
            1
            );
        expect_equal(
            test.merge.vcf.with.pgs.no.missing$missing.snp.data,
            NULL
            );

        # check that merge.vcf.with.pgs returns the correct number of columns
        expect_equal(
            ncol(test.merge.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data),
            8
            );
        expect_equal(
            ncol(test.merge.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            8
            );

        # check that merge.vcf.with.pgs returns the correct columns
        expect_true(
            all(colnames(test.merge.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data) %in% c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta'))
            );
        expect_true(
            all(colnames(test.merge.vcf.with.pgs.missing$merged.vcf.with.pgs.data) %in% c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta'))
            );

        # check that merge.vcf.with.pgs returns the correct values
        expect_equal(
            test.merge.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr2', 'chr3', 'chr4')
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr2', 'chr3', 'chr4')
            );
        expect_equal(
            test.merge.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$POS,
            c(1, 2, 3, 4)
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$merged.vcf.with.pgs.data$POS,
            c(1, 2, 3, 4)
            );
        # check that merge.vcf.with.pgs returns the correct values for missing SNPs
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$CHROM,
            'chr4'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$POS,
            4
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$REF,
            'NA'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$ALT,
            'NA'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$Indiv,
            'NA'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$gt_GT_alleles,
            'NA'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$effect_allele,
            'G'
            );
        expect_equal(
            test.merge.vcf.with.pgs.missing$missing.snp.data$beta,
            4
            );
        }
    );

test_that(
    'merge.vcf.with.pgs works correctly with real data', {
        test.vcf.data <- import.vcf('data/HG001_GRCh38_1_22_v4.2.1_benchmark_in_PGS003378_hmPOS_GRCh38_slop10_duplicated-sample.vcf.gz')
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');

        expect_no_error(
            merge.vcf.with.pgs(
                vcf.data = test.vcf.data$dat,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            );
        }
    )
