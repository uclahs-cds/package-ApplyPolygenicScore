test_that(
    'combine.vcf.with.pgs correctly checks inputs', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID  = c('rs1', 'rs2', 'rs3', 'rs4'),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that inputs are data.frames
        expect_error(
            combine.vcf.with.pgs(
                vcf.data = 'not a data.frame',
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must be a data.frame'
            );
        expect_error(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = 'not a data.frame'
                ),
            'pgs.weight.data must be a data.frame'
            );

        # check required columns
        # missing POS
        expect_error(
            combine.vcf.with.pgs(
                vcf.data = data.frame(
                    CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
                    foo = c(1, 2, 3, 4),
                    ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                    REF = c('A', 'T', 'C', 'G'),
                    ALT = c('T', 'A', 'G', 'C')
                    ),
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must contain columns named CHROM and POS'
            );
        # missing CHROM
        expect_error(
            combine.vcf.with.pgs(
                vcf.data = data.frame(
                    foo = c('chr1', 'chr1', 'chr1', 'chr1'),
                    POS = c(1, 2, 3, 4),
                    ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                    REF = c('A', 'T', 'C', 'G'),
                    ALT = c('T', 'A', 'G', 'C')
                    ),
                pgs.weight.data = test.pgs.weight.data
                ),
            'vcf.data must contain columns named CHROM and POS'
            );
        # missing ID
        expect_warning(
            combine.vcf.with.pgs(
                vcf.data = data.frame(
                    CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
                    POS = c(1, 2, 3, 4),
                    foo = c('rs1', 'rs2', 'rs3', 'rs4'),
                    REF = c('A', 'T', 'C', 'G'),
                    ALT = c('T', 'A', 'G', 'C')
                    ),
                pgs.weight.data = test.pgs.weight.data
                ),
            'ID column not found in VCF data. Merging by rsID will not be possible.'
            );

        expect_warning(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = data.frame(
                    CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
                    POS = c(1, 2, 3, 4),
                    foo = c(1, 2, 3, 4),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
                ),
            'ID column not found in PGS weight data. Merging by rsID will not be possible.'
            );

        }
    )

test_that(
    'combine.vcf.with.pgs stops when all PGS SNPs are missing from the VCF', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(5, 6, 7, 8),
            ID  = c('rs5', 'rs6', 'rs7', 'rs8'),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        expect_error(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = test.pgs.weight.data
                ),
            'All PGS SNPs are missing from the VCF, terminating merge.'
            );
        }
    );

test_that(
    'combine.vcf.with.pgs correctly formats CHROM column', {

        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID  = c('rs1', 'rs2', 'rs3', 'rs4'),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        chr.prefix.vcf.no.prefix.pgs <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data,
            pgs.weight.data = data.frame(
                CHROM = c('1', '1', '1', '1'),
                POS = c(1, 2, 3, 4),
                ID  = c('rs1', 'rs2', 'rs3', 'rs4'),
                effect_allele = c('A', 'T', 'C', 'G'),
                beta = c(1, 2, 3, 4)
                )
            )$merged.vcf.with.pgs.data;
        no.prefix.vcf.chr.prefix.pgs <- combine.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('1', '1', '1', '1'),
                POS = c(1, 2, 3, 4),
                ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = test.pgs.weight.data
            )$merged.vcf.with.pgs.data;
        numeric.sex.vcf.alpha.sex.pgs <- combine.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('23', '23', '24', '24'),
                POS = c(1, 2, 3, 4),
                ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = data.frame(
                    CHROM = c('X', 'X', 'Y', 'Y'),
                    POS = c(1, 2, 3, 4),
                    ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                    effect_allele = c('A', 'T', 'C', 'G'),
                    beta = c(1, 2, 3, 4)
                    )
            )$merged.vcf.with.pgs.data;
        alpha.sex.vcf.numeric.sex.pgs <- combine.vcf.with.pgs(
            vcf.data = data.frame(
                CHROM = c('X', 'X', 'Y', 'Y'),
                POS = c(1, 2, 3, 4),
                ID = c('rs1', 'rs2', 'rs3', 'rs4'),
                REF = c('A', 'T', 'C', 'G'),
                ALT = c('T', 'A', 'G', 'C')
                ),
                pgs.weight.data = data.frame(
                    CHROM = c('23', '23', '24', '24'),
                    POS = c(1, 2, 3, 4),
                    ID = c('rs1', 'rs2', 'rs3', 'rs4'),
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
    'combine.vcf.with.pgs correctly formats output', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr1', 'chr1', 'chr1'),
            POS = c(1, 2, 3, 5),
            ID  = c('rs1', 'rs2', 'rs3', 'rs5'),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that no error is thrown when inputs are correct
        expect_no_error(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data,
                pgs.weight.data = test.pgs.weight.data
                )
            );

        test.combine.vcf.with.pgs <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data,
            pgs.weight.data = test.pgs.weight.data
            );

        # check that combine.vcf.with.pgs returns a list
        expect_type(
            test.combine.vcf.with.pgs,
            'list'
            );

        # check that combine.vcf.with.pgs returns a data.frame
        expect_s3_class(
            test.combine.vcf.with.pgs$merged.vcf.with.pgs.data,
            'data.frame'
            );

        # check that combine.vcf.with.pgs returns a data.frame
        expect_s3_class(
            test.combine.vcf.with.pgs$missing.snp.data,
            'data.frame'
            );

        }
    );

test_that(
    'combine.vcf.with.pgs correctly performs merge', {
        test.vcf.data.no.missing <- data.frame(
            CHROM = c('chr1', 'chr2', 'chr3', 'chr4'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.vcf.data.missing <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 2, 5),
            ID = c('rs1', 'rs3', 'rs2', 'rs5'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.vcf.data.missing.locus.matching.rsid <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 3, 6), #rs2/chr2:3 is missing by POS but present by rsID, #rs5/chr4:6 is missing by rsID AND by POS
            ID = c('rs1', 'rs3', 'rs2', 'rs5'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.vcf.data.missing.locus.matching.rsid.with.semicolons <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 3, 6),
            ID = c('rs1;rs1', 'rs3', 'rs2;rsB', 'rs5;rsC'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.pgs.weight.data <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 2, 4),
            ID  = c('rs1', 'rs3', 'rs2', 'rs4'),
            effect_allele = c('A', 'T', 'C', 'G'),
            beta = c(1, 2, 3, 4)
            );

        # check that combine.vcf.with.pgs throws a warning for missing SNPs
        expect_warning(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data.missing,
                pgs.weight.data = test.pgs.weight.data
                ),
            'PGS is missing 1 SNPs from VCF'
            );

        expect_warning(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data.missing.locus.matching.rsid,
                pgs.weight.data = test.pgs.weight.data
                ),
            'PGS is missing 1 SNPs from VCF'
            );

        expect_warning(
            combine.vcf.with.pgs(
                vcf.data = test.vcf.data.missing.locus.matching.rsid.with.semicolons,
                pgs.weight.data = test.pgs.weight.data
                ),
            'PGS is missing 1 SNPs from VCF'
            );

        test.combine.vcf.with.pgs.no.missing <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data.no.missing,
            pgs.weight.data = test.pgs.weight.data
            );

        test.combine.vcf.with.pgs.missing <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data.missing,
            pgs.weight.data = test.pgs.weight.data
            );

        test.combine.vcf.with.pgs.missing.locus.matching.rsid <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data.missing.locus.matching.rsid,
            pgs.weight.data = test.pgs.weight.data
            );

        test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data.missing.locus.matching.rsid.with.semicolons,
            pgs.weight.data = test.pgs.weight.data
            );

        # check that combine.vcf.with.pgs returns the correct number of rows
        expect_equal(
            nrow(test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data),
            4
            );

        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing$missing.snp.data),
            1
            );
        expect_equal(
            test.combine.vcf.with.pgs.no.missing$missing.snp.data,
            NULL
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data),
            1
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data),
            1
            );

        # check that combine.vcf.with.pgs returns the correct number of columns
        expect_equal(
            ncol(test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data),
            11
            );
        expect_equal(
            ncol(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            11
            );
        expect_equal(
            ncol(test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data),
            11
            );
        expect_equal(
            ncol(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data),
            11
            );

        # check that combine.vcf.with.pgs returns the correct columns
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'ID.pgs', 'ID.vcf', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data))
            );
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'ID.pgs', 'ID.vcf', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data))
            );
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'ID.pgs', 'ID.vcf', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data))
            );
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'ID.pgs', 'ID.vcf', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data))
            );

        # check that combine.vcf.with.pgs returns the correct values
        expect_equal(
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr2', 'chr3', 'chr4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr2', 'chr3', 'chr4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr3', 'chr2', 'chr4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data$CHROM,
            c('chr1', 'chr3', 'chr2', 'chr4')
            );

        expect_equal(
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$POS,
            c(1, 2, 3, 4)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data$POS,
            c(1, 2, 3, 4)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data$POS,
            c(1, 3, 3, 4)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data$POS,
            c(1, 3, 3, 4)
            );

        expect_equal(
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$ID.pgs,
            c('rs1', 'rs2', 'rs3', 'rs4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data$ID.pgs,
            c('rs1', 'rs2', 'rs3', 'rs4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data$ID.pgs,
            c('rs1', 'rs3', 'rs2', 'rs4')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data$ID.pgs,
            c('rs1', 'rs3', 'rs2', 'rs4')
            );

        expect_equal(
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$ID.pgs,
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$ID.vcf
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data$ID.vcf,
            c('rs1', 'rs2', 'rs3', NA)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data$ID.vcf,
            c('rs1', 'rs3', 'rs2', NA)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data$ID.vcf,
            c('rs1;rs1', 'rs3', 'rs2;rsB', NA)
            );

        expect_equal(
            test.combine.vcf.with.pgs.no.missing$merged.vcf.with.pgs.data$merge.strategy,
            rep('genomic coordinate', 4)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data$merge.strategy,
            c(rep('genomic coordinate', 3), 'rsID')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$merged.vcf.with.pgs.data$merge.strategy,
            c('genomic coordinate', 'genomic coordinate', 'rsID','rsID')
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$merged.vcf.with.pgs.data$merge.strategy,
            c('genomic coordinate', 'genomic coordinate', 'rsID','rsID')
            );

        # check that combine.vcf.with.pgs returns the correct values for missing SNPs
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$CHROM.pgs,
            'chr4'
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$CHROM.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$POS.pgs,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$POS.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$ID.pgs,
            'rs4'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$effect_allele,
            'G'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$beta,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$REF)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$ALT)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$Indiv)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing$missing.snp.data$gt_GT_alleles)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$effect_allele,
            'G'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing$missing.snp.data$beta,
            4
            );

        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$CHROM.pgs,
            'chr4'
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$CHROM.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$POS.pgs,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$POS.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$ID.pgs,
            'rs4'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$effect_allele,
            'G'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$beta,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$REF)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$ALT)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$Indiv)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid$missing.snp.data$gt_GT_alleles)
            );

        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$CHROM.pgs,
            'chr4'
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$CHROM.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$POS.pgs,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$POS.vcf)
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$ID.pgs,
            'rs4'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$effect_allele,
            'G'
            );
        expect_equal(
            test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$beta,
            4
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$REF)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$ALT)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$Indiv)
            );
        expect_true(
            is.na(test.combine.vcf.with.pgs.missing.locus.matching.rsid.with.semicolons$missing.snp.data$gt_GT_alleles)
            );

        }
    );

test_that(
    'combine.vcf.with.pgs correctly performs merge with no rsID', {
        test.vcf.data <- data.frame(
            CHROM = c('chr1', 'chr2', 'chr3', 'chr4'),
            POS = c(1, 2, 3, 4),
            ID = c('rs1', 'rs2', 'rs3', 'rs4'),
            REF = c('A', 'T', 'C', 'G'),
            ALT = c('T', 'A', 'G', 'C'),
            Indiv = c('sample1', 'sample2', 'sample3', 'sample4'),
            gt_GT_alleles = c('A/T', 'T/A', 'C/G', 'G/C')
            );
        test.vcf.data.missing <- data.frame(
            CHROM = c('chr1', 'chr3', 'chr2', 'chr4'),
            POS = c(1, 3, 2, 5),
            ID = c('rs1', 'rs3', 'rs2', 'rs5'),
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

        test.combine.vcf.with.pgs <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data,
            pgs.weight.data = test.pgs.weight.data
            );

        test.combine.vcf.with.pgs.missing <- combine.vcf.with.pgs(
            vcf.data = test.vcf.data.missing,
            pgs.weight.data = test.pgs.weight.data
            );

        # check that combine.vcf.with.pgs returns the correct number of rows
        expect_equal(
            nrow(test.combine.vcf.with.pgs$merged.vcf.with.pgs.data),
            4
            );
        expect_equal(
            nrow(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            4
            );

        # check that combine.vcf.with.pgs returns the correct number of columns
        expect_equal(
            ncol(test.combine.vcf.with.pgs$merged.vcf.with.pgs.data),
            10
            );
        expect_equal(
            ncol(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data),
            10
            );

        # check that combine.vcf.with.pgs returns the correct columns
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'ID', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs$merged.vcf.with.pgs.data))
            );
        expect_true(
            all(c('CHROM', 'POS', 'REF', 'ALT', 'ID', 'Indiv', 'gt_GT_alleles', 'effect_allele', 'beta', 'merge.strategy') %in% colnames(test.combine.vcf.with.pgs.missing$merged.vcf.with.pgs.data))
            );

    }
)

test_that(
    'combine.vcf.with.pgs works correctly with real data', {
        vcf.import.split.name <- 'split.wide.vcf.matrices';
        vcf.import.long.name <- 'combined.long.vcf.df';

        test.vcf.data <- import.vcf('data/HG001_GIAB.vcf.gz', long.format = TRUE);
        vcf.long.format <- test.vcf.data[[vcf.import.long.name]]$dat;
        vcf.wide.format <- test.vcf.data[[vcf.import.split.name]];
        test.pgs.weight.data <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');

        expect_no_error(
            combine.vcf.with.pgs(
                vcf.data = vcf.long.format,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            );
        expect_no_error(
            combine.vcf.with.pgs(
                vcf.data = vcf.wide.format$vcf.fixed.fields,
                pgs.weight.data = test.pgs.weight.data$pgs.weight.data
                )
            );
        }
    )
