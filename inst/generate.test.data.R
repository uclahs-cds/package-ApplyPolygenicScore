# create data frame of imported PGS weight file header (for testing)
PGS00662.metadata <- data.frame(
    key = c(
        'format_version',
        'pgs_id',
        'pgs_name',
        'trait_reported',
        'trait_mapped',
        'trait_efo',
        'genome_build',
        'variants_number',
        'weight_type',
        'pgp_id',
        'citation',
        'HmPOS_build',
        'HmPOS_date',
        'HmPOS_match_chr',
        'HmPOS_match_pos'
        ),
    value = c(
        '2.0',
        'PGS000662',
        'GRS.PCa.269',
        'Prostate Cancer',
        'prostate carcinoma',
        'EFO_0001663',
        'GRCh37',
        '269',
        'beta',
        'PGP000122',
        'Conti DV, Darst BF et al. Nat Genet (2020). doi:10.1038/s41588-020-00748-0',
        'GRCh38',
        '2022-07-29',
        '{"True": null, "False": null}',
        '{"True": null, "False": null}'
        )
    );


# save .Rda file of import function test data
import.test.data <- list(
    PGS00662.metadata = PGS00662.metadata
    );

save(
    import.test.data,
    file = 'tests/testthat/data/import.test.data.Rda'
    );

# create data frame of test PGS coordinates (for testing)
tiny.pgs.test.data <- data.frame(
    CHROM = c('1', '2', 'chr3', 'chr4', 'chr23', 'chr24', 'chrX', 'chrY'),
    POS = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7),
    effect_allele = c('A', 'T', 'C', 'G', 'A', 'T', 'C', 'G'),
    effect_weight = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, NA, NA),
    beta = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, NA, NA)
    );

# save .Rda file of tiny PGS test data
save(
    tiny.pgs.test.data,
    file = 'tests/testthat/data/tiny.pgs.test.data.Rda'
    );

# create list of data frames of test BED coordinates (for testing)
tiny.bed.test.data <- simple.test.input <- list(
    pgs1 = data.frame(
        chr = c('chr1', 'chr1', 'chr2', 'chr2'),
        start = c(1, 2, 3, 4),
        end = c(2, 3, 4, 5),
        overlap = c('no overlap with pgs2', 'overlap with pgs2', 'overlap with pgs2', 'no overlap with pgs2')
        ),
    pgs2 = data.frame(
        chr = c('chr1', 'chr1', 'chr2', 'chr3'),
        start = c(5, 2, 3, 4),
        end = c(6, 3, 4, 5),
        overlap = c('no overlap with pgs1', 'overlap with pgs1', 'overlap with pgs1', 'no overlap with pgs1')
        )
    );

# save .Rda file of tiny BED test data
save(
    tiny.bed.test.data,
    file = 'tests/testthat/data/tiny.bed.test.data.Rda'
    );
