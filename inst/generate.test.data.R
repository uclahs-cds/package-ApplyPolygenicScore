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

# create BED file for VCF filtration
pgs.weights <- ApplyPolygenicScore::import.pgs.weight.file(
    pgs.weigh.path <- 'tests/testthat/data/PGS003378_hmPOS_GRCh38.txt',
    use.harmonized.data = TRUE
    );
pgs.bed <- ApplyPolygenicScore::convert.pgs.to.bed(
    pgs.weight.data = pgs.weights$pgs.weight.data,
    chr.prefix = TRUE,
    numeric.sex.chr = FALSE,
    slop = 10
    );
write.table(
    x = pgs.bed,
    file = 'tests/testthat/data/PGS003378_hmPOS_GRCh38_slop10.bed',
    sep = '\t',
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
    );

# create simple PGS weight and VCF data for PGS application testing
simple.pgs.application.test.data <- list(
    pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2'),
        POS = c(1, 2),
        effect_allele = c('T', 'T'),
        beta = c(1.0, 1.0)
        ),
    vcf.data = data.frame(
        CHROM = c('chr1', 'chr1', 'chr2', 'chr2'),
        POS = c(1, 1, 2, 2),
        REF = c('A', 'A', 'T', 'T'),
        ALT = c('T', 'T', 'A', 'A'),
        Indiv = c('sample1', 'sample2', 'sample1', 'sample2'),
        gt_GT_alleles = c('A/A', 'T/T', 'A/T', 'T/A')
        )
    );
save(
    simple.pgs.application.test.data,
    file = 'tests/testthat/data/simple.pgs.application.test.data.Rda'
    );

# create simple VCF data for testing multiallelic site handling
merged.multiallelic.site.test.data <- list(
    merged.multiallelic.vcf.data = data.frame(
        CHROM = c('chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2', 'chr3', 'chr3', 'chr3'),
        # merged multiallelic site at chr2:2 with betas provided
        # merged multiallelic site at chr3:3 with no betas provided
        POS = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
        REF = c('T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T'),
        # three possible alleles at chr2:2 (T, A, C)
        ALT = c('A', 'A', 'A', 'A,C', 'A,C', 'A,C', 'A,G', 'A,G', 'A,G'),
        Indiv = c('sample1', 'sample2', 'sample3', 'sample1', 'sample2', 'sample3', 'sample1', 'sample2', 'sample3'),
        # multiallelic genotypes at chr2:2
        # sample 1 is heterozygous for alt allele A
        # sample 2 is heterozygous for alt allele C
        # sample 3 is homozygous for alt allele C
        # multiallelic genotypes at chr3:3
        # sample 1 is heterozygous for alt allele A
        # sample 2 is homozygous for alt allele G
        # sample 3 is homozygous for alt allele G
        gt_GT_alleles = c('T/T', 'T/A', 'A/A', 'T/A', 'T/C', 'C/C', 'T/A', 'T/G', 'G/G')
        ),
    ref.as.single.risk.allele.multiallelic.pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2', 'chr3'),
        POS = c(1, 2, 3),
        effect_allele = c('T', 'T', 'T'),
        beta = c(1.0, 1.0, 1.0)
        ),
    alt.as.single.risk.allele.multiallelic.pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2', 'chr3'),
        POS = c(1, 2, 3),
        effect_allele = c('A', 'A', 'A'),
        beta = c(1.0, 1.0, 1.0)
        ),
    alt.as.two.risk.alleles.multiallelic.pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2', 'chr2', 'chr3'),
        POS = c(1, 2, 2, 3),
        effect_allele = c('A', 'A', 'C', 'A'),
        beta = c(1.0, 1.0, 0.5, 1.0)
        ),
    ref.and.alt.as.two.risk.alelles.multiallelic.pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2', 'chr2', 'chr3'),
        POS = c(1, 2, 2, 3),
        effect_allele = c('A', 'A', 'T', 'A'),
        beta = c(1.0, 1.0, 0.5, 1.0)
        )
    );

save(
    merged.multiallelic.site.test.data,
    file = 'tests/testthat/data/merged.multiallelic.site.test.data.Rda'
    );

# create simple VCF data for testing missing site handling
missing.genotype.test.data <- list(
    missing.genotype.vcf.data = data.frame(
        CHROM = c('chr1', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2', 'chr2', 'chr3', 'chr3', 'chr3', 'chr3'),
        POS = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
        REF = c('T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T'),
        ALT = c('A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'),
        Indiv = c('sample1', 'sample2', 'sample3', 'sample4', 'sample1', 'sample2', 'sample3', 'sample4', 'sample1', 'sample2', 'sample3', 'sample4'),
        gt_GT_alleles = c('T/T', 'T/A', './.', 'T/A', 'T/A', 'A/A', '.', 'T/A', 'T/T', 'T/A', NA, NA)
        ),
    missing.genotype.pgs.weight.data = data.frame(
        CHROM = c('chr1', 'chr2', 'chr3', 'chr4'),
        POS = c(1, 2, 3, 4),
        effect_allele = c('A', 'A', 'A', 'A'),
        beta = c(1.0, 1.0, 1.0, 1.0)
        )
    );
save(
    missing.genotype.test.data,
    file = 'tests/testthat/data/missing.genotype.test.data.Rda'
    );

# create data for testing phenotype related functionality
# generate random genotype data
set.seed(123);
n.samples <- 10;
n.variants <- 10;
genotypes <- rmultinom(n = n.samples, size = 2, prob = seq(0, 1, length.out = n.variants));
genotypes <- as.vector(t(genotypes));
# convert to ref and alt letter alleles
ref.alleles <- rep('T', length(genotypes));
ref.alleles[genotypes == 2] <- 'A';
alt.alleles <- rep('T', length(genotypes));
alt.alleles[genotypes == 1] <- 'A';
alt.alleles[genotypes == 2] <- 'A';


phenotype.test.data <- list(
    # vcf with 10 samples and 10 variants
    vcf.data = data.frame(
        CHROM = paste0('chr', rep(1:10, each = 10)),
        POS = rep(1:10, each = 10),
        REF = rep('T', 100),
        ALT = rep('A', 100),
        Indiv = rep(paste0('sample', 1:n.samples), n.variants),
        gt_GT_alleles = paste0(ref.alleles, '/', alt.alleles)
        ),
    pgs.weight.data = data.frame(
        CHROM = paste0('chr', 1:10),
        POS = 1:10,
        effect_allele = 'A',
        beta = rnorm(10) # random beta values
        ),
    phenotype.data = data.frame(
        sample = paste0('sample', 1:n.samples),
        continuous.phenotype = rnorm(n.samples), # random phenotype values
        binary.phenotype = rbinom(n.samples, size = 1, prob = 0.5) # random binary phenotype values
        )
    );
save(
    phenotype.test.data,
    file = 'tests/testthat/data/phenotype.test.data.Rda'
    );
