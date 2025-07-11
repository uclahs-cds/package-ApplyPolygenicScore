# function that converts test data in long format into wide format

convert.long.vcf.to.wide.vcf <- function(long.vcf) {
    fixed.colnames <- c('CHROM', 'POS', 'ID', 'REF', 'ALT');
    fixed.data <- long.vcf[, fixed.colnames];
    fixed.data$variant.id <- paste0(fixed.data$CHROM, ':', fixed.data$POS);
    fixed.data <- unique(fixed.data); # remove duplicates
    #fixed.data$allele.matrix.row.index <- seq_len(nrow(fixed.data));

    variant.id <- paste0(long.vcf$CHROM, ':', long.vcf$POS);
    allele.matrix <- as.matrix(get.variant.by.sample.matrix(
        long.data = long.vcf,
        variant.id = variant.id,
        value.var = 'gt_GT_alleles'
        ));

    # sort allele matrix in order of corresponding fixed data variant id
    allele.matrix <- allele.matrix[match(fixed.data$variant.id, rownames(allele.matrix)), ];

    # save allele matrix row indices in fixed data
    fixed.data$allele.matrix.row.index <- seq_len(nrow(allele.matrix));

    # format output according to import.vcf standards
    output <- list(
        genotyped.alleles = allele.matrix,
        vcf.fixed.fields = fixed.data
        );

    return(output);

    }

initialize.defaults <- function() {
    phenotype.data = NULL
    phenotype.analysis.columns = NULL
    correct.strand.flips = TRUE
    remove.ambiguous.allele.matches = FALSE
    max.strand.flips = 0
    remove.mismatched.indels = FALSE
    output.dir = NULL
    file.prefix = NULL
    missing.genotype.method = 'mean.dosage'
    use.external.effect.allele.frequency = FALSE
    n.percentiles = NULL
    analysis.source.pgs = NULL
    validate.inputs.only = FALSE
    }