#' @title Merge VCF with PGS
#' @description Match PGS SNPs to corresponding VCF information by genomic coordinates via merge.
#' The merge is a left outer join, meaning that all PGS SNPs will be kept even if they are missing from the VCF, and all VCF SNPs that are not a component of the PGS will be dropped.
#' @param vcf.data A data.frame containing VCF data.
#' @param pgs.weight.data A data.frame containing PGS data.
#' @return A list containing a data.frame of merged VCF and PGS data and a data.frame of PGS SNPs missing from the VCF.
#' @examples
#' # Example VCF
#' vcf.path <- system.file(
#'  'extdata',
#'  'HG001_GIAB.vcf.gz',
#'  package = 'ApplyPolygenicScore',
#'  mustWork = TRUE
#'  );
#' vcf.import <- import.vcf(vcf.path);
#'
#' # Example pgs weight file
#' pgs.weight.path <- system.file(
#'  'extdata',
#'  'PGS000662_hmPOS_GRCh38.txt.gz',
#'  package = 'ApplyPolygenicScore',
#'  mustWork = TRUE
#'  );
#' pgs.weight.import <- import.pgs.weight.file(pgs.weight.path);
#'
#' merge.data <- merge.vcf.with.pgs(
#'  vcf.data = vcf.import$dat,
#'  pgs.weight.data = pgs.import$pgs.weight.data,
#'  );
#' @export
merge.vcf.with.pgs <- function(vcf.data, pgs.weight.data) {

    # check that inputs are data.frames
    if (!is.data.frame(vcf.data)) {
        stop('vcf.data must be a data.frame');
        }

    if (!is.data.frame(pgs.weight.data)) {
        stop('pgs.weight.data must be a data.frame');
        }

    # check that inputs contain required columns for mergeing
    required.vcf.columns <- c('CHROM', 'POS');
    required.pgs.columns <- c('CHROM', 'POS');
    if (!all(required.vcf.columns %in% colnames(vcf.data))) {
        stop('vcf.data must contain columns named CHROM and POS');
        }
    if (!all(required.pgs.columns %in% colnames(pgs.weight.data))) {
        stop('pgs.weight.data must contain columns named CHROM and POS');
        }

    # match CHROM notation in pgs data to vcf data formatting
    chr.prefix <- grepl('^chr', vcf.data$CHROM[1]);
    numeric.sex.chr <- any(grepl('23$', vcf.data$CHROM));

    pgs.weight.data$CHROM <- ApplyPolygenicScore::format.chromosome.notation(
        chromosome = pgs.weight.data$CHROM,
        chr.prefix = chr.prefix,
        numeric.sex.chr = numeric.sex.chr
        );

    # "left outer join" on CHROM and POS columns
    # this merge keeps all pgs SNPs even if they are missing from the VCF
    # and drops all VCF SNPs that are not a component of the pgs
    merged.vcf.with.pgs.data <- merge(
        x = pgs.weight.data,
        y = vcf.data,
        by = c('CHROM', 'POS'),
        all.x = TRUE
        );

    # check for pgs SNPs missing from the VCF data
    missing.pgs.snp.index <- is.na(merged.vcf.with.pgs.data$REF);
    if (any(missing.pgs.snp.index)) {
        warning(paste('PGS is missing', sum(missing.pgs.snp.index)), ' SNPs from VCF');
        missing.snp.data <- merged.vcf.with.pgs.data[missing.pgs.snp.index, ];
        } else {
        missing.snp.data <- NULL;
        }

    output <- list(
        merged.vcf.with.pgs.data = merged.vcf.with.pgs.data,
        missing.snp.data = missing.snp.data
        );

    return(output);
    }
