# A utility function for vcfR::vcfR2tidy() that handles a bug in vcfR
# that causes an error when the VCF file contains no INFO fields.
check.for.no.info.fields <- function(vcf.vcfR) {
    # If no INFO fields present, add dummy INFO line
    INFO.meta.lines <- grepl('^##INFO', vcf.vcfR@meta);
    if (!any(INFO.meta.lines)) {
            dummy.INFO.line <- '##INFO=<ID=DUMMY>';
            vcf.vcfR@meta[length(vcf.vcfR@meta) + 1] <- dummy.INFO.line;
            warning('No INFO fields present in VCF file. Adding dummy INFO field "DUMMY" to avoid vcfR import errors.');
            }

    return(vcf.vcfR);
    }

# site coordinates should be CHROM:POS format
check.vcf.for.split.multiallelic.sites <- function(site.coordinates) {
    # check for duplicate CHROM:POS entries
    if (anyDuplicated(site.coordinates)) {
        dupes <- duplicated(site.coordinates);
        stop(
            paste0('Split multiallelic site detected at ', site.coordinates[dupes], '. Please merge multiallelic sites before importing.')
            );
        }
    }

#' @title Import VCF file
#' @description A wrapper for the VCF import function in the vcfR package that formats VCF data for PGS application with \code{apply.polygenic.score()}.
#' @param vcf.path A character string indicating the path to the VCF file to be imported.
#' @param info.fields A character vector indicating the INFO fields to be imported.
#' @param format.fields A character vector indicating the FORMAT fields to be imported.
#' @param verbose A logical indicating whether verbose output should be printed by vcfR.
#' @return A list containing a tibble of VCF meta data and a tibble data.frame containing the parsed VCF data in long form.
#' @examples
#' # Example VCF
#' vcf <- system.file(
#'     'extdata',
#'     'HG001_GIAB.vcf.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' vcf.data <- import.vcf(vcf.path = vcf);
#' @export
import.vcf <- function(vcf.path, long.format = FALSE, info.fields = NULL, format.fields = NULL, verbose = FALSE) {

    # check that vcf.path exists
    if (!file.exists(vcf.path)) {
        stop(paste0(vcf.path, ' does not exist.'));
        }

    # check file size and throw warning if very large
    file.size.gb <- file.size(vcf.path) / (1024^3); # convert to GB
    if (file.size.gb > 2) {
        warning('input VCF file is large, which may cause memory issues during import.')
        }

    # import VCF file with vcfR
    vcf.vcfR <- vcfR::read.vcfR(file = vcf.path, convertNA = TRUE, verbose = verbose);

    # check for split multiallelic sites
    site.coordinates <- paste0(vcf.vcfR@fix[, 'CHROM'], ':', vcf.vcfR@fix[, 'POS']);
    check.vcf.for.split.multiallelic.sites(site.coordinates);

    if (long.format == TRUE) {
        # check for no INFO fields vcfR bug
        vcf.vcfR <- check.for.no.info.fields(vcf.vcfR);
        # convert to long form tibble dataframe w/ vcfR
        long.vcf <- vcfR::vcfR2tidy(
            x = vcf.vcfR,
            single_frame = TRUE,
            info_fields = info.fields,
            format_fields = format.fields,
            verbose = verbose
            );
        } else {
            long.vcf = NULL;
        }

    # Create a sample by allele matrix
    gt.alleles <- vcfR::extract.gt(vcf.vcfR, return.alleles = TRUE);

    # Save VCF row-wise fixed information
    extraction.cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT');
    # Convert to data table due to future merging compatibility
    vcf.fixed <- data.table::as.data.table(vcf.vcfR@fix[ , extraction.cols, drop = FALSE]);

    # remove vcfR object to free memory
    rm(vcf.vcfR);

    # add column that stores corresponding row numbers in allele matrix
    vcf.fixed$allele.matrix.row.index <- seq_len(nrow(vcf.fixed));


    # assemble output object
    output <- list(
        split.wide.vcf.matrices = list(
            genotyped.alleles = gt.alleles,
            vcf.fixed.fields = vcf.fixed
            ),
        combined.long.vcf.df = long.vcf
        );

    return(output);
    }
