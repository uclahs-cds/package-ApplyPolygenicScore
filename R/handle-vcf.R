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
#' @param long.format A logical indicating whether the VCF import should be converted into long format (one row per sample-variant combination)
#' @param info.fields A character vector indicating the INFO fields to be imported, only applicable when long format is \code{TRUE}.
#' @param format.fields A character vector indicating the FORMAT fields to be imported, only applicable when long format is \code{TRUE}.
#' @param verbose A logical indicating whether verbose output should be printed by vcfR.
#' @return A list of two elements containing imported VCF information in wide format and in long format if requested.
#'
#' \strong{Output Structure}
#'
#' The outputed list contains the following elements:
#' \itemize{
#' \item split.wide.vcf.matrices: A list with two elements: a \code{data.table} of fixed VCF fields and a \code{matrix} of genotyped alleles.
#' \item combined.long.vcf.df: Default is \code{NULL} otherwise if \code{long.format == TRUE} a list with two elements inherited from vcfR: a data frame meta data from the VCF header and a data frame of all requested VCF fields (including INFO and FORMAT fields) in long format. Number of rows is equal to the number of samples times the number of sites in the VCF.
#' }
#'
#' The \code{split.wide.vcf.matrices} list contains the following elements:
#' \itemize{
#' \item genotyped.alleles: A matrix of genotyped alleles (e.g. "A/C"). Rows are unique sites and columns are unique samples in the input VCF.
#' \item vcf.fixed.fields: A data table of the following fixed (not varying by sample) VCF fields: CHROM, POS, ID, REF, ALT. Also one additional column \code{allele.matrix.row.index} indicating the corresponding row in \code{genotyped.alleles}
#' }
#'
#' The \code{combined.long.vcf.df} list contains the following elements:
#' \itemize{
#' \item meta: A data frame of meta data parsed from the VCF header
#' \item dat: A data frame of all default VCF fields and all requested INFO and FORMAT fields in long format. Number of rows is equal to the number of unique samples times the number of unique sites in the VCF.
#' }
#'
#' The wide format is intended to efficiently contain the bare minimum information required for PGS application.
#' It intentionally excludes much of the additional information included in a typical VCF, and splits off genotypes into a separate matrix for easy manipulation.
#' If users wish to maintain additional information in the INFO and FORMAT fields for e.g. variant filtering, the long format allows this.
#' However, the long format requires substantially more memory to store, and is not recommended for large input files.
#'
#' @examples
#' # Example VCF
#' vcf <- system.file(
#'     'extdata',
#'     'HG001_GIAB.vcf.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' vcf.data <- import.vcf(vcf.path = vcf, long.format = TRUE);
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
            long.vcf <- NULL;
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
