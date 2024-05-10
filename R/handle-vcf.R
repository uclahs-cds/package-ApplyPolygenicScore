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

check.vcf.for.split.multiallelic.sites <- function(vcf.vcfR) {
    # check for duplicate CHROM:POS entries
    site.coordinates <- paste(vcf.vcfR@fix[ ,'CHROM'], vcf.vcfR@fix[ ,'POS'], sep = ':');
    split.site.coordinates.index <- duplicated(site.coordinates);
    if (any(split.site.coordinates.index)) {
        stop(
            paste0('Split multiallelic site detected at ', site.coordinates[split.site.coordinates.index], '. Please merge multiallelic sites before importing.')
            );
        }
    }

#' @title Import VCF file
#' @description Import a VCF file into a long form data.frame using vcfR package.
#' @param vcf.path A character string indicating the path to the VCF file to be imported.
#' @param info.fields A character vector indicating the INFO fields to be imported.
#' @param format.fields A character vector indicating the FORMAT fields to be imported.
#' @param verbose A logical indicating whether verbose output should be printed by vcfR.
#' @return A list containing a tibble of VCF meta data and a tibble data.frame containing the parsed VCF data in long form.
#' @examples
#' # Example VCF
#' vcf <- system.file(
#'  'extdata',
#'  'HG001_GIAB.vcf.gz',
#'  package = 'ApplyPolygenicScore',
#'  mustWork = TRUE
#'  );
#' vcf.data <- import.vcf(vcf.path = vcf);
#' @export
import.vcf <- function(vcf.path, info.fields = NULL, format.fields = NULL, verbose = FALSE) {
    # import VCF file with vcfR
    vcf.vcfR <- vcfR::read.vcfR(file = vcf.path, convertNA = TRUE, verbose = verbose);

    # check for no INFO fields vcfR bug
    vcf.vcfR <- check.for.no.info.fields(vcf.vcfR);

    # check for split multiallelic sites
    check.vcf.for.split.multiallelic.sites(vcf.vcfR);

    # convert to long form tibble dataframe w/ vcfR
    long.vcf <- vcfR::vcfR2tidy(
        x = vcf.vcfR,
        single_frame = TRUE,
        info_fields = info.fields,
        format_fields = format.fields,
        verbose = verbose
        );

    return(long.vcf);
    }
