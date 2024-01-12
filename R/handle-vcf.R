# A utility wrapper function for vcfR::vcfR2tidy() that handles a bug in vcfR
# that causes an error when the VCF file contains no INFO fields.
get.long.form.vcf.vcfR <- function(vcf.vcfR, info.fields = NULL, format.fields = NULL) {
  # check for vcfR empty INFO bug:
  INFO.meta.lines <- grepl('^##INFO', vcf.vcfR@meta);
  if (any(!INFO.meta.lines)) {
      dummy.INFO.line <- '##INFO=<ID=AF>';
      vcf.vcfR@meta[length(vcf.vcfR@meta) + 1] <- dummy.INFO.line;
      }

  # convert to long form dataframe w/ vcfR
  long.vcf <- vcfR::vcfR2tidy(
      x = vcf.vcfR,
      single_frame = TRUE,
      info_fields = info.fields,
      format_fields = format.fields,
      );

  return(long.vcf);
  }

#' @title Import VCF file
#' @description Import a VCF file into a long form data.frame using vcfR package.
#' @param vcf.path A character string indicating the path to the VCF file to be imported.
#' @param info.fields A character vector indicating the INFO fields to be imported.
#' @param format.fields A character vector indicating the FORMAT fields to be imported.
#' @param verbose A logical indicating whether verbose output should be printed by vcfR.
#' @return A long form data.frame containing the VCF data.
#' @export
import.vcf <- function(vcf.path, info.fields = NULL, format.fields = NULL, verbose = FALSE) {
  # import VCF file with vcfR
  vcf.vcfR <- vcfR::read.vcfR(file = vcf.path, convertNA = TRUE, verbose = verbose);

  # convert to long form dataframe w/ vcfR
  long.vcf <- get.long.form.vcf.vcfR(
      vcf.vcfR = vcf.vcfR,
      info.fields = info.fields,
      format.fields = format.fields
      );

  return(long.vcf);
  }
