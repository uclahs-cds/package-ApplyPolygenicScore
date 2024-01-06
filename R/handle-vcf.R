### GET.LONG.FORM.VCF.VCFR ########################################################
# Wrapper for vcfR2tidy which accounts for the bug in this function that prevents
# processing of VCF files with empty INFO columns.
get.long.form.vcf.vcfR <- function(vcf.vcfR, info.fields, format.fields) {
  # check for vcfR empty INFO bug:
  INFO.meta.lines <- grepl('^##INFO', vcf.vcfR@meta);
  if (any(!INFO.meta.lines)) {
      dummy.INFO.line <- '##INFO=<ID=AF>';
      vcf.vcfR@meta[length(vcf.vcfR@meta) + 1] <- dummy.INFO.line;
      }

  # convert to long form dataframe w/ vcfR
  long.vcf <- vcfR2tidy(
      x = vcf.vcfR,
      single_frame = TRUE,
      info_fields = info.fields,
      format_fields = format.fields,
      );

  return(long.vcf);
  }

import.vcf <- function(input, info.fields, format.fields, verbose = FALSE) {
  # import VCF file with vcfR
  vcf.vcfR <- vcfR::read.vcfR(input, convertNA = TRUE, verbose = verbose);

  # convert to long form dataframe w/ vcfR
  long.vcf <- get.long.form.vcf.vcfR(
      vcf.vcfR = vcf.vcfR,
      info.fields = info.fields,
      format.fields = format.fields,
      );

  return(long.vcf);
  }
