#' Genome In a Bottle (GIAB) VCF formatted data.
#'
#' A subset of the genotype callset provided by the GIAB consortium for public use as a benchmarking dataset.
#' These files contain genotype information from the NA12878/HG001 HapMap individual.
#' These variants were aligned and called against the GRCh38 reference genome.
#' To provide a more manageable dataset, we have subsetted the variants to the intersection of the GIAB callset
#' and the harmonized component variants of the prostate cancer polygenic score PGS PGS003378 (accessible through the PGS Catalog)
#' plus 10bp of flanking sequence on either side of the variant.
#' Typically, a cohort intended for polygenic scoring is combined into one VCF, thus we have added a second sample to the VCF by
#' merging the single-sample GIAB VCF with itself, producing genotype data for HG001 and its duplicate HG001.2.
#'
#' @format ## `*_GIAB.vcf.gz`
#' g-zipped VCF file containing 234 lines of metadata in the header indicated by `##`, one line of column header
#' indicated by `#` and 63 lines of genetic information. The VCF has 9 standard columns and 1 or 2 additional sample columns.
#'
#' \describe{
#'  \item{CHROM}{Chromosome}
    #' \item{POS}{Position}
    #' \item{ID}{Variant ID}
    #' \item{REF}{Reference allele}
    #' \item{ALT}{Alternate allele}
    #' \item{QUAL}{Quality score}
    #' \item{FILTER}{Filter status}
    #' \item{INFO}{Variant information}
    #' \item{FORMAT}{Genotype format}
    #' \item{HG001}{Genotype data for sample HG001}
    #' \item{HG001.2}{Genotype data for sample HG001.2, a duplicate of HG001}
#' }
#' @source Genome in a Bottle Consortium. https://www.nist.gov/programs-projects/genome-bottle
#' FTP: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/
"GIAB_VCF"
