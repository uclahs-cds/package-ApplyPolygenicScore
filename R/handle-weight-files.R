#' @title Check PGS weight file columns
#' @description Check that a PGS weight file contains the required columns for PGS application with \code{apply.polygenic.score}.
#' @param pgs.weight.colnames A character vector of column names.
#' @param harmonized A logical indicating whether the presence of harmonized columns should be checked.
#' @return A logical indicating whether the file contains the required columns.
#' @export
check.pgs.weight.columns <- function(pgs.weight.colnames, harmonized = TRUE) {
    required.generic.columns <- c('chr_name', 'chr_position', 'effect_allele', 'effect_weight');
    required.harmonized.columns <- c('hm_chr', 'hm_pos');

    if (harmonized) {
        required.columns <- c(required.generic.columns, required.harmonized.columns);
        } else {
        required.columns <- required.generic.columns;
        }

    if (!all(required.columns %in% pgs.weight.colnames)) {
        stop('The following required columns are missing from the PGS weight file: ', paste(setdiff(required.columns, pgs.weight.colnames), collapse = ', '));
        }

    rsid.columns <- c('rsID', 'hm_rsID');
    if (any(!(rsid.columns %in% pgs.weight.colnames))) {
        warning('rsID or hm_rsID column not found in PGS weight file. Merging by rsID will not be possible.');
        }

    return(TRUE);
    }

# utility function for checking if a file is gzipped and opening a connection accordingly
open.input.connection <- function(input) {
    # check if file is zipped
    if (grepl('.gz$', input)) {
        input.connection <- gzfile(input);
        } else {
        input.connection <- file(input);
        }

    return(input.connection);
    }

#' @title Parse PGS input file header
#' @description Parse metadata from a PGS input file header.
#' @param pgs.weight.path A character string indicating the path to the pgs weight file.
#' @return A data frame containing the metadata from the file header.
#' @examples
#' # Example pgs weight file
#' pgs.weight.path <- system.file(
#'     'extdata',
#'     'PGS000662_hmPOS_GRCh38.txt.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' parse.pgs.input.header(pgs.weight.path);
#' @export
parse.pgs.input.header <- function(pgs.weight.path) {

    # check that pgs.weight.path exists
    if (!file.exists(pgs.weight.path)) {
        stop(paste0(pgs.weight.path, ' does not exist.'));
        }

    # open file connection
    input.connection <- open.input.connection(input = pgs.weight.path);

    file <- readLines(input.connection);

    # close file connection
    close(input.connection);

    # identify file header
    file.header <- file[grepl('^#', file)];
    # remove leading '#' from file header
    file.header <- gsub('^#', '', file.header);
    file.header <- trimws(file.header);

    # remove all lines that still have leading '#'
    file.header <- file.header[!grepl('^#', file.header)];

    # split each line on '=' and convert to data frame
    file.header <- as.data.frame(do.call(rbind, strsplit(file.header, '=')));
    colnames(file.header) <- c('key', 'value');
    return(file.header);
    }

# function for importing a PGS weight file formatted according to PGS catalog guidelines
#' @title Import PGS weight file
#' @description Import a PGS weight file formatted according to PGS catalog guidelines, and prepare for PGS application with \code{apply.polygenic.score()}.
#' @param pgs.weight.path A character string indicating the path to the pgs weight file.
#' @param use.harmonized.data A logical indicating whether the file should be formatted to indicate harmonized data columns for use in future PGS application.
#' @return A list containing the file metadata and the weight data.
#' @examples
#' # Example pgs weight file
#' pgs.weight.path <- system.file(
#'     'extdata',
#'     'PGS000662_hmPOS_GRCh38.txt.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' import.pgs.weight.file(pgs.weight.path);
#'
#' # Note, harmonized data is used by default. To disable set `use.harmonized.data = FALSE`
#' import.pgs.weight.file(pgs.weight.path, use.harmonized.data = FALSE);
#' @export
import.pgs.weight.file <- function(pgs.weight.path, use.harmonized.data = TRUE) {

    # check that pgs.weight.path exists
    if (!file.exists(pgs.weight.path)) {
        stop(paste0(pgs.weight.path, ' does not exist.'));
        }

    # parse file header
    file.metadata <- parse.pgs.input.header(pgs.weight.path = pgs.weight.path);

    # open file connection
    input.connection <- open.input.connection(input = pgs.weight.path);

    # read in data frame
    pgs.weight.data <- read.table(
        file = input.connection,
        header = TRUE,
        sep = '\t',
        comment.char = '#',
        stringsAsFactors = FALSE
        );

    # check that required columns are present
    check.pgs.weight.columns(pgs.weight.colnames = colnames(pgs.weight.data), harmonized = use.harmonized.data);

    # check for duplicate variants
    if (any(duplicated(paste(pgs.weight.data$chr_name, pgs.weight.data$chr_position, pgs.weight.data$effect_allele, sep = ':')))) {
        stop('Duplicate variants are present in the PGS weight file. Please remove duplicate variants.');
        }

    # check for multiallelic variants
    if (any(duplicated(paste(pgs.weight.data$chr_name, pgs.weight.data$chr_position, sep = ':')))) {
        warning('Duplicate variants with different effect alleles detected in the PGS weight file. Please ensure these are true multiallelic sites.');
        }

    # check if file is harmonized and format columns accordingly
    if (use.harmonized.data) {

        # label harmonized data columns with standardized names
        # if harmonized rsID column is present and not empty, use it
        if ('hm_rsID' %in% colnames(pgs.weight.data)) {
            if (!all(is.na(pgs.weight.data$hm_rsID))) {
                pgs.weight.data$ID <- pgs.weight.data$hm_rsID;
                }
            }

        pgs.weight.data$CHROM <- pgs.weight.data$hm_chr;
        pgs.weight.data$POS <- pgs.weight.data$hm_pos;

        } else {

        # label non-harmonized data columns with standardized names
        if ('rsID' %in% colnames(pgs.weight.data)) {
            if (!all(is.na(pgs.weight.data$rsID == ''))) {
                pgs.weight.data$ID <- pgs.weight.data$rsID;
                }
            }

        pgs.weight.data$CHROM <- pgs.weight.data$chr_name;
        pgs.weight.data$POS <- pgs.weight.data$chr_position;

        }

    # extract weight format from file metadata key 'weight_type'
    weight.format <- file.metadata$value[file.metadata$key == 'weight_type'];

    # check if weight format is provided (NR by default)
    if (weight.format == 'NR') {
        # report a warning that weight format was Not Reported
        warning('Weight format was not reported in the PGS file header. Assuming beta weights.');
        }

    if (weight.format == 'NR' | grepl('beta', weight.format, ignore.case = TRUE)) {
        pgs.weight.data$beta <- as.numeric(pgs.weight.data$effect_weight);

        } else if (grepl('OR|HR', weight.format, ignore.case = TRUE)) {
        # if weight format is OR/HR, convert to numeric and then to beta
        pgs.weight.data$beta <- log(as.numeric(pgs.weight.data$effect_weight));
        # report a warning that OR/HR weights were converted to beta
        warning('OR/HR weights were converted to beta weights.');

        } else {
        # if weight format is not recognized, throw an error
        stop('Weight format is not recognized. Please specify whether weights are betas or OR/HR.');
        }

    result <- list(
        file.metadata = file.metadata,
        pgs.weight.data = pgs.weight.data
        );

    return(result);
    }
