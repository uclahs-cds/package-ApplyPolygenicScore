
# function that checks for the presence of required columns in a PGS weight file
check.pgs.weight.columns <- function(x, harmonized = TRUE) {
    required.generic.columns <- c('chr_name', 'chr_position', 'effect_allele', 'effect_weight');
    required.harmonized.columns <- c('hm_chr', 'hm_pos');

    if (harmonized) {
        required.columns <- c(required.generic.columns, required.harmonized.columns);
        } else {
        required.columns <- required.generic.columns;
        }

    if (!all(required.columns %in% x)) {
        stop('The following required columns are missing from the PGS weight file: ', paste(setdiff(required.columns, x), collapse = ', '));
        }
    
    return(TRUE);
    }


# function for parsing metadata from a file header that is indicated by '#' or '##'
parse.pgs.input.header <- function(input) {
    # check if file is zipped
    if (grepl('.gz$', input)) {
        # unzip file
        input.connection <- gzfile(input);
        } else {
        # open file connection
        input.connection <- file(input);
        }

    file <- readLines(input);

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
import.pgs.weight.file <- function(input, use.harmonized.data = TRUE) {

    # parse file header
    file.metadata <- parse.pgs.input.header(input = input);

    # check if file is zipped
    if (grepl('.gz$', input)) {
        # unzip file
        input.connection <- gzfile(input);
        } else {
        # open file connection
        input.connection <- file(input);
        }
    
    # read in data frame
    pgs.weight.data <- read.table(
        file = input.connection,
        header = TRUE,
        sep = '\t',
        comment.char = '#',
        stringsAsFactors = FALSE
        );

    # check that required columns are present
    check.pgs.weight.columns(x = colnames(pgs.weight.data), harmonized = use.harmonized.data);

    # check if file is harmonized and format columns accordingly
    if (use.harmonized.data) {

        # label harmonized data columns with standardized names
        pgs.weight.data$CHROM <- pgs.weight.data$hm_chr;
        pgs.weight.data$POS <- pgs.weight.data$hm_pos;
        format.harmonized.columns = TRUE
        } else {

        # label non-harmonized data columns with standardized names
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
