# function for parsing metadata from a file header that is indicated by '#' or '##'
parse.pgs.input.header <- function(input) {
    # check if file is zipped
    if (grepl('.gz$', input)) {
        # unzip file
        input <- gunzip(input);
        }

    file <- readLines(input);
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

    # close file connection
    close(input.connection);

    # check if file is harmonized and format columns accordingly
    if (use.harmonized.data) {
        if (any(file.metadata$file.header$key == 'HmPOS_build')) {
            # label harmonized data columns with standardized names
            pgs.weight.data$CHROM <- pgs.weight.data$hm_chr;
            pgs.weight.data$POS <- pgs.weight.data$hm_pos;
            format.harmonized.columns = TRUE
            }
        else {
            # throw an error if harmonized data is requested but not available
            stop('Harmonized data is not available for this PGS');
            }
        } else {
        # label non-harmonized data columns with standardized names
        pgs.weight.data$CHROM <- pgs.weight.data$chr_name;
        pgs.weight.data$POS <- pgs.weight.data$chr_position;
        }

    result <- list(
        pgs.weight.data = pgs.weight.data,
        file.metadata = file.metadata
        );

    return(result);

    }
