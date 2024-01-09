# a utility function for adding slop to BED coordinates
add.slop <- function(bed, slop) {
    bed$start <- bed$start - slop;
    bed$end <- bed$end + slop;

    # check for negative start coordinates, replace with 0, and issue a warning
    if (any(bed$start < 0)) {
        bed$start[bed$start < 0] <- 0;
        warning('Slop caused negative start coordinates; replacing with 0.');
        }

    return(bed);
    }

#' @title Convert PGS data to BED format
#' @description Convert imported and formatted PGS compnent SNP coordinate data to BED format.
#' @param pgs.weight.data A data.frame containing SNP coordinate data with standardized CHROM and POS columns.
#' @param chr.prefix A logical indicating whether the 'chr' prefix should be used when formatting chromosome name.
#' @param numeric.sex.chr A logical indicating whether the sex chromosomes should be formatted numerically, as opposed to alphabetically.
#' @param slop An integer indicating the number of base pairs to add to the BED interval on either side.
#' @return A data.frame containing the PGS component SNP coordinate data in BED format and any other columns provided in pgs.weight.data.
#' @export
convert.pgs.to.bed <- function(pgs.weight.data, chr.prefix = TRUE, numeric.sex.chr = FALSE, slop = 0) {
    # check that data is a data.frame
    if (!is.data.frame(pgs.weight.data)) {
        stop('data must be a data.frame');
        }

    # check that data has CHROM and POS columns
    if (!all(c('CHROM', 'POS') %in% colnames(pgs.weight.data))) {
        stop('data must have CHROM and POS columns');
        }

    # check that slop is a non-negative integer
    if (!is.numeric(slop) | slop < 0) {
        stop('slop must be a non-negative integer');
        }

    # convert CHROM to default format (no 'chr' prefix, alphabetic sex chromosomes)
	pgs.weight.data$CHROM <- gsub('chr', '', pgs.weight.data$CHROM);
    pgs.weight.data$CHROM <- gsub('23', 'X', pgs.weight.data$CHROM);
    pgs.weight.data$CHROM <- gsub('24', 'Y', pgs.weight.data$CHROM);

    # apply requested CHROM formatting
    if (chr.prefix) {
        pgs.weight.data$CHROM <- paste0('chr', pgs.weight.data$CHROM);
        }

    if (numeric.sex.chr) {
        pgs.weight.data$CHROM <- gsub('X', '23', pgs.weight.data$CHROM);
        pgs.weight.data$CHROM <- gsub('Y', '24', pgs.weight.data$CHROM);
        }

    ## assemble BED file ##
	# 0-index coordinates
	pgs.bed <- data.frame(
		chr = pgs.weight.data$CHROM,
		start = pgs.weight.data$POS - 1,
		end = pgs.weight.data$POS
		);
    # check for negative start coordinates, report an error
    if (any(pgs.bed$start < 0)) {
        stop('0-indexing caused negative start coordinates.');
        }

    # add slop
    if (slop > 0) {
        pgs.bed <- add.slop(bed = pgs.bed, slop = slop);
        }

	# concat with the rest of the prs columns
	pgs.bed <- cbind(pgs.bed, subset(pgs.weight.data, select = -c(CHROM, POS)));

	return(pgs.bed);
	}

#' @title Merge PGS BED files
#' @description Merge overlapping PGS coordinates in multiple BED files.
#' @param pgs.bed.list A named list of data.frames containing PGS coordinates in BED format.
#' @param add.annotation.data A logical indicating whether an additional annotation data column should be added to the annotation column.
#' @param annotation.column.index An integer indicating the index of the column in the data frames in pgs.bed.list that should be added to the annotation column.
#' @param slop An integer indicating the number of base pairs to add to the BED interval on either side.
#' @return A data.frame containing the merged PGS coordinates in BED format with an extra annotation column containing the name of the PGS and data from one additional column optionally selected by the user.
#' @export
merge.pgs.bed <- function(pgs.bed.list, add.annotation.data = FALSE, annotation.column.index = 4, slop = 0) {

    # check that pgs.bed.list is a named list
    if (!is.list(pgs.bed.list) | is.null(names(pgs.bed.list))) {
        stop('pgs.bed.list must be a named list');
        }

    # check that all elements of pgs.bed.list are data.frames
    if (!all(sapply(pgs.bed.list, is.data.frame))) {
        stop('all elements of pgs.bed.list must be data.frames');
        }
    
    # check that all elements of pgs.bed.list have the same column names
    if (!all(sapply(pgs.bed.list, function(x) all(colnames(x) == colnames(pgs.bed.list[[1]]))))) {
        stop('all elements of pgs.bed.list must have the same column names');
        }
    
    # check that each element of pgs.bed.list is in BED format (chr, start, end)
    if (!all(sapply(pgs.bed.list, function(x) all( c('chr', 'start', 'end') %in% colnames(x))))) {
        stop('all elements of pgs.bed.list must have columns named chr, start, and end');
        }
    
    # check that all intervals specified in pgs.bed.list are one bp in length
    if (!all(sapply(pgs.bed.list, function(x) all(x$start == x$end - 1)))) {
        stop('all intervals specified in pgs.bed.list must represent one SNP and be one bp in length');
        }


    if (add.annotation.data) {
        # check that annotation.column.index is within the range of the number of columns in the data.frames in pgs.bed.list
        if (annotation.column.index > ncol(pgs.bed.list[[1]])) {
            stop('annotation.column.index must be within the range of the number of columns in the data.frames in pgs.bed.list');
            }
    }

    # Annotate each PGS BED file with the name of the PGS
    for (i in 1:length(pgs.bed.list)) {
        pgs.bed.list[[i]]$annotation <- names(pgs.bed.list)[i];
        }

    # concatenate all BED files in list
    concatenated.bed <- do.call(rbind, pgs.bed.list);

    # add requested annotation data to annotation column

    if (add.annotation.data) {
        concatenated.bed$annotation <- paste(concatenated.bed[ , annotation.column.index], concatenated.bed$annotation, sep = '|');
        }

    # sort by chromosome and position
    concatenated.bed <- concatenated.bed[order(concatenated.bed$chr, concatenated.bed$start), ];

    # merge overlapping SNP intervals and combine overlapping annotations using aggregate
    merged.bed <- aggregate(annotation ~ chr + start + end, data = concatenated.bed, FUN = paste, collapse = ',');

    # sort by chromosome and position
    merged.bed <- merged.bed[order(merged.bed$chr, merged.bed$start), ];

    # add slop
    if (slop > 0) {
        merged.bed <- add.slop(bed = merged.bed, slop = slop);
        }

    # return merged BED file
    return(merged.bed);
    }
