### CONVERT.PGS.TO.BED ###################################################
# convert imported and formatted PRS files to BED format
# data = data.frame, SNP info with standardized CHROM and POS columns
# chr.prefix = Bool, whether to add the 'chr' prefix to CHROM name
# chr.X = Bool, whether to use 'X' for X-chromosome notation. If false, use numeric '23'.
# slop = int, extra distance in bp to add to the interval on either side
convert.pgs.to.bed <- function(data, chr.prefix = TRUE, numeric.sex.chr = FALSE, slop = 0) {
    # check that data is a data.frame
    if (!is.data.frame(data)) {
        stop('data must be a data.frame');
        }

    # check that data has CHROM and POS columns
    if (!all(c('CHROM', 'POS') %in% colnames(data))) {
        stop('data must have CHROM and POS columns');
        }

    # check that slop is a non-negative integer
    if (!is.numeric(slop) | slop < 0) {
        stop('slop must be a non-negative integer');
        }

    # convert CHROM to default format (no 'chr' prefix, alphabetic sex chromosomes)
	data$CHROM <- gsub('chr', '', data$CHROM);
    data$CHROM <- gsub('23', 'X', data$CHROM);
    data$CHROM <- gsub('24', 'Y', data$CHROM);

    # apply requested CHROM formatting
    if (chr.prefix) {
        data$CHROM <- paste0('chr', data$CHROM);
        }

    if (numeric.sex.chr) {
        data$CHROM <- gsub('X', '23', data$CHROM);
        data$CHROM <- gsub('Y', '24', data$CHROM);
        }

    ## assemble BED file ##
	# 0-index coordinates
	prs.bed <- data.frame(
		chr = data$CHROM,
		start = data$POS - 1,
		end = data$POS
		);
    # check for negative start coordinates, report an error
    if (any(prs.bed$start < 0)) {
        stop('0-indexing caused negative start coordinates.');
        }

    # add slop
    if (slop > 0) {
        prs.bed$start <- prs.bed$start - slop;
        prs.bed$end <- prs.bed$end + slop;

        # check for negative start coordinates, replace with 0, and issue a warning
        if (any(prs.bed$start < 0)) {
            prs.bed$start[prs.bed$start < 0] <- 0;
            warning('Slop caused negative start coordinates; replacing with 0.');
            }
        }

	# concat with the rest of the prs columns
	prs.bed <- cbind(prs.bed, subset(data, select = -c(CHROM, POS)));

	return(prs.bed);
	}

### MERGE.PGS.BED ##################################################
# pgs.bed.list = named list of pgs coordinates in BED format
merge.pgs.bed <- function(pgs.bed.list, slop = 0) {

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
    if (!all(sapply(pgs.bed.list, function(x) all(colnames(x) %in% c('chr', 'start', 'end'))))) {
        stop('all elements of pgs.bed.list must have columns named chr, start, and end');
        }
    
    # check that all intervals specified in pgs.bed.list are one bp in length
    if (!all(sapply(pgs.bed.list, function(x) all(x$start == x$end - 1)))) {
        stop('all intervals specified in pgs.bed.list must represent one SNP and be one bp in length');
        }

    # Annotate each PGS BED file with the name of the PGS
    for (i in 1:length(pgs.bed.list)) {
        pgs.bed.list[[i]]$annotation <- names(pgs.bed.list)[i];
        }

    # concatenate all BED files in list
    concatenated.bed <- do.call(rbind, pgs.bed.list);

    # add snpID column to annotation
    concatenated.bed$annotation <- paste(concatenated.bed$rsID, concatenated.bed$annotation, sep = '|');

    # sort by chromosome and position
    concatenated.bed <- concatenated.bed[order(concatenated.bed$chr, concatenated.bed$start), ];

    # merge overlapping SNP intervals and combine overlapping annotations using aggregate
    merged.bed <- aggregate(annotation ~ chr + start + end, data = concatenated.bed, FUN = paste, collapse = ',');

    # sort by chromosome and position
    merged.bed <- merged.bed[order(merged.bed$chr, merged.bed$start), ];

    # add slop
    merged.bed$start <- merged.bed$start - slop;
    merged.bed$end <- merged.bed$end + slop;

    # return merged BED file
    return(merged.bed);
    }
