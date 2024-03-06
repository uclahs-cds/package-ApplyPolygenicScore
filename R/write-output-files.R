# utility for constructing a filename
generate.filename <- function(project.stem, file.core, extension, file.date = Sys.Date()) {

	# build up the filename piece-wise
	file.name <- paste(project.stem, sep = '_');
	file.name <- paste(file.name, file.core, sep = '_');
	file.name <- paste(file.name, extension, sep = '.');

	# now add the date, if requested
	if (file.date != FALSE) {
		file.name <- paste(file.date, file.name, sep = '_');
		}
	
	return(file.name);

	}

#' @title Write PGS per sample table to file
#' @description Write PGS per sample summary data table to tab separated text file.
#' @param per.sample.pgs.summary.data data.frame of PGS per sample data
#' @param output.path character string of the path to write the PGS per sample table
#' @param file.prefix character string of the file prefix to use for the PGS per sample table
#' @export
write.per.sample.pgs.table <- function(per.sample.pgs.summary.data, output.path, file.prefix = NULL) {
    # check that input is a data.frame
    if (!is.data.frame(per.sample.pgs.summary.data)) {
        stop('pgs.data must be a data.frame');
        }

    # check that file.prefix is a character string
    if (!is.null(file.prefix) && !is.character(file.prefix)) {
        stop('file.prefix must be a character string');
        }

    # check that output.path is a directory
    if (!file.exists(output.path)) {
        stop(paste('output path', output.path, 'not found'));
        }

    # generate filename
    filename.for.per.sample.pgs.summary.data <- generate.filename(
        project.stem = file.prefix,
        file.core = 'per-sample-pgs-summary',
        extension = 'txt'
        );

    # write PGS per sample data to file
    write.table(
        x = per.sample.pgs.summary.data,
        file = file.path(output.path, filename.for.per.sample.pgs.summary.data),
        sep = '\t',
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
        );
    }
