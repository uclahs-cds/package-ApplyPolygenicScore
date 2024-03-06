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

    # check that output.path is a directory
    if (!is.character(output.path)) {
        stop('output.path must be a character string');
        }

    # check that file.prefix is a character string
    if (!is.null(file.prefix) && !is.character(file.prefix)) {
        stop('file.prefix must be a character string');
        }

    # generate filename
    filename.for.per.sample.pgs.summary.data <- BoutrosLab.utilities::generate.filename(
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
