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

# utility for validating inputs to file writing functions
validate.file.write.inputs <- function(data, output.dir, file.prefix) {
    # check that input is a data.frame
    if (!is.data.frame(data)) {
        stop('input data must be a data.frame');
        }

    # check that file.prefix is a character string
    if (!is.null(file.prefix) && !is.character(file.prefix)) {
        stop('file.prefix must be a character string');
        }

    # check that output.dir is a directory
    if (!file.exists(output.dir)) {
        stop(paste('output path', output.dir, 'not found'));
        }
    }

# utility for writing pgs summary data to a file
write.per.sample.pgs.table <- function(per.sample.pgs.summary.data, output.dir, file.prefix = NULL) {
    # validate inputs
    validate.file.write.inputs(per.sample.pgs.summary.data, output.dir, file.prefix);

    # generate filename
    filename.for.per.sample.pgs.summary.data <- generate.filename(
        project.stem = file.prefix,
        file.core = 'per-sample-pgs-summary',
        extension = 'txt'
        );

    # write PGS per sample data to file
    write.table(
        x = per.sample.pgs.summary.data,
        file = file.path(output.dir, filename.for.per.sample.pgs.summary.data),
        sep = '\t',
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
        );
    }

# utility for writing pgs regression output to a file
write.pgs.regression.table <- function(pgs.regression.output, output.dir, file.prefix = NULL) {
    # validate inputs
    validate.file.write.inputs(pgs.regression.output, output.dir, file.prefix);

    # generate filename
    filename.for.pgs.regression.output <- generate.filename(
        project.stem = file.prefix,
        file.core = 'pgs-regression-output',
        extension = 'txt'
        );

    # write PGS regression output to file
    write.table(
        x = pgs.regression.output,
        file = file.path(output.dir, filename.for.pgs.regression.output),
        sep = '\t',
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
        );
    }

# Utility function for writing apply.polygenic.score outputs to file
#' @title Write apply.polygenic.score output to file
#' @description A utility function that writes the two data frames outputted by apply.polygenic.score to two tab-delimited text files.
#' @param apply.polygenic.score.output list of two data frames: pgs.output and regression.output
#' @param output.dir character string of the path to write both output files
#' @param file.prefix character string of the file prefix to use for both output files
#' @export
write.apply.polygenic.score.output.to.file <- function(apply.polygenic.score.output, output.dir, file.prefix = NULL) {
    # validate data input
    if (!is.list(apply.polygenic.score.output)) {
        stop('apply.polygenic.score.output must be a list');
        }
    if (!all(c('pgs.output', 'regression.output') %in% names(apply.polygenic.score.output))) {
        stop('apply.polygenic.score.output must contain pgs.output and regression.output');
        }

    write.per.sample.pgs.table(
        per.sample.pgs.summary.data = apply.polygenic.score.output$pgs.output,
        output.dir = output.dir,
        file.prefix = file.prefix
        );

    if (!is.null(apply.polygenic.score.output$regression.output)) {
        write.pgs.regression.table(
            pgs.regression.output = apply.polygenic.score.output$regression.output,
            output.dir = output.dir,
            file.prefix = file.prefix
            );
        }

    }
