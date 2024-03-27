# utility function for checking PGS plotting inputs
plotting.input.checks <- function(pgs.data, phenotype.columns, filname.prefix = NULL) {
    # Check that pgs.data is a data.frame
    if (!is.data.frame(pgs.data)) {
        stop("pgs.data must be a data.frame");
        }

    # Check that phenotype.columns is a character vector
    if (!is.null(phenotype.columns) && !is.character(phenotype.columns)) {
        stop("phenotype.columns must be a character vector");
        }

    # Check that phenotype.columns are in pgs.data
    if (!all(phenotype.columns %in% colnames(pgs.data))) {
        stop("phenotype.columns must be in pgs.data");
        }

    # Identify possible PGS columns
    recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
    if (!any(recognized.pgs.colnames %in% colnames(pgs.data))) {
        stop("No recognized PGS columns found in pgs.data");
        }
    }

# utility function that formats PGS data for plotting in density plots
split.pgs.by.phenotype <- function(pgs, phenotype.data) {
    output <- lapply(
        X = phenotype.data,
        FUN = function(phenotype.column) {
            split.data <- aggregate(
                x = pgs,
                by = list(phenotype.column),
                FUN = c
                );
            split.data.list <- split.data$x;
            names(split.data.list) <- split.data$Group.1;
            return(split.data.list);
            }
        );
    return(output);
    }

#' @title Plot PGS Density
#' @description Create density plots for PGS data
#' @param pgs.data data.frame with PGS data
#' @param phenotype.columns character vector of phenotype columns in pgs.data to plot (optional)
#' @param output.dir character directory to save output plots
#' @param filename.prefix character prefix for output filenames
#' @param file.extension character file extension for output plots
#' @param width numeric width of output plot
#' @param height numeric height of output plot
#' @return multipanel plot object
plot.pgs.density <- function(
    pgs.data,
    phenotype.columns,
    output.dir = getwd(),
    filename.prefix = NULL,
    file.extension = 'png',
    width = 10,
    height = 10
    ) {
    # check input
    plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, filname.prefix = filename.prefix);    

    recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
    pgs.columns <- colnames(pgs.data)[colnames(pgs.data) %in% recognized.pgs.colnames];

    # identify categorical phenotype variables for plotting
    if (!is.null(phenotype.columns)) {
        phenotype.data <- pgs.data[ , phenotype.columns];
        phenotype.index.by.type <- classify.variable.type(data = phenotype.data);
        phenotype.data.for.plotting <- subset(phenotype.data, select = phenotype.index.by.type$binary | phenotype.index.by.type$other);
        }

    # Plotting
    global.cex <- 1.5;
    pgs.density.plots <- list();
    pgs.density.by.phenotype.plots <- list();
    for (pgs.column in pgs.columns) {
        ### Single Density Plots ###
        pgs.data.for.plotting <- data.frame(pgs.data[ , pgs.column]);

        pgs.density.plots[[pgs.column]] <- BoutrosLab.plotting.general::create.densityplot(
            x = pgs.data.for.plotting,
            ylab.label = 'Density',
            main = pgs.column,
            main.cex = global.cex,
            ylab.cex = global.cex,
            yaxis.cex = global.cex,
            xaxis.cex = global.cex,
            lwd = 2
            );

        ### Density Plots by Phenotype ###
        if (!is.null(phenotype.columns)) {
            pgs.by.phenotype <- split.pgs.by.phenotype(pgs = pgs.data[ , pgs.column], phenotype.data = phenotype.data.for.plotting);
            for (phenotype in names(pgs.by.phenotype)) {
                pgs.data.for.plotting <- pgs.by.phenotype[[phenotype]];

                pgs.density.by.phenotype.plots[[paste0(pgs.column,'_',phenotype)]] <- BoutrosLab.plotting.general::create.densityplot(
                    x = pgs.data.for.plotting,
                    xlab.label = phenotype,
                    ylab.label = 'Density',
                    main = '',
                    main.cex = 0,
                    ylab.cex = global.cex,
                    xlab.cex = global.cex,
                    yaxis.cex = global.cex,
                    xaxis.cex = global.cex,
                    lwd = 2,
                    col = BoutrosLab.plotting.general::default.colours(length(pgs.data.for.plotting))
                    );
                }

            }
        
        }

    # create filename
    if (is.null(filename.prefix)) {
        filename.prefix <- 'ApplyPolygenicScore-Plot';
        }
    filename.for.density.multiplot <- generate.filename(
        project.stem = filename.prefix,
        file.core = 'pgs-density',
        extension = file.extension
        );

    # assemble multipanel plot
    if (length(pgs.density.by.phenotype.plots) != 0) {
        density.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
            plot.objects = c(pgs.density.plots, pgs.density.by.phenotype.plots),
            filename = file.path(output.dir, filename.for.density.multiplot),
            layout.height = 1 + length(pgs.density.by.phenotype.plots),
            layout.width = length(pgs.density.plots),
            main = '',
            main.cex = 0,
            width = width,
            height = height
            );
        } else {
            density.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
                plot.objects = pgs.density.plots,
                filename = file.path(output.dir, filename.for.density.multiplot),
                layout.height = 1,
                layout.width = length(pgs.density.plots),
                main = '',
                main.cex = 0,
                width = width,
                height = height
                );

        }
    return(density.multipanel); # this returns null when filename is provided to create.multipanelplot

    }
