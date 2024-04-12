# utility function for checking PGS plotting inputs
pgs.distribution.plotting.input.checks <- function(pgs.data, phenotype.columns, output.dir) {
    # Check that pgs.data is a data.frame
    if (!is.data.frame(pgs.data)) {
        stop('pgs.data must be a data.frame');
        }

    # validate phenotype.columns
    if (!is.null(phenotype.columns)) {
        if (!is.character(phenotype.columns)) {
            stop('phenotype.columns must be a character vector');
            }
        if (!all(phenotype.columns %in% colnames(pgs.data))) {
            stop('phenotype.columns must be a subset of the column names in pgs.data');
            }
        # check for collisions with recognized PGS columns
        recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
        if (any(phenotype.columns %in% recognized.pgs.colnames)) {
            stop('phenotype.columns cannot contain recognized PGS column names');
            }
        }

    # Identify possible PGS columns
    recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
    if (!any(recognized.pgs.colnames %in% colnames(pgs.data))) {
        stop('No recognized PGS columns found in pgs.data');
        }

    # validate output.dir
    if (!is.null(output.dir) && !dir.exists(output.dir)) {
        stop(paste0(output.dir), ' does not exist');
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
#' @description Plot density curves of PGS data from automatically detected PGS columns outputted by apply.polygenic.score (PGS, PGS.with.replaced.missing, PGS.with.normalized.missing).
#' If phenotype columns are provided, multiple density curves are plotted for automatically detected categories for each categorical variable.
#' @param pgs.data data.frame with PGS data
#' @param phenotype.columns character vector of phenotype columns in pgs.data to plot (optional)
#' @param output.dir character directory to save output plots
#' @param filename.prefix character prefix for output filenames
#' @param file.extension character file extension for output plots
#' @param tidy.titles logical whether to reformat PGS plot titles to remove periods
#' @param width numeric width of output plot in inches
#' @param height numeric height of output plot in inches
#' @param xaxes.cex numeric cex for all x-axis labels in inches
#' @param yaxes.cex numeric cex for all y-axis labels in inches
#' @param titles.cex numeric cex for all plot titles in inches
#' @param border.padding numeric padding for plot borders
#' @return if no output.dir is provided, a multipanel lattice plot object is returned, otherwise a plot is written to the indicated path and NULL is returned.
#' @export
plot.pgs.density <- function(
    pgs.data,
    phenotype.columns = NULL,
    output.dir = NULL,
    filename.prefix = NULL,
    file.extension = 'png',
    tidy.titles = FALSE,
    width = 10,
    height = 10,
    xaxes.cex = 1.5,
    yaxes.cex = 1.5,
    titles.cex = 1.5,
    border.padding = 1
    ) {
    # check input
    pgs.distribution.plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, output.dir = output.dir);

    # identify PGS columns
    recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
    pgs.columns <- colnames(pgs.data)[colnames(pgs.data) %in% recognized.pgs.colnames];

    # identify categorical phenotype variables for plotting
    if (!is.null(phenotype.columns)) {
        phenotype.data <- subset(pgs.data, select = phenotype.columns);
        phenotype.index.by.type <- classify.variable.type(data = phenotype.data);
        phenotype.data.for.plotting <- subset(phenotype.data, select = phenotype.index.by.type$binary | phenotype.index.by.type$other);
        }

    # Plotting
    pgs.density.plots <- list();
    pgs.density.by.phenotype.plots <- list();
    # iterate over PGS inputs
    for (pgs.column in pgs.columns) {
        ### Single Density Plots ###
        pgs.data.for.plotting <- data.frame(pgs.data[ , pgs.column]);

        # tidy titles
        if (tidy.titles) {
            pgs.column.main <- gsub(pattern = '\\.', replacement = ' ', x = pgs.column);
            } else {
                pgs.column.main <- pgs.column;
            }

        pgs.density.plots[[pgs.column]] <- BoutrosLab.plotting.general::create.densityplot(
            x = pgs.data.for.plotting,
            ylab.label = NULL,
            main = pgs.column.main,
            main.cex = titles.cex,
            yaxis.cex = yaxes.cex,
            xaxis.cex = xaxes.cex,
            lwd = 2
            );

        ### Density Plots by Phenotype ###
        if (!is.null(phenotype.columns)) {
            pgs.by.phenotype <- split.pgs.by.phenotype(pgs = pgs.data[ , pgs.column], phenotype.data = phenotype.data.for.plotting);

            # remove phenotype categories containing fewer than 2 samples
            pgs.by.phenotype <- lapply(
                X = pgs.by.phenotype,
                FUN = function(x) {
                    x <- x[sapply(x, length) > 1];
                    }
                );

            # iterate over phenotype variables
            for (phenotype in names(pgs.by.phenotype)) {
                pgs.data.for.plotting <- pgs.by.phenotype[[phenotype]];

                pgs.density.by.phenotype.plots[[paste0(pgs.column,'_',phenotype)]] <- BoutrosLab.plotting.general::create.densityplot(
                    x = pgs.data.for.plotting,
                    xlab.label = phenotype,
                    ylab.label = NULL,
                    main = '',
                    main.cex = 0,
                    xlab.cex = titles.cex,
                    yaxis.cex = yaxes.cex,
                    xaxis.cex = xaxes.cex,
                    lwd = 2,
                    col = BoutrosLab.plotting.general::default.colours(length(pgs.data.for.plotting)),
                    # Legend
                    legend = list(
                        inside = list(
                            fun = lattice::draw.key,
                            args = list(
                                key = list(
                                    lines = list(
                                        col = BoutrosLab.plotting.general::default.colours(length(pgs.data.for.plotting)),
                                        lwd = 2
                                        ),
                                    text = list(
                                        lab = names(pgs.data.for.plotting)
                                        ),
                                    padding.text = c(0,5,0),
                                    cex = 1
                                    )
                                ),
                            x = 0.65,
                            y = 0.97,
                            draw = FALSE
                            )
                        )
                    );
                }

            }

        }

    # organize filename if plot writing requested
    if (!is.null(output.dir)) {

        if (is.null(filename.prefix)) {
            filename.prefix <- 'ApplyPolygenicScore-Plot';
            }
        # construct multipanel plot
        filename.for.density.multiplot <- generate.filename(
            project.stem = filename.prefix,
            file.core = 'pgs-density',
            extension = file.extension
            );

        output.path <- file.path(output.dir, filename.for.density.multiplot);
        } else {
            output.path <- NULL;
        }

    # assemble multipanel plot
    if (length(pgs.density.by.phenotype.plots) != 0) {
        density.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
            plot.objects = c(pgs.density.plots, pgs.density.by.phenotype.plots),
            filename = output.path,
            layout.height = 1 + length(pgs.by.phenotype),
            layout.width = length(pgs.density.plots),
            main = '',
            main.cex = 0,
            ylab.label = 'Density',
            ylab.cex = titles.cex,
            ylab.axis.padding = -5,
            width = width,
            height = height,
            x.spacing = 1.5,
            y.spacing = 0,
            left.padding = border.padding,
            right.padding = border.padding,
            bottom.padding = border.padding,
            top.padding = border.padding
            );
        } else {
            density.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
                plot.objects = pgs.density.plots,
                filename = output.path,
                layout.height = 1,
                layout.width = length(pgs.density.plots),
                main = '',
                main.cex = 0,
                ylab.label = 'Density',
                ylab.cex = titles.cex,
                ylab.axis.padding = -5,
                width = width,
                height = height,
                x.spacing = 1.5,
                y.spacing = 0,
                left.padding = border.padding,
                right.padding = border.padding,
                bottom.padding = border.padding,
                top.padding = border.padding
                );

        }
    return(density.multipanel); # this returns null when filename is provided to create.multipanelplot

    }

#' @title Plot PGS Scatterplots
#' @description Create scatterplots for PGS data with continuous phenotype variables
#' @param pgs.data data.frame with PGS and phenotype data
#' @param phenotype.columns character vector of continuous phenotype columns in pgs.data to plot
#' @param output.dir character directory to save output plots
#' @param filename.prefix character prefix for output filenames
#' @param file.extension character file extension for output plots
#' @param compute.correlation logical whether to compute correlation between PGS and phenotype and display in plot
#' @param corr.legend.corner numeric vector indicating the corner of the correlation legend e.g. c(0,1) for top left
#' @param corr.legend.cex numeric cex for correlation legend
#' @param width numeric width of output plot in inches
#' @param height numeric height of output plot in inches
#' @param xaxes.cex numeric cex for x-axis labels
#' @param yaxes.cex numeric cex for y-axis labels
#' @param titles.cex numeric cex for plot titles
#' @param point.cex numeric cex for plot points
#' @return multipanel plot object
#' @export
plot.pgs.with.continuous.phenotype <- function(
    pgs.data,
    phenotype.columns,
    output.dir = NULL,
    filename.prefix = NULL,
    file.extension = 'png',
    tidy.titles = FALSE,
    compute.correlation = TRUE,
    corr.legend.corner = c(0,1),
    corr.legend.cex = 1.5,
    width = 10,
    height = 10,
    xaxes.cex = 1.5,
    yaxes.cex = 1.5,
    titles.cex = 1.5,
    point.cex = 0.75,
    border.padding = 1
    ) {
    # check input
    pgs.distribution.plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, output.dir = output.dir);

    recognized.pgs.colnames <- c('PGS', 'PGS.with.replaced.missing', 'PGS.with.normalized.missing');
    pgs.columns <- colnames(pgs.data)[colnames(pgs.data) %in% recognized.pgs.colnames];

    # identify continuous phenotype variables for plotting
    phenotype.data <- subset(pgs.data, select = phenotype.columns);
    phenotype.index.by.type <- classify.variable.type(data = phenotype.data);
    phenotype.data.for.plotting <- subset(phenotype.data, select = phenotype.index.by.type$continuous);

    if (length(phenotype.data.for.plotting) == 0) {
        stop('No continuous phenotype variables detected');
        }

    # Plotting
    pgs.scatterplots <- list();
    for (phenotype in colnames(phenotype.data.for.plotting)) {

        for (pgs.column in pgs.columns) {

            # handle tidy titles
            if (tidy.titles) {
                pgs.column.label <- gsub(pattern = '\\.', replacement = ' ', x = pgs.column);
                } else {
                    pgs.column.label <- pgs.column;
                }

            # handle correlation key
            if (compute.correlation) {
                correlation.legend <- list(
                    inside = list(
                        fun = lattice::draw.key,
                        args = list(
                            key = get.corr.key(
                                x = pgs.data[ , pgs.column],
                                y = phenotype.data.for.plotting[ , phenotype],
                                label.items = c('spearman', 'spearman.p'),
                                alpha.background = 0,
                                key.cex = corr.legend.cex
                                )
                            ),
                    #x = 0.04,
                    #y = 0.97,
                    corner = corr.legend.corner
                    )
                );
                } else {
                    correlation.legend <- NULL;
                }

            pgs.scatterplots[[paste0(pgs.column,'_',phenotype)]] <- BoutrosLab.plotting.general::create.scatterplot(
                formula = as.formula(paste0(phenotype, ' ~ ', pgs.column)),
                data = pgs.data,
                type = 'p',
                cex = point.cex,
                xlab.label = pgs.column.label,
                ylab.label = phenotype,
                main = '',
                main.cex = 0,
                # Correlation Legend
                legend = correlation.legend,
                ylab.cex = titles.cex,
                xlab.cex = titles.cex,
                yaxis.cex = yaxes.cex,
                xaxis.cex = xaxes.cex
                );
            }

        }

    # organize filename if plot writing requested
    if (!is.null(output.dir)) {

        if (is.null(filename.prefix)) {
            filename.prefix <- 'ApplyPolygenicScore-Plot';
            }
        # construct multipanel plot
        filename.for.scatterplot.multiplot <- generate.filename(
            project.stem = filename.prefix,
            file.core = 'pgs-scatter',
            extension = file.extension
            );

        output.path <- file.path(output.dir, filename.for.scatterplot.multiplot);
        } else {
            output.path <- NULL;
        }

    # assemble multipanel plot
    scatterplot.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
        plot.objects = pgs.scatterplots,
        filename = output.path,
        layout.height = length(phenotype.data.for.plotting),
        layout.width = length(pgs.columns),
        main = '',
        main.cex = 0,
        width = width,
        height = height,
        x.spacing = 1.5,
        y.spacing = 0,
        left.padding = border.padding,
        right.padding = border.padding,
        bottom.padding = border.padding,
        top.padding = border.padding
        );

    return(scatterplot.multipanel);
    }
