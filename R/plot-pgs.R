# utility function for checking PGS plotting inputs
plotting.input.checks <- function(pgs.data, phenotype.columns, filname.prefix = NULL) {
    # Check that pgs.data is a data.frame
    if (!is.data.frame(pgs.data)) {
        stop("pgs.data must be a data.frame");
        }

    # Check that phenotype.columns is a character vector
    if (!is.character(phenotype.columns)) {
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

plot.pgs.density <- function(pgs.data, phenotype.columns, filename.prefix = NULL) {
    # check input
    plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, filname.prefix = filename.prefix);    

    pgs.columns <- colnames(pgs.data)[colnames(pgs.data) %in% recognized.pgs.colnames];

    # identify categorical phenotype variables for plotting
    phenotype.data <- pgs.data[ , phenotype.columns];
    phenotype.index.by.type <- classify.variable.type(data = phenotype.data);
    phenotype.data.for.plotting <- data.frame(phenotype.data[ , phenotype.index.by.type$binary | phenotype.index.by.type$other]);

    # Plot density of each PGS as separate plot
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
        pgs.by.phenotype <- split.pgs.by.phenotype(pgs = pgs.data[ , pgs.column], phenotype.data = phenotype.data.for.plotting);
        for (phenotype in names(pgs.by.phenotype)) {
            pgs.data.for.plotting <- pgs.by.phenotype[[phenotype]];

            pgs.density.by.phenotype.plots[[pgs.column]][[phenotype]] <- BoutrosLab.plotting.general::create.densityplot(
                x = pgs.data.for.plotting,
                xlab.label = pgs.column,
                ylab.label = 'Density',
                main = phenotype,
                main.cex = global.cex,
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