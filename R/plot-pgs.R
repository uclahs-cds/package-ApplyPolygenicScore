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

plot.pgs.density <- function(pgs.data, phenotype.columns, filename.prefix = NULL) {
    # check input
    plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, filname.prefix = filename.prefix);    

    pgs.columns <- colnames(pgs.data)[colnames(pgs.data) %in% recognized.pgs.colnames];

    # Plot density of each PGS as separate plot
    global.cex <- 1.5;
    pgs.density.plots <- list();
    for (pgs.column in pgs.columns) {
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

        }

}