rank.plotting.input.checks <- function(pgs.data, phenotype.columns) {
    # check pgs.data
    if (!is.data.frame(pgs.data)) {
        stop('pgs.data must be a data frame');
        }

    # validate phenotype.columns
    if (is.null(phenotype.columns)) {
        stop('phenotype.columns must be specified');
        }
    if (!is.character(phenotype.columns)) {
        stop('phenotype.columns must be a character vector');
        }
    if (!all(phenotype.columns %in% colnames(pgs.data))) {
        stop('phenotype.columns must be a subset of the column names in pgs.data');
        }

    # check for required columns
    required.columns <- c('Indiv', 'percentile', 'decile', 'quartile', 'n.missing.genotypes');
    if (!all(required.columns %in% colnames(pgs.data))) {
        stop('pgs.data must contain columns for Indiv, percentile, decile, and quartile');
        }

    }

plot.pgs.rank <- function(
    pgs.data,
    phenotype.columns,
    output.dir = getwd(),
    filename.prefix = NULL,
    file.extension = 'png',
    width = 8,
    height = 8,
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    titles.cex = 1.5
    ) {

    # check input
    rank.plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns);

    # factor Indiv by perentile rank
    pgs.data$Indiv <- factor(pgs.data$Indiv, levels = pgs.data$Indiv[rev(order(pgs.data$percentile))]);

    phenotype.data <- subset(pgs.data, select = phenotype.columns);
    # identify phenotype variable type
    phenotype.index.by.type <- classify.variable.type(data = phenotype.data);

    # Plotting
    # percentile rank barplot
    rank.barplot <- BoutrosLab.plotting.general::create.barplot(
        formula = percentile ~ Indiv,
        data = pgs.data,
        ylimits = c(0, 1.05),
        yat = seq(0, 1, 0.2),
        xlab.label = '',
        ylab.label = 'Percentile',
        main = ''
        # main.cex = titles.cex,
        # ylab.cex = titles.cex,
        # xaxis.cex = xaxis.cex,
        # yaxis.cex = yaxis.cex
        );

    # missing genotypes barplot
    missing.genotype.count.ymax <- ifelse(max(pgs.data$n.missing.genotypes) == 0, 1, max(pgs.data$n.missing.genotypes));
    missing.genotypes.barplot <- BoutrosLab.plotting.general::create.barplot(
        formula = n.missing.genotypes ~ Indiv,
        data = pgs.data,
        ylimits = c(0, missing.genotype.count.ymax * 1.05),
        yat = seq(0, missing.genotype.count.ymax, by = 1),
        xlab.label = '',
        ylab.label = 'Missing Genotypes',
        main = ''
        # main.cex = titles.cex,
        # ylab.cex = titles.cex,
        # xaxis.cex = xaxis.cex,
        # yaxis.cex = yaxis.cex
        );

    # percentile covariate bars
    



    }