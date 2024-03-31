# utility function for color handling
# Given a vector of feature names and given a named color scheme
# (named vector of colors with names being feature names)
# return a vector of colors corresponding to each feature.
create.feature.color.vector <- function(features, color.scheme) {
  # Replace features with corresponding color
  feature.colors <- sapply(
    X = features,
    FUN = function(feature) {
        color.scheme[feature]
        }
    );
  names(feature.colors) <- features;
  return(feature.colors);
  }

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
        xlab.label = '',
        ylab.label = 'Missing Genotypes',
        main = ''
        # main.cex = titles.cex,
        # ylab.cex = titles.cex,
        # xaxis.cex = xaxis.cex,
        # yaxis.cex = yaxis.cex
        );

    # percentile covariate bars

    # quantile color schemes
    decile.color.scheme <- c(
        paste0('grey', seq(1,100, length.out = 10))
        );
    names(decile.color.scheme) <- as.character(1:10);

    quartile.color.scheme <- c(
        decile.color.scheme[seq(1, 10, length.out = 4)]
        );
    names(quartile.color.scheme) <- as.character(1:4);

    if ('n.percentile' %in% colnames(pgs.data)) {
        percentile.color.scheme <- c(
            paste0('grey', round(seq(1,100, length.out = length(unique(pgs.data$n.percentile)))))
            );
        names(percentile.color.scheme) <- as.character(unique(pgs.data$n.percentile));
        percentile.color.scheme.list <- list(decile = decile.color.scheme, quartile = quartile.color.scheme, n.percentile = percentile.color.scheme);
        percentile.covariate.data <- subset(pgs.data, select = c('decile', 'quartile', 'n.percentile'));
        } else {
            percentile.covariate.data <- subset(pgs.data, select = c('decile', 'quartile'));
            percentile.color.scheme.list <- list(decile = decile.color.scheme, quartile = quartile.color.scheme);
            }

    percentile.covariate.df <- sapply(
        X = colnames(percentile.covariate.data),#1:length(percentile.covariate.data),
        FUN = function(x) {
            covariate.values <- percentile.covariate.data[ , x];
            color.vector <- create.feature.color.vector(
                features = covariate.values,
                color.scheme = percentile.color.scheme.list[[x]]
                );
            return(color.vector)
            }
        );
    percentile.covariate.df <- data.frame(t(percentile.covariate.df));
    colnames(percentile.covariate.df) <- pgs.data$Indiv;
    # order by percentile
    percentile.covariate.df <- percentile.covariate.df[ , order(pgs.data$percentile)];

    }
