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
                features = as.character(covariate.values),
                color.scheme = percentile.color.scheme.list[[x]]
                );
            return(color.vector)
            }
        );
    percentile.covariate.df <- data.frame(t(percentile.covariate.df));
    colnames(percentile.covariate.df) <- pgs.data$Indiv;
    # order by percentile
    percentile.covariate.df <- percentile.covariate.df[ , rev(order(pgs.data$percentile))];

    # percentile covariate heatmap
    percentile.covariate.heatmap <- BoutrosLab.plotting.general::create.heatmap(
        x = percentile.covariate.df,
        input.colours = TRUE,
        clustering.method = 'none',
        same.as.matrix = TRUE,
        print.colour.key = FALSE,
        yaxis.lab = rownames(percentile.covariate.df),
        ylab.cex = 1
        # main.cex = titles.cex,
        # ylab.cex = titles.cex,
        # xaxis.cex = xaxis.cex,
        # yaxis.cex = yaxis.cex
        );

    max.binary.colors <- sum(phenotype.index.by.type$binary | phenotype.index.by.type$continuous);
    if (max.binary.colors > 0) {
        binary.color.schemes <- default.colours(
            number.of.colours = rep(2, max.binary.colors + 1),
            palette = rep('binary', max.binary.colors + 1)
            );
        # remove black and white from binary color schemes
        binary.color.schemes[[1]] <- NULL;
        binary.color.schemes.start.index <- 1;
        }

    categorical.phenotype.heatmap <- NULL;
    continuous.phenotype.heatmap <- NULL;

    if (any(phenotype.index.by.type$binary) | any(phenotype.index.by.type$other)) {
        binary.phenotype.df <- NULL;
        other.phenotype.df <- NULL;

        if (any(phenotype.index.by.type$binary)) {
            binary.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$binary);

            for (i in 1:length(binary.color.schemes)) {
                names(binary.color.schemes[[i]]) <- as.character(sort(unique(binary.phenotype.data[ , i])));
                }

            names(binary.color.schemes) <- colnames(binary.phenotype.data);

            binary.phenotype.df <- sapply(
                X = colnames(binary.phenotype.data),
                FUN = function(x) {
                    phenotype.values <- binary.phenotype.data[ , x];
                    color.vector <- create.feature.color.vector(
                        features = as.character(phenotype.values),
                        color.scheme = binary.color.schemes[[x]]
                        );
                    return(color.vector)
                    }
                );
            binary.phenotype.df <- data.frame(t(binary.phenotype.df));
            colnames(binary.phenotype.df) <- pgs.data$Indiv;
            # order by percentile
            binary.phenotype.df <- binary.phenotype.df[ , rev(order(pgs.data$percentile))];

            binary.color.schemes.start.index <- binary.color.schemes.start.index + ncol(binary.phenotype.data);

            # # binary phenotype covariate heatmap
            # binary.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
            #     x = data.frame(binary.phenotype.df),
            #     input.colours = TRUE,
            #     clustering.method = 'none',
            #     same.as.matrix = TRUE,
            #     print.colour.key = FALSE,
            #     yaxis.lab = rownames(binary.phenotype.df),
            #     ylab.cex = 1
            #     # main.cex = titles.cex,
            #     # ylab.cex = titles.cex,
            #     # xaxis.cex = xaxis.cex,
            #     # yaxis.cex = yaxis.cex
            #     );

            }

        if (any(phenotype.index.by.type$other)) {
            other.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$other);
            number.of.categories <- lapply(
                X = other.phenotype.data,
                FUN = function(x) {
                    length(unique(x))
                    }
                );
            max.colors <- 12;
            all.qual.colors <- default.colours(number.of.colors <- max.colors, palette = 'qual');
            other.color.schemes <- list();
            start.palette <- 1;
            for (i in 1:length(number.of.categories)) {
                other.color.schemes[[i]] <- all.qual.colors[start.palette:(number.of.categories[[i]] + start.palette - 1)];
                start.palette <- start.palette + number.of.categories[[i]];
                }

            for (i in 1:length(other.color.schemes)) {
                names(other.color.schemes[[i]]) <- as.character(sort(unique(other.phenotype.data[ , i])));
                }

            names(other.color.schemes) <- colnames(other.phenotype.data);

            other.phenotype.df <- sapply(
                X = colnames(other.phenotype.data),
                FUN = function(x) {
                    phenotype.values <- other.phenotype.data[ , x];
                    color.vector <- create.feature.color.vector(
                        features = as.character(phenotype.values),
                        color.scheme = other.color.schemes[[x]]
                        );
                    return(color.vector)
                    }
                );
            other.phenotype.df <- data.frame(t(other.phenotype.df));
            colnames(other.phenotype.df) <- pgs.data$Indiv;

            # order by percentile
            other.phenotype.df <- other.phenotype.df[ , rev(order(pgs.data$percentile))];

            # # other phenotype covariate heatmap
            # other.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
            #     x = data.frame(other.phenotype.df),
            #     input.colours = TRUE,
            #     clustering.method = 'none',
            #     same.as.matrix = TRUE,
            #     print.colour.key = FALSE,
            #     yaxis.lab = rownames(other.phenotype.df),
            #     ylab.cex = 1
            #     # main.cex = titles.cex,
            #     # ylab.cex = titles.cex,
            #     # xaxis.cex = xaxis.cex,
            #     # yaxis.cex = yaxis.cex
            #    );
            }

        # combine binary and categorical phenotype covariates into one heatmap
        all.category.phenotype.df <- rbind(binary.phenotype.df, other.phenotype.df);

        # plot binary and categorical phenotype covariate heatmap
        categorical.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
            x = all.category.phenotype.df,
            input.colours = TRUE,
            clustering.method = 'none',
            same.as.matrix = TRUE,
            print.colour.key = FALSE,
            yaxis.lab = rownames(all.category.phenotype.df),
            ylab.cex = 1
            # main.cex = titles.cex,
            # ylab.cex = titles.cex,
            # xaxis.cex = xaxis.cex,
            # yaxis.cex = yaxis.cex
            );

        }

    # continuous phenotype covariate heatmap
    if (any(phenotype.index.by.type$continuous)) {
        continuous.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$continuous);

        continuous.color.schemes <- binary.color.schemes[binary.color.schemes.start.index:length(binary.color.schemes)]

        continuous.phenotypes.df <- sapply(
            X = 1:length(continuous.phenotype.data),
            FUN = function(x) {
                phenotype.values <- continuous.phenotype.data[ , x];
                color.scheme <- colorRampPalette(continuous.color.schemes[[x]])(length(unique(phenotype.values)));
                names(color.scheme) <- as.character(sort(unique(phenotype.values)));
                color.vector <- create.feature.color.vector(
                    features = as.character(phenotype.values),
                    color.scheme = color.scheme
                    );
                }
            );
        continuous.phenotypes.df <- data.frame(t(continuous.phenotypes.df));
        colnames(continuous.phenotypes.df) <- pgs.data$Indiv;
        rownames(continuous.phenotypes.df) <- colnames(continuous.phenotype.data);
        # order by percentile
        continuous.phenotypes.df <- continuous.phenotypes.df[ , rev(order(pgs.data$percentile))];

        # continuous phenotype covariate heatmap
        continuous.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
            x = continuous.phenotypes.df,
            input.colours = TRUE,
            clustering.method = 'none',
            same.as.matrix = TRUE,
            print.colour.key = FALSE,
            yaxis.lab = rownames(continuous.phenotypes.df),
            ylab.cex = 1
            # main.cex = titles.cex,
            # ylab.cex = titles.cex,
            # xaxis.cex = xaxis.cex,
            # yaxis.cex = yaxis.cex
            );

        }

    if (is.null(filename.prefix)) {
        filename.prefix <- 'ApplyPolygenicScore-Plot';
        }
    # construct multipanel plot
    filename.for.rank.multiplot <- generate.filename(
        project.stem = filename.prefix,
        file.core = 'pgs-rank-plot',
        extension = file.extension
        );

    plot.list <- list(
        missing.genotypes.barplot,
        rank.barplot,
        percentile.covariate.heatmap,
        categorical.phenotype.heatmap,
        continuous.phenotype.heatmap
        );

    plot.heights <- rep(1, length(plot.list));
    plot.heights[1] <- 2;
    plot.heights[2] <- 5;

    multipanel.plot <- BoutrosLab.plotting.general::create.multipanelplot(
        plot.objects = plot.list,
        filename = file.path(output.dir, filename.for.rank.multiplot),
        main = '',
        main.cex = 0,
        layout.height = length(plot.list),
        layout.width = 1,
        plot.objects.heights = plot.heights,
        y.spacing = -1,
        width = width,
        height = height
        );

    }
