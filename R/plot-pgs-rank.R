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

# utility function for assembling a data frame of colors for heatmap plotting
# Given a data frame of covariate data and a list of color schemes with matching colnames/names
# return a data frame with all covariate data replaced by colors indicated by color schemes.
assemble.heatmap.colors <- function(covariate.data, color.scheme.list) {
    covariate.df <- sapply(
        X = colnames(covariate.data),
        FUN = function(x) {
            covariate.values <- covariate.data[ , x];
            color.vector <- create.feature.color.vector(
                features = as.character(covariate.values),
                color.scheme = color.scheme.list[[x]]
                );
            return(color.vector)
            }
        );
    }

# utility function for formatting coordinates of missing values in a heatmap matrix
get.na.coordinates.for.heatmap <- function(data) {
    na.boolean <- is.na(data);
    rev.na.boolean <- na.boolean[nrow(na.boolean):1, , drop = FALSE];
    na.col.coordinates <- which(na.boolean, arr.ind = TRUE);
    na.row.coordinates <- which(rev.na.boolean, arr.ind = TRUE);
    return(list(col = na.col.coordinates[ ,'col'], row = na.row.coordinates[ ,'row']));
    }

rank.plotting.input.checks <- function(pgs.data, phenotype.columns, output.dir) {
    # check pgs.data
    if (!is.data.frame(pgs.data)) {
        stop('pgs.data must be a data frame');
        }

    # validate phenotype.columns
    if (!is.null(phenotype.columns)) {
        if (!is.character(phenotype.columns)) {
            stop('phenotype.columns must be a character vector');
            }
        if (!all(phenotype.columns %in% colnames(pgs.data))) {
            stop('phenotype.columns must be a subset of the column names in pgs.data');
            }
    }

    # check for required columns
    required.columns <- c('Indiv', 'percentile', 'decile', 'quartile', 'n.missing.genotypes');
    if (!all(required.columns %in% colnames(pgs.data))) {
        stop('pgs.data must contain columns for Indiv, percentile, decile, quartile, and n.missing.genotypes');
        }

    # validate output.dir
    if (!is.null(output.dir) && !dir.exists(output.dir)) {
        stop(paste0(output.dir), ' does not exist');
        }

    }

plot.pgs.rank <- function(
    pgs.data,
    phenotype.columns = NULL,
    output.dir = NULL,
    filename.prefix = NULL,
    file.extension = 'png',
    width = 8,
    height = 8,
    xaxis.cex = 1.2,
    yaxis.cex = 1,
    titles.cex = 1.2,
    border.padding = 1
    ) {

    FILL.COLOR <- 'grey';

    # check input
    rank.plotting.input.checks(pgs.data = pgs.data, phenotype.columns = phenotype.columns, output.dir = output.dir);

    # factor Indiv by perentile rank
    pgs.data$Indiv <- factor(pgs.data$Indiv, levels = pgs.data$Indiv[order(pgs.data$percentile)]);

    # Plot percentile rank barplot
    rank.barplot <- BoutrosLab.plotting.general::create.barplot(
        formula = percentile ~ Indiv,
        data = pgs.data,
        ylimits = c(0, 1.05),
        yat = seq(0, 1, 0.2),
        xaxis.rot = 90,
        xlab.label = '',
        ylab.label = 'PGS Percentile',
        main = '',
        main.cex = 0,
        ylab.cex = titles.cex,
        xaxis.cex = xaxis.cex,
        yaxis.cex = yaxis.cex
        );

    # Plot missing genotypes barplot
    # handle plot limits in case where there are no missing genotypes
    missing.genotype.count.ymax <- ifelse(max(pgs.data$n.missing.genotypes) == 0, 1, max(pgs.data$n.missing.genotypes));
    missing.genotypes.barplot <- BoutrosLab.plotting.general::create.barplot(
        formula = n.missing.genotypes ~ Indiv,
        data = pgs.data,
        ylimits = c(0, missing.genotype.count.ymax * 1.05),
        xlab.label = '',
        ylab.label = 'Missing GT',
        xaxis.lab = '',
        yat = 'auto',
        main = '',#'Missing Genotypes',
        main.cex = 0,
        ylab.cex = titles.cex,
        xaxis.cex = 0,
        yaxis.cex = yaxis.cex
        );

    ## Begin Percentile Covariate Heatmap Assembly ##
    # Assemble covariate heatmap for deciles, quartiles, and user-defined percentiles

    # assign percentiles to shades of grey
    decile.color.scheme <- c(
        paste0('grey', seq(100,1, length.out = 10))
        );
    names(decile.color.scheme) <- as.character(1:10);

    quartile.color.scheme <- c(
        decile.color.scheme[seq(1, 10, length.out = 4)]
        );
    names(quartile.color.scheme) <- as.character(1:4);

    user.defined.percentile.column.index <- grep('percentile\\.[0-9]', colnames(pgs.data));
    if (length(user.defined.percentile.column.index) == 1) {
        # assign user-defined percentiles to shades of grey
        percentile.color.scheme <- c(
            paste0('grey', round(seq(100,1, length.out = length(unique(na.omit(pgs.data[ ,user.defined.percentile.column.index]))))))
            );
        names(percentile.color.scheme) <- as.character(sort(unique(na.omit(pgs.data[ ,user.defined.percentile.column.index]))));
        # assemble all percentiles color schemes and data
        percentile.color.scheme.list <- list(decile.color.scheme, quartile.color.scheme, percentile.color.scheme);
        names(percentile.color.scheme.list) <- c('decile', 'quartile', colnames(pgs.data)[user.defined.percentile.column.index]);
        percentile.covariate.data <- subset(pgs.data, select = c('decile', 'quartile', colnames(pgs.data)[user.defined.percentile.column.index]));
        } else {
            # assemble all percentiles color schemes and data
            percentile.covariate.data <- subset(pgs.data, select = c('decile', 'quartile'));
            percentile.color.scheme.list <- list(decile = decile.color.scheme, quartile = quartile.color.scheme);
            }

    # replace covariate data with colors
    percentile.covariate.df <- assemble.heatmap.colors(
        covariate.data = percentile.covariate.data,
        color.scheme.list = percentile.color.scheme.list
        );

    # reorient data frame to match barplot layout (samples on x-axis)
    percentile.covariate.df <- data.frame(t(percentile.covariate.df));
    colnames(percentile.covariate.df) <- pgs.data$Indiv;
    # order by percentile
    percentile.covariate.df <- percentile.covariate.df[ , order(pgs.data$percentile)];

    # save NA coordinates for labeling
    percentile.cov.na.coords <- get.na.coordinates.for.heatmap(percentile.covariate.df);
    # replace NA values with a color
    percentile.covariate.df[is.na(percentile.covariate.df)] <- FILL.COLOR;

    # Plot percentile covariate heatmap
    percentile.covariate.heatmap <- BoutrosLab.plotting.general::create.heatmap(
        x = percentile.covariate.df,
        input.colours = TRUE,
        clustering.method = 'none',
        same.as.matrix = TRUE,
        print.colour.key = FALSE,
        # missing value handling
        cell.text = 'NA',
        col.pos = percentile.cov.na.coords$col,
        row.pos = percentile.cov.na.coords$row,
        text.col = 'white',
        text.cex = 0.5,
        yaxis.lab = rownames(percentile.covariate.df),
        #ylab.cex = title.cex,
        # main.cex = titles.cex,
        # ylab.cex = titles.cex,
        xaxis.cex = xaxis.cex,
        yaxis.cex = yaxis.cex
        );

        percentile.covariates.legend <- list(list(
            title = 'Percentiles',
            colours = c(decile.color.scheme[1], decile.color.scheme[10]),
            labels = c('lowest', 'highest'),
            continuous = TRUE
            ));
        names(percentile.covariates.legend) <- 'legend';

    ## End Percentile Covariate Heatmap Assembly ##

    ## Begin Phenotype Covariate Heatmap Assembly ##

    # initialize variables for phenotype covariates as NULL since they are optional
    categorical.phenotype.heatmap <- NULL;
    continuous.phenotype.heatmap <- NULL;
    binary.covariates.legend <- NULL;
    categorical.covariates.legend <- NULL;
    continuous.covariates.legend <- NULL;

    if (!is.null(phenotype.columns)) {
        phenotype.data <- subset(pgs.data, select = phenotype.columns);
        # identify phenotype variable type
        phenotype.index.by.type <- classify.variable.type(data = phenotype.data);

        # retrieve binary color schemes sufficient for plotting all binary and continuous phenotypes
        max.binary.colors <- sum(phenotype.index.by.type$binary | phenotype.index.by.type$continuous);
        if (max.binary.colors > 0) {
            binary.color.schemes <- BoutrosLab.plotting.general::default.colours(
                number.of.colours = rep(2, max.binary.colors + 1),
                palette = rep('binary', max.binary.colors + 1)
                );
            # remove black and white from binary color schemes (always returned first by default.colours())
            binary.color.schemes[[1]] <- NULL;
            binary.color.schemes.start.index <- 1;
            }


        # assemble binary and categorical phenotype covariates in one heatmap
        if (any(phenotype.index.by.type$binary) | any(phenotype.index.by.type$other)) {
            binary.phenotype.df <- NULL;
            other.phenotype.df <- NULL;

            # binary phenotype covariate data
            if (any(phenotype.index.by.type$binary)) {
                binary.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$binary);

                # extract the required number of color schemes for binary covariates from pre-generated list
                binary.covariate.color.schemes <- binary.color.schemes[1:ncol(binary.phenotype.data)];
                for (i in 1:length(binary.covariate.color.schemes)) {
                    names(binary.covariate.color.schemes[[i]]) <- as.character(sort(unique(na.omit(binary.phenotype.data[ , i]))));
                    }

                names(binary.covariate.color.schemes) <- colnames(binary.phenotype.data);

                # replace covariate data with colors
                binary.phenotype.df <- assemble.heatmap.colors(
                    covariate.data = binary.phenotype.data,
                    color.scheme.list = binary.covariate.color.schemes
                    );
                # reorient data frame to match barplot layout (samples on x-axis)
                binary.phenotype.df <- data.frame(t(binary.phenotype.df));
                colnames(binary.phenotype.df) <- pgs.data$Indiv;
                # order by percentile
                binary.phenotype.df <- binary.phenotype.df[ , order(pgs.data$percentile)];
                # update binary color schemes start index for continuous phenotypes
                binary.color.schemes.start.index <- binary.color.schemes.start.index + ncol(binary.phenotype.data);

                # build legend for binary covariates
                binary.covariates.legend <- lapply(
                    X = 1:length(binary.covariate.color.schemes),
                    FUN = function(x) {
                        list(
                            title = names(binary.covariate.color.schemes)[x],
                            colours = binary.covariate.color.schemes[[x]],
                            labels = names(binary.covariate.color.schemes[[x]])
                            );
                        }
                    );
                names(binary.covariates.legend) <- rep('legend', length(binary.covariates.legend));
                }

            # categorical phenotype covariate data
            if (any(phenotype.index.by.type$other)) {
                other.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$other);

                # count number of categories in each covariate
                number.of.categories <- lapply(
                    X = other.phenotype.data,
                    FUN = function(x) {
                        length(unique(na.omit(x)))
                        }
                    );
                total.categories <- sum(unlist(number.of.categories));

                # retrieve colors for qualitative data that can be divided up into color schemes
                # only 12 distinct colors are available in the default color palette
                max.colors <- 12;
                suppressWarnings( #suppress grey scale incompatibility warnings
                    all.qual.colors <- default.colours(number.of.colors <- max.colors, palette = 'qual')
                    );
                # if there are more categories than colors, extend the size of the color palette by repeating colors
                if (total.categories > max.colors) {
                    all.qual.colors <- rep(all.qual.colors, ceiling(total.categories / max.colors));
                    }

                # assemble a color scheme for each categorical variable
                other.color.schemes <- list();
                start.palette <- 1;
                for (i in 1:length(number.of.categories)) {
                    other.color.schemes[[i]] <- all.qual.colors[start.palette:(number.of.categories[[i]] + start.palette - 1)];
                    start.palette <- start.palette + number.of.categories[[i]];
                    }

                for (i in 1:length(other.color.schemes)) {
                    names(other.color.schemes[[i]]) <- as.character(sort(unique(na.omit(other.phenotype.data[ , i]))));
                    }

                names(other.color.schemes) <- colnames(other.phenotype.data);

                # replace covariate data with colors
                other.phenotype.df <- assemble.heatmap.colors(
                    covariate.data = other.phenotype.data,
                    color.scheme.list = other.color.schemes
                    );
                # reorient data frame to match barplot layout (samples on x-axis)
                other.phenotype.df <- data.frame(t(other.phenotype.df));
                colnames(other.phenotype.df) <- pgs.data$Indiv;

                # order by percentile
                other.phenotype.df <- other.phenotype.df[ , order(pgs.data$percentile)];

                # legend
                categorical.covariates.legend <- lapply(
                    X = 1:length(other.color.schemes),
                    FUN = function(x) {
                        list(
                            title = names(other.color.schemes)[x],
                            colours = other.color.schemes[[x]],
                            labels = names(other.color.schemes[[x]])
                            );
                        }
                    );
                names(categorical.covariates.legend) <- rep('legend', length(categorical.covariates.legend));

                }

            # combine binary and categorical phenotype covariates into one heatmap
            all.category.phenotype.df <- rbind(binary.phenotype.df, other.phenotype.df);

            # save NA coordinates for labeling
            cat.phen.na.coords <- get.na.coordinates.for.heatmap(all.category.phenotype.df);
            # replace NA values with a color
            all.category.phenotype.df[is.na(all.category.phenotype.df)] <- FILL.COLOR;

            # plot binary and categorical phenotype covariate heatmap
            categorical.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
                x = all.category.phenotype.df,
                input.colours = TRUE,
                clustering.method = 'none',
                same.as.matrix = TRUE,
                print.colour.key = FALSE,
                cell.text = 'NA',
                col.pos = cat.phen.na.coords$col,
                row.pos = cat.phen.na.coords$row,
                text.col = 'white',
                text.cex = 0.5,
                yaxis.lab = NULL,
                ylab.cex = 0
                # main.cex = titles.cex,
                # ylab.cex = titles.cex,
                # xaxis.cex = xaxis.cex,
                # yaxis.cex = yaxis.cex
                );

            }

        # assemble continuous phenotype covariates in one heatmap
        if (any(phenotype.index.by.type$continuous)) {
            continuous.phenotype.data <- subset(phenotype.data, select = phenotype.index.by.type$continuous);

            # retreive remaining binary color schemes
            continuous.color.schemes <- binary.color.schemes[binary.color.schemes.start.index:length(binary.color.schemes)]
            names(continuous.color.schemes) <- colnames(continuous.phenotype.data);

            # assemble color ramps for continuous covariates and replace covariate data with colors
            continuous.phenotypes.df <- sapply(
                X = 1:length(continuous.phenotype.data),
                FUN = function(x) {
                    phenotype.values <- continuous.phenotype.data[ , x];
                    color.scheme <- colorRampPalette(continuous.color.schemes[[x]])(length(unique(na.omit(phenotype.values))));
                    names(color.scheme) <- as.character(sort(unique(na.omit(phenotype.values))));
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
            continuous.phenotypes.df <- continuous.phenotypes.df[ , order(pgs.data$percentile)];

            # save NA coordinates for labeling
            cont.pheno.na.coords <- get.na.coordinates.for.heatmap(continuous.phenotypes.df);

            # replace NA values with a color
            continuous.phenotypes.df[is.na(continuous.phenotypes.df)] <- FILL.COLOR;


            continuous.phenotype.heatmap <- BoutrosLab.plotting.general::create.heatmap(
                x = continuous.phenotypes.df,
                input.colours = TRUE,
                clustering.method = 'none',
                same.as.matrix = TRUE,
                print.colour.key = FALSE,
                cell.text = 'NA',
                col.pos = cont.pheno.na.coords$col,
                row.pos = cont.pheno.na.coords$row,
                text.col = 'white',
                text.cex = 0.5,
                yaxis.lab = NULL,
                ylab.cex = 1
                # main.cex = titles.cex,
                # ylab.cex = titles.cex,
                # xaxis.cex = xaxis.cex,
                # yaxis.cex = yaxis.cex
                );

            # legend
            continuous.covariates.legend <- lapply(
                X = 1:length(continuous.color.schemes),
                FUN = function(x) {
                    list(
                        title = names(continuous.color.schemes)[x],
                        colours = continuous.color.schemes[[x]],
                        labels = round(c(min(continuous.phenotype.data[ , x], na.rm = TRUE), max(continuous.phenotype.data[ , x], na.rm = TRUE)), 1),
                        continuous = TRUE
                        );
                    }
                );
            names(continuous.covariates.legend) <- rep('legend', length(continuous.covariates.legend));

            }
    }
    ## End Phenotype Covariate Heatmap Assembly ##

    # organize filename if plot writing requested
    if (!is.null(output.dir)) {

        if (is.null(filename.prefix)) {
            filename.prefix <- 'ApplyPolygenicScore-Plot';
            }
        # construct multipanel plot
        filename.for.rank.multiplot <- generate.filename(
            project.stem = filename.prefix,
            file.core = 'pgs-rank-plot',
            extension = file.extension
            );

        output.path <- file.path(output.dir, filename.for.rank.multiplot);
        } else {
            output.path <- NULL;
        }

    

    # assemble plot legend
    cov.legends <- c(
        percentile.covariates.legend,
        binary.covariates.legend,
        categorical.covariates.legend,
        continuous.covariates.legend
        );

    cov.legend.grob <- BoutrosLab.plotting.general::legend.grob(
        cov.legends
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
        filename = output.path,
        main = '',
        main.cex = 0,
        layout.height = length(plot.list),
        layout.width = 1,
        plot.objects.heights = plot.heights,
        y.spacing = -2,
        ylab.axis.padding = -6,
        legend = list(right = list(fun = cov.legend.grob)),
        width = width,
        height = height,
        right.padding = border.padding,
        left.padding = border.padding,
        top.padding = border.padding,
        bottom.padding = border.padding
        );

    return(multipanel.plot);
    }
