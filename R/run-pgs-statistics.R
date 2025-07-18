#' @title get.pgs.percentiles
#' @description Calculate percentiles and report decile and quartile ranks for a vector of polygenic scores
#' @param pgs numeric vector of polygenic scores
#' @param n.percentiles integer number of percentiles to calculate (optional)
#' @return data frame with columns for percentile, decile, quartile, and optional n.percentiles
#' @examples
#' x <- rnorm(100);
#' get.pgs.percentiles(x, n.percentiles = 20);
#' @export
get.pgs.percentiles <- function(pgs, n.percentiles = NULL) {

    # check that pgs is numeric
    if (!is.numeric(pgs)) {
        stop('pgs must be a numeric vector');
        }

    # check that n.percentiles is a mathematical integer
    if (!is.null(n.percentiles) && (n.percentiles %% 1 != 0)) {
        stop('n.percentiles must be an integer');
        }

    # save the original order of the PGS
    pgs.percentile.data <- data.frame(pgs = pgs, order = 1:length(pgs));

    # calculate percentiles
    pgs.percentile.data$percentile <- sapply(
        X = pgs,
        FUN = function(x) {
            sum(x >= pgs, na.rm = TRUE) / length(pgs);
            }
        );

    # sort pgs percentiles in descending order to facilitate decile and quartile calculations
    pgs.percentile.data <- pgs.percentile.data[order(pgs.percentile.data$percentile, decreasing = TRUE),];

    # compute which decile each PGS belongs to
    pgs.percentile.data$decile <- cut(pgs.percentile.data$percentile, breaks = seq(0, 1, 0.1), labels = FALSE);

    # compute which quartile each pgs belongs to
    pgs.percentile.data$quartile <- cut(pgs.percentile.data$percentile, breaks = seq(0, 1, 0.25), labels = FALSE);

    if (!is.null(n.percentiles)) {
        # calculate user-specified percentiles
        pgs.percentile.data$percentile.X <- cut(pgs.percentile.data$percentile, breaks = seq(0, 1, 1 / n.percentiles), labels = FALSE);
        # replace column name with the number of percentiles
        colnames(pgs.percentile.data) <- gsub('percentile.X', paste0('percentile.', n.percentiles), colnames(pgs.percentile.data));
        }

    # go back to the original order
    pgs.percentile.data <- pgs.percentile.data[order(pgs.percentile.data$order),];

    # remove extra columns
    pgs.percentile.data$order <- NULL;
    pgs.percentile.data$pgs <- NULL;

    return(pgs.percentile.data);
    }

# utility function for identifying data as continuous or binary for analysis and plotting purposes
classify.variable.type <- function(data, continuous.threshold = 4) {
    # identify continuous and binary variables
    continuous.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            ('numeric' == class(x) | 'integer' == class(x)) & continuous.threshold < length(unique(na.omit(x)));
            }
        );
    binary.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            2 == length(unique(na.omit(x)));
            }
        );

    other.vars.index <- !continuous.vars.index & !binary.vars.index;

    return(list(
        continuous = continuous.vars.index,
        binary = binary.vars.index,
        other = other.vars.index
        ));
}

# utility function that runs linear and logistic regression on a polygenic score and a set of phenotypes
#' @title Run linear and logistic regression on a polygenic score and a set of phenotypes
#' @description Phenotype data variables are automatically classified as continuous or binary and a simple linear regression or logistic regression, respectively, is run between the polygenic score and each phenotype.
#' Categorical phenotypes with more than two category are ignored.
#' If a binary variable is not formatted as a factor, it is converted to a factor using \code{as.factor()} defaults. For logistic regression, the first level is classified as "failure" and the second "success" by \code{glm()} defaults.
#' @param pgs numeric vector of polygenic scores
#' @param phenotype.data data.frame of phenotypes
#' @return data frame with columns for phenotype, model, beta, se, p.value, r.squared, and AUC
#' @examples
#' set.seed(200);
#' pgs <- rnorm(200, 0, 1);
#' phenotype.data <- data.frame(
#'     continuous.pheno = rnorm(200, 1, 1),
#'     binary.pheno = sample(c(0, 1), 200, replace = TRUE)
#'     );
#' run.pgs.regression(pgs, phenotype.data);
#' @export
run.pgs.regression <- function(pgs, phenotype.data) {

    # initialize conditional outputs
    linear.model.aggregated <- NULL;
    logistic.model.aggregated <- NULL;

    # identify continuous and binary phenotypes
    variable.index.by.type <- classify.variable.type(phenotype.data);

    # run linear regression on continuous phenotypes
    if (any(variable.index.by.type$continuous)) {
        continuous.data <- subset(phenotype.data, select = variable.index.by.type$continuous);
        linear.model <- lapply(
            X = continuous.data,
            FUN = function(x) {
                return(summary(lm(x ~ pgs, data = phenotype.data)));
                }
            );

        # aggregate results in a data frame
        linear.model.aggregated <- lapply(
            X = linear.model,
            FUN = function(x) {
                coeff.index <- if (nrow(x$coefficients) == 1) NA else 'pgs';
                data.frame(
                    beta = x$coefficients[coeff.index, 'Estimate'],
                    se = x$coefficients[coeff.index, 'Std. Error'],
                    p.value = x$coefficients[coeff.index, 'Pr(>|t|)'],
                    r.squared = x$r.squared,
                    AUC = NA
                    );
                }
            );
        linear.model.aggregated <- do.call(rbind, linear.model.aggregated);
        linear.model.aggregated <- data.frame(
            phenotype = names(linear.model),
            model = 'linear.regression',
            linear.model.aggregated
            );

    }


    # run logistic regression on binary phenotypes
    if (any(variable.index.by.type$binary)) {

        binary.data <- subset(phenotype.data, select = variable.index.by.type$binary);

        # ensure binary data is formatted as factors
        binary.data <- lapply(
            X = binary.data,
            FUN = function(x) {
                if (!is.factor(x)) {
                    return(as.factor(x));
                    } else {
                    return(x);
                    }
                }
            );

        logistic.model <- lapply(
            X = binary.data,
            FUN = function(x) {
                return(glm(x ~ pgs, data = binary.data, family = binomial));
                }
            );

        logistic.model.summary <- lapply(
            X = logistic.model,
            FUN = summary
            );

        logistic.auc <- lapply(
            X = logistic.model,
            FUN = function(x) {
                predictions <- predict(x, type = 'response'); # get predicted probabilities on the training set
                # auc is a bit wordy by default
                auc <- suppressMessages(pROC::auc(x$y, predictions)); # compute area under the curve on training set
                }
            );


        logistic.model.aggregated <- lapply(
            X = logistic.model.summary,
            FUN = function(x) {
                coeff.index <- if (nrow(x$coefficients) == 1) NA else 'pgs';
                data.frame(
                    beta = x$coefficients[coeff.index, 'Estimate'],
                    se = x$coefficients[coeff.index, 'Std. Error'],
                    p.value = x$coefficients[coeff.index, 'Pr(>|z|)'],
                    r.squared = NA
                    );
                }
            );
        logistic.model.aggregated <- do.call(rbind, logistic.model.aggregated);
        logistic.model.aggregated$AUC <- unlist(logistic.auc);
        logistic.model.aggregated <- data.frame(
            phenotype = names(logistic.model),
            model = 'logistic.regression',
            logistic.model.aggregated
            );
    }

    all.model.results <- rbind(linear.model.aggregated, logistic.model.aggregated);

    return(all.model.results);
    }

analyze.pgs.binary.predictiveness <- function(
    data,
    pgs.columns,
    phenotype.columns,
    covariate.columns = NULL,
    phenotype.type = c('binary', 'continuous'),
    cutoff.threshold = NULL,
    output.dir = NULL,
    filename.prefix = NULL,
    file.extension = 'png',
    width = 8,
    height = 8,
    xaxes.cex = 1.5,
    yaxes.cex = 1.5,
    titles.cex = 1.5
    ) {
    ## Validate Inputs ##

    if (!is.data.frame(data)) {
        stop('`data` must be a data frame.');
        }

    if (!all(pgs.columns %in% names(data))) {
        stop('Not all specified `pgs.columns` found in the data frame.');
        }
    if (!all(sapply(data[, pgs.columns, drop = FALSE], is.numeric))) {
        stop('All `pgs.columns` must be numeric.');
        }

    if (!all(phenotype.columns %in% names(data))) {
        stop(paste0('Not all specified `phenotype.columns` found in the data frame.'));
        }

    if (!is.null(covariate.columns) && !all(covariate.columns %in% names(data))) {
        stop('Not all specified `covariate.columns` found in the data frame.');
        }

    if (!is.character(phenotype.type)) {
        stop('`phenotype.type` must be either "binary" or "continuous"')
        }

    if (length(phenotype.type) > 1) {
        stop('`phenotype.type` must be either "binary" or "continuous"')
        }

     if (!phenotype.type %in% c('binary', 'continuous')) {
        stop('`phenotype.type` must be either "binary" or "continuous"')
        }

    # Validate phenotype types consistency
    for (pheno.col in phenotype.columns) {
        if (phenotype.type == 'binary') {
            if (!is.factor(data[[pheno.col]]) && !all(unique(na.omit(data[[pheno.col]])) %in% c(0, 1))) {
                stop(paste0('Phenotype column \'', pheno.col, '\' is specified as binary but is not a factor or 0/1 numeric. Convert to factor.'));
                }
            if (is.factor(data[[pheno.col]]) && nlevels(data[[pheno.col]]) != 2) {
                stop(paste0('Binary phenotype column \'', pheno.col, '\' must have exactly two levels.'));
                }
        } else { # phenotype.type == 'continuous'
            if (!is.numeric(data[[pheno.col]])) {
                stop(paste0('Phenotype column \'', pheno.col, '\' is specified as continuous but is not numeric.'));
                }
            if (is.null(cutoff.threshold)) {
                stop('`cutoff.threshold` must be provided for `continuous` phenotype.type.');
                }
            if (is.list(cutoff.threshold) && !pheno.col %in% names(cutoff.threshold)) {
                stop(paste0('`cutoff.threshold` must contain an entry for \'', pheno.col, '\' when it is a named list.'));
                }
            if (!is.list(cutoff.threshold) && !is.numeric(cutoff.threshold)) {
                stop('`cutoff.threshold` must be numeric or a named list for `continuous` phenotype.type.');
                }
            }
        }

    # Prepare results data frame
    results.df <- data.frame(
        phenotype = character(),
        PGS = character(),
        AUC = numeric(),
        OR = numeric(),
        OR.Lower.CI = numeric(),
        OR.Upper.CI = numeric(),
        p.value = numeric(),
        stringsAsFactors = FALSE
        );

    # Prepare plot variable
    roc.plot <- NULL;

    # Store ROC data for plotting
    roc.data.for.plotting <- list();

    # Outer loop, iterate through each phenotype
    for (current.phenotype in phenotype.columns) {
        # message(paste0('Processing phenotype: ', current.phenotype));

        temp.data <- data; # Work on a fresh copy for each phenotype

        pheno.var <- temp.data[[current.phenotype]];

        if (phenotype.type == 'continuous') {
            # parse cutoff threshold for this phenotype
            current.cutoff <- if (is.list(cutoff.threshold)) cutoff.threshold[[current.phenotype]] else cutoff.threshold;
            # check that valid threshold was provided
            if (is.null(current.cutoff) || !is.numeric(current.cutoff)) {
                stop(paste0('Cutoff threshold for phenotype \'', current.phenotype, '\' is missing or not numeric.'));
                }

            # convert continuous pheno to binary; below threshold = 0, above threshold = 1
            pheno.var <- as.integer(pheno.var >= current.cutoff);
            temp.data[[current.phenotype]] <- factor(pheno.var, levels = c(0, 1));

        } else { # phenotype.type == 'binary'
            if (!is.factor(pheno.var)) {
                warning(paste0('Phenotype column \'', current.phenotype, '\' is not a factor. Attempting to convert to factor.'));
                pheno.var <- as.factor(pheno.var);
                temp.data[[current.phenotype]] <- pheno.var;
                }
            if (nlevels(pheno.var) > 2) {
                stop(paste0('Phenotype column \'', current.phenotype, '\' has more than two levels. Please reconfigure.'))
                }
            }

        # Inner loop: iterate through each PGS for the current phenotype
        for (pgs in pgs.columns) {
            # Define columns relevant for the current model
            model.cols <- c(current.phenotype, pgs);
            if (!is.null(covariate.columns)) {
                model.cols <- c(model.cols, covariate.columns);
                }

            # Filter data for complete cases within the current model
            clean.data <- temp.data[complete.cases(temp.data[, model.cols, drop = FALSE]), ];

            # If no complete cases, skip to next PGS
            if (nrow(clean.data) == 0) {
                warning(paste0('No complete cases for PGS \'', pgs, '\' and phenotype \'', current.phenotype, '\' after NA removal. Skipping.'));
                results.df <- rbind(results.df, data.frame(
                    phenotype = current.phenotype, PGS = pgs, AUC = NA, OR = NA,
                    OR.Lower.CI = NA, OR.Upper.CI = NA, p.value = NA
                    ));
                next; # Skip to next PGS
                }

             # Construct formula for logistic regression
            formula.str <- paste(current.phenotype, '~', pgs);
            if (!is.null(covariate.columns)) {
                formula.str <- paste(formula.str, '+', paste(covariate.columns, collapse = ' + '));
                }
            model.formula <- as.formula(formula.str);

            # Run logistic regression, gracefully move on if it doesn't work
            model <- tryCatch({
                glm(model.formula, data = clean.data, family = binomial(link = 'logit'));
                }, error = function(e) {
                warning(paste0('Logistic regression failed for PGS \'', pgs, '\' and phenotype \'', current.phenotype, '\': ', e$message,
                               '. Skipping AUC and OR calculation for this PGS/phenotype pair.'));
                return(NULL);
                });

            # If model succeds, compute AUC
            if (!is.null(model)) {
                # Calculate AUC
                predicted.probs <- predict(model, type = 'response');

                roc.obj <- tryCatch({
                    roc(response = clean.data[[current.phenotype]],
                        predictor = predicted.probs,
                        levels = levels(clean.data[[current.phenotype]]), # Ensure correct order of levels
                        direction = '<', # If higher predicted.prob means higher risk, use '<'
                        ci = TRUE, # Calculate confidence interval for AUC
                        quiet = TRUE
                        );
                }, error = function(e) {
                    warning(paste0('ROC curve calculation failed for PGS \'', pgs, '\' and phenotype \'', current.phenotype, '\': ', e$message,
                                   '. Skipping ROC and AUC for this PGS/phenotype pair.'));
                    return(NULL);
                    });

            # If ROC succeeds, extract ROC coordinates
            if (!is.null(roc.obj)) {
                    current.auc <- auc(roc.obj);

                    # Extract FPR and TPR for plotting
                    roc.coords <- coords(roc.obj, 'all', ret = c('fpr', 'tpr'), drop = FALSE);
                    roc.coords$PGS <- pgs;
                    roc.coords$phenotype <- current.phenotype;
                    roc.data.for.plotting[[paste0(pgs, '_', current.phenotype)]] <- roc.coords;

                    # Get Odds Ratios and CIs
                    coefs <- coef(model);
                    vcov.matrix <- vcov(model);
                    std.errors <- sqrt(diag(vcov.matrix));
                    z.values <- coefs / std.errors;
                    p.values <- 2 * pnorm(-abs(z.values));

                    # Check if the PGS term exists in the model coefficients
                    if (pgs %in% names(coefs)) {

                        if (is.na(coefs[pgs]) || is.infinite(coefs[pgs])) {
                            warning(paste0('PGS term \'', pgs, '\' in model for phenotype \'', current.phenotype, '\' is NA or Inf (e.g., due to no variance or perfect separation). Skipping AUC and OR calculation for this PGS/phenotype pair.'));
                            results.df <- rbind(results.df, data.frame(
                                phenotype = current.phenotype, PGS = pgs, AUC = NA, OR = NA,
                                OR.Lower.CI = NA, OR.Upper.CI = NA, p.value = NA
                                ));
                            next; # Skip to next PGS
                            }

                        pgs.index <- which(names(coefs) == pgs);
                        current.beta <- coefs[pgs.index];
                        current.se <- std.errors[pgs.index];
                        current.p.value <- p.values[pgs.index];

                        # 95% Confidence Interval for beta (log odds)
                        # For a 95% CI, Z-score is approx 1.96
                        lower.beta.ci <- current.beta - 1.96 * current.se;
                        upper.beta.ci <- current.beta + 1.96 * current.se;

                        # Exponentiate to get Odds Ratio and its CIs
                        current.or <- exp(current.beta);
                        current.or.lower.ci <- exp(lower.beta.ci);
                        current.or.upper.ci <- exp(upper.beta.ci);

                        results.df <- rbind(results.df, data.frame(
                            phenotype = current.phenotype,
                            PGS = pgs,
                            AUC = current.auc,
                            OR = current.or,
                            OR.Lower.CI = current.or.lower.ci,
                            OR.Upper.CI = current.or.upper.ci,
                            p.value = current.p.value
                            ));
                    } else { # if no pgs term found in coefficients
                        warning(paste0('PGS term \'', pgs, '\' not found in model summary for OR calculation for phenotype \'', current.phenotype, '\'. This can happen if the PGS is collinear with other predictors or has no variance.'));
                        results.df <- rbind(results.df, data.frame(
                                phenotype = current.phenotype,
                                PGS = pgs,
                                AUC = current.auc,
                                OR = NA,
                                OR.Lower.CI = NA,
                                OR.Upper.CI = NA,
                                p.value = NA
                            ));
                        } 
                    } else { # if ROC doesn't compute
                    results.df <- rbind(results.df, data.frame(
                        phenotype = current.phenotype,
                        PGS = pgs,
                        AUC = NA,
                        OR = NA,
                        OR.Lower.CI = NA,
                        OR.Upper.CI = NA,
                        p.value = NA
                        ));
                    }
            } else { # If model itself failed
               results.df <- rbind(results.df, data.frame(
                    phenotype = current.phenotype,
                    PGS = pgs,
                    AUC = NA,
                    OR = NA,
                    OR.Lower.CI = NA,
                    OR.Upper.CI = NA,
                    p.value = NA
                    ));
                }
            } # End PGS loop
        } # End phenotype loop

    ## PLOTTING ##
    roc.plot <- NULL;

    if (length(roc.data.for.plotting) > 0) {
        # Combine all ROC data into a single data frame for plotting
        plot.data <- do.call(rbind, roc.data.for.plotting);

        plot.objects.list <- list();
        panel.counter <- 1;

        # Pre-calculate colors for all unique PGSs (to maintain consistent coloring across plots)
        max.colors = 12;
        if (length(pgs.columns) < max.colors) {
            pgs.line.colours <- BoutrosLab.plotting.general::default.colours(length(pgs.columns));
        } else {
            if (requireNamespace('scales', quietly = TRUE)) {
                pgs.line.colours <- scales::hue_pal()(length(pgs.columns));
            } else {
                stop('Too many PGS columns selected, must be 12 or fewer');
                }
        }
        names(pgs.line.colours) <- pgs.columns;

        # Outer loop for plotting: one plot per phenotype
        for (pheno.name in phenotype.columns) {
            # Filter data for the current phenotype, including all PGSs
            current.pheno.plot.data <- plot.data[plot.data$phenotype == pheno.name, ];

            # Ensure 'PGS' is a factor in this subset for grouping
            current.pheno.plot.data$PGS <- factor(current.pheno.plot.data$PGS, levels = pgs.columns);

            # If there's no data for this phenotype (e.g., all NAs filtered out earlier)
            if (nrow(current.pheno.plot.data) == 0) {
                warning(paste0('No plot data for phenotype \'', pheno.name, '\'. Creating empty plot.'));
                plot.objects.list[[panel.counter]] <- BoutrosLab.plotting.general::create.scatterplot(
                    formula = 0 ~ 0, # Dummy formula
                    data = data.frame(fpr = 0, tpr = 0),
                    xlab.label = 'False Positive Rate (1 - Specificity)',
                    ylab.label = 'True Positive Rate (Sensitivity)',
                    main = paste0('ROC: ', pheno.name, ' (No Data)'),
                    main.cex = titles.cex,
                    xat = seq(0, 1, 0.2), yat = seq(0, 1, 0.2),
                    xlimits = c(0, 1), ylimits = c(0, 1),
                    xaxis.cex = xaxes.cex, yaxis.cex = yaxes.cex,
                    xlab.cex = titles.cex, ylab.cex = titles.cex
                    );
                } else { # If there is data for this phenotype

                # Construct legend items for this plot
                legend.items <- list();
                for (pgs.name.in.plot in levels(current.pheno.plot.data$PGS)) {
                    # Get AUC for the specific PGS and phenotype
                    auc.val.row <- results.df[results.df$phenotype == pheno.name & results.df$PGS == pgs.name.in.plot, 'AUC'];
                    auc.val <- if (length(auc.val.row) > 0) round(auc.val.row[1], 3) else NA;
                    legend.label <- paste0(pgs.name.in.plot, ' (AUC = ', ifelse(is.na(auc.val), 'N/A', auc.val), ')');

                    legend.items[[length(legend.items) + 1]] <- list(
                        color = pgs.line.colours[pgs.name.in.plot],
                        label = legend.label
                        );
                    } # End inner legend loop

                # Extract all colors and labels into single vectors as required by 'key'
                all.legend.colors <- sapply(legend.items, `[[`, 'color');
                all.legend.labels <- sapply(legend.items, `[[`, 'label');

                pgs.legend <- list(
                        inside = list(
                            fun = lattice::draw.key,
                            args = list(
                                key = list(
                                    lines = list(
                                        col = all.legend.colors,
                                        lwd = 2, # Assuming constant lwd for all lines
                                        size = 1.5
                                        ),
                                    text = list(
                                        lab = all.legend.labels, # All labels in one vector
                                        cex = 0.8
                                        ),
                                    columns = 1
                                    )
                                ),
                            corner = c(0.05, 0.95), # Top-left corner for legend
                            x = 0.05, y = 0.95 # Specific coordinates within panel for legend
                            )
                    );

                current.plot.object <- BoutrosLab.plotting.general::create.scatterplot(
                    formula = tpr ~ fpr,
                    data = current.pheno.plot.data,
                    groups = current.pheno.plot.data$PGS, # Group by PGS
                    type = 'l', # Line plot for ROC curve
                    col = pgs.line.colours[levels(current.pheno.plot.data$PGS)], # Colors mapped to groups
                    legend = pgs.legend,
                    lwd = 2,
                    xlab.label = 'False Positive Rate (1 - Specificity)',
                    ylab.label = 'True Positive Rate (Sensitivity)',
                    main = paste0('ROC Curve: ', pheno.name), # Title for individual panel
                    main.cex = titles.cex,
                    xat = seq(0, 1, 0.2),
                    yat = seq(0, 1, 0.2),
                    xlimits = c(0, 1),
                    ylimits = c(0, 1),
                    xaxis.cex = xaxes.cex,
                    yaxis.cex = yaxes.cex,
                    xlab.cex = titles.cex,
                    ylab.cex = titles.cex,
                    add.xyline = TRUE,
                    xyline.col = 'grey50',
                    xyline.lwd = 1,
                    xyline.lty = 2
                    );

                plot.objects.list[[panel.counter]] <- current.plot.object;

                } # End condition there is data for this phenotype
            panel.counter <- panel.counter + 1;

            } # End outer plotting loop

        # organize filename if plot writing requested
        if (!is.null(output.dir)) {

            if (is.null(filename.prefix)) {
                filename.prefix <- 'ApplyPolygenicScore-Plot';
                }
            # construct multipanel plot
            filename.for.roc.multiplot <- generate.filename(
                project.stem = filename.prefix,
                file.core = 'pgs-roc-curves',
                extension = file.extension
                );

            output.path <- file.path(output.dir, filename.for.roc.multiplot);
            } else {
                output.path <- NULL;
            }

        # Create a single multipanel plot for all ROC curves
        roc.multipanel <- BoutrosLab.plotting.general::create.multipanelplot(
            plot.objects = plot.objects.list,
            filename = output.path,
            layout.height = length(phenotype.columns), # One row per phenotype
            layout.width = 1, # Two columns of plots
            main.cex = 0,
            width = width,
            height = height,
            x.spacing = 0.5,
            y.spacing = 0.5
            );

        roc.plot <- roc.multipanel;
    } else {# If no plotting data exists
        warning('No ROC curves could be generated for plotting')
        } 

    return(list(
        results.df = results.df,
        roc.plot = roc.plot
        ));

    } # End function