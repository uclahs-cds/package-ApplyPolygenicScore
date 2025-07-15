# Handling CRAN warnings for data.table syntax:
if (getRversion() >= '2.15.1') utils::globalVariables(c(
    'model', 'phenotype', 'se', 'p.value',
    'r.squared', 'AUC'
    ));

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
    data <- data.table::as.data.table(data);
    # identify continuous and binary variables
    continuous.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            ('numeric' == class(x) | 'integer' == class(x)) & (continuous.threshold < data.table::uniqueN(x, na.rm = TRUE));
            }
        );
    binary.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            (2 == data.table::uniqueN(x, na.rm = TRUE));
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
    # Ensure phenotype.data is a data.table for efficient operations
    data.table::setDT(phenotype.data);

    # initialize conditional outputs
    linear.model.aggregated <- NULL;
    logistic.model.aggregated <- NULL;

    # identify continuous and binary phenotypes
    variable.index.by.type <- classify.variable.type(phenotype.data);

    # run linear regression on continuous phenotypes
    if (any(variable.index.by.type$continuous)) {
        continuous.data <- phenotype.data[, .SD, .SDcols = variable.index.by.type$continuous];
        linear.model <- lapply(
            X = colnames(continuous.data), # Iterate over column names for robustness
            FUN = function(pheno.name) {
                # Create a temporary data.table for lm to find the variables 'y' and 'pgs.val'
                temp.data <- data.table::data.table(
                    y = continuous.data[[pheno.name]], # Extract column content by name
                    pgs.val = pgs                      # The global pgs vector
                    );
                return(summary(lm(y ~ pgs.val, data = temp.data)));
                }
            );
        names(linear.model) <- colnames(continuous.data); # Reassign names to the list elements


        # aggregate results in a data frame
        linear.model.aggregated <- lapply(
            X = linear.model,
            FUN = function(x) {
                # Ensure 'pgs.val' is in the coefficients table; if only intercept, handle NA
                coeff.index <- if ("pgs.val" %in% rownames(x$coefficients)) "pgs.val" else NA;
                data.table::data.table( # Create data.table directly
                    beta = x$coefficients[coeff.index, 'Estimate'],
                    se = x$coefficients[coeff.index, 'Std. Error'],
                    p.value = x$coefficients[coeff.index, 'Pr(>|t|)'],
                    r.squared = x$r.squared,
                    AUC = NA_real_ # Use NA_real_ for numeric NA
                    );
                }
            );
        # Combine list of data.tables into one using rbindlist, adding phenotype column
        linear.model.aggregated <- data.table::rbindlist(linear.model.aggregated, idcol = "phenotype");
        linear.model.aggregated[, model := 'linear.regression']; # Add model column by reference
        
        # Ensure column order
        linear.model.aggregated <- linear.model.aggregated[, .(phenotype, model, beta, se, p.value, r.squared, AUC)];

    }


    # Run logistic regression on binary phenotypes
    if (any(variable.index.by.type$binary)) {

        # Select binary phenotype columns as a data.table using .SDcols
        binary.data <- phenotype.data[, .SD, .SDcols = variable.index.by.type$binary];

        # Ensure binary data is formatted as factors using data.table::set for in-place modification
        for (pheno.name in colnames(binary.data)) {
            if (!is.factor(binary.data[[pheno.name]])) {
                data.table::set(binary.data, j = pheno.name, value = as.factor(binary.data[[pheno.name]]));
                }
            }

        # Run glm for each binary phenotype
        logistic.model <- lapply(
            X = colnames(binary.data), # Iterate over column names
            FUN = function(pheno.name) {
                # Create a temporary data.table for glm to find the variables 'y' and 'pgs.val'
                temp.data <- data.table::data.table(
                    y = binary.data[[pheno.name]], # Extract column content by name
                    pgs.val = pgs                      # The global pgs vector
                    );
                return(glm(y ~ pgs.val, data = temp.data, family = binomial));
                }
            );
        names(logistic.model) <- colnames(binary.data); # Reassign names

        logistic.model.summary <- lapply(
            X = logistic.model,
            FUN = summary
            );

        # Calculate AUC for logistic models
        logistic.auc <- lapply(
            X = logistic.model,
            FUN = function(x) {
                predictions <- predict(x, type = 'response'); # Get predicted probabilities on the training set
                # Compute area under the curve on training set. Ensure pROC package is available.
                return(suppressMessages(pROC::auc(x$y, predictions)));
                }
            );


        # Aggregate logistic model results into a data.table
        logistic.model.aggregated <- lapply(
            X = logistic.model.summary,
            FUN = function(x) {
                # Ensure 'pgs.val' is in the coefficients table; if only intercept, handle NA
                coeff.index <- if ("pgs.val" %in% rownames(x$coefficients)) "pgs.val" else NA;
                data.table::data.table( # Create data.table directly
                    beta = x$coefficients[coeff.index, 'Estimate'],
                    se = x$coefficients[coeff.index, 'Std. Error'],
                    p.value = x$coefficients[coeff.index, 'Pr(>|z|)'],
                    r.squared = NA_real_ # Use NA_real_ for numeric NA
                    );
                }
            );
        # Combine list of data.tables into one using rbindlist, adding phenotype column
        logistic.model.aggregated <- data.table::rbindlist(logistic.model.aggregated, idcol = "phenotype");
        logistic.model.aggregated[, model := 'logistic.regression']; # Add model column by reference

        # Add AUC column by reference
        logistic.model.aggregated[, AUC := unlist(logistic.auc)];
        
        # Ensure column order
        logistic.model.aggregated <- logistic.model.aggregated[, .(phenotype, model, beta, se, p.value, r.squared, AUC)];
    }


        # Combine all model results (linear and logistic)
        # Using fill = TRUE to handle cases where columns might be missing (e.g., AUC for linear, r.squared for logistic)
        all.model.results <- data.table::rbindlist(
            list(linear.model.aggregated, logistic.model.aggregated),
            fill = TRUE
        );

    return(as.data.frame(all.model.results));
    }
