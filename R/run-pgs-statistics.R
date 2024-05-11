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

    sorted.pgs <- sort(pgs, decreasing = TRUE)
    pgs.percentiles <- sapply(
        X = sorted.pgs,
        FUN = function(x) {
            sum(pgs >= x, na.rm = TRUE) / length(pgs);
            }
        );

    # convert to data frame
    pgs.percentile.data <- data.frame(
        percentile = pgs.percentiles
        );

    # compute which decile each PGS belongs to
    pgs.percentile.data$decile <- cut(pgs.percentiles, breaks = seq(0, 1, 0.1), labels = FALSE);

    # compute which quartile each pgs belongs to
    pgs.percentile.data$quartile <- cut(pgs.percentiles, breaks = seq(0, 1, 0.25), labels = FALSE);

    if (!is.null(n.percentiles)) {
        # calculate user-specified percentiles
        pgs.percentile.data$percentile.X <- cut(pgs.percentiles, breaks = seq(0, 1, 1 / n.percentiles), labels = FALSE);
        # replace column name with the number of percentiles
        colnames(pgs.percentile.data) <- gsub('percentile.X', paste0('percentile.', n.percentiles), colnames(pgs.percentile.data));
        }

    # go back to the original order
    pgs.percentile.data <- pgs.percentile.data[order(pgs),];
    return(pgs.percentile.data);
    }

# utility function for identifying data as continuous or binary for analysis and plotting purposes
classify.variable.type <- function(data, continuous.threshold = 4) {
    # identify continuous and binary variables
    continuous.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            'numeric' == class(x) & continuous.threshold < length(unique(na.omit(x)));
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
#' If a binary variable is not formatted as a factor, it will be converted to a factor using as.factor() defaults. For logistic regression, the first level is classified as "failure" and the second "success" by glm defaults.
#' @param pgs numeric vector of polygenic scores
#' @param phenotype.data data.frame of phenotypes
#' @return data frame with columns for phenotype, model, beta, se, p.value, r.squared, and AUC
#' @examples
#' set.seed(200);
#' pgs <- rnorm(200, 0, 1);
#' phenotype.data <- data.frame(
#'  continuous.pheno = rnorm(200, 1, 1),
#'  binary.pheno = sample(c(0, 1), 200, replace = TRUE)
#'  );
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
