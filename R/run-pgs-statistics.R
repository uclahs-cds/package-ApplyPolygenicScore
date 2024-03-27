#' @title get.pgs.percentiles
#' @description Calculate percentiles and report decile and quartile ranks for a vector of polygenic scores
#' @param pgs numeric vector of polygenic scores
#' @param n.percentiles integer number of percentiles to calculate (optional)
#' @return data frame with columns for percentile, decile, quartile, and optional n.percentiles
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
classify.variable.type <- function(data) {
    # identify continuous and binary variables
    continuous.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            'numeric' == class(x) & 2 < length(unique(x));
            }
        );
    binary.vars.index <- sapply(
        X = data,
        FUN = function(x) {
            2 == length(unique(x));
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
run.pgs.regression <- function(pgs, phenotype.data) {

    # identify continuous and binary phenotypes
    variable.index.by.type <- classify.variable.type(phenotype.data);

    # run linear regression on continuous phenotypes
    continuous.data <- subset(phenotype.data, select = variable.index.by.type$continuous);
    linear.model <- lapply(
        X = continuous.data,
        FUN = function(x) {
            summary(lm(x ~ pgs, data = phenotype.data));
            }
        );

    # run logistic regression on binary phenotypes
    binary.data <- subset(phenotype.data, select = variable.index.by.type$binary);
    binary.data <- lapply(
        X = binary.data,
        FUN = function(x) {
            factor(x, levels = c(0, 1));
            }
        );

    logistic.model <- lapply(
        X = binary.data,
        FUN = function(x) {
            summary(glm(x ~ pgs, data = binary.data, family = binomial));
            }
        );

    # aggregate results in a data frame
    linear.model.aggregated <- lapply(
        X = linear.model,
        FUN = function(x) {
            data.frame(
                beta = x$coefficients['pgs', 'Estimate'],
                se = x$coefficients['pgs', 'Std. Error'],
                p.value = x$coefficients['pgs', 'Pr(>|t|)'],
                r.squared = x$r.squared
                );
            }
        );
    linear.model.aggregated <- do.call(rbind, linear.model.aggregated);
    linear.model.aggregated <- data.frame(
        phenotype = names(linear.model),
        model = 'linear.regression',
        linear.model.aggregated
        );

    logistic.model.aggregated <- lapply(
        X = logistic.model,
        FUN = function(x) {
            data.frame(
                beta = x$coefficients['pgs', 'Estimate'],
                se = x$coefficients['pgs', 'Std. Error'],
                p.value = x$coefficients['pgs', 'Pr(>|z|)'],
                r.squared = NA
                );
            }
        );
    logistic.model.aggregated <- do.call(rbind, logistic.model.aggregated);
    logistic.model.aggregated <- data.frame(
        phenotype = names(logistic.model),
        model = 'logistic.regression',
        logistic.model.aggregated
        );

    all.model.results <- rbind(linear.model.aggregated, logistic.model.aggregated);

    return(all.model.results);
    }
