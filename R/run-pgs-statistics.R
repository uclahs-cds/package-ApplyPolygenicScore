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
            sum(pgs >= x) / length(pgs);
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
