# function that returns the percentile of each PGS within the PGS distribution
# also returns Decile (10% bins) and Quantile (25% bins) of the PGS distribution
get.pgs.percentiles <- function(pgs, n.percentiles = NULL) {
    sorted.pgs <- sort(pgs, decreasing = TRUE)
    pgs.percentiles <- sapply(
        X = sorted.pgs,
        FUN = function(x) {
            sum(pgs >= x) / length(pgs)
            }
        )

    # convert to data frame
    pgs.percentile.data <- data.frame(
        percentile = pgs.percentiles
        );

    # compute which decile each PGS belongs to
    pgs.percentile.data$decile <- cut(pgs.percentiles, breaks = seq(0, 1, 0.1), labels = FALSE);

    # compute which quartile each pgs belongs to
    pgs.percentile.data$quartile <- cut(pgs.percentiles, breaks = seq(0, 1, 0.25), labels = FALSE);

    if (!is.null(n.percentiles)) {
        pgs.percentile.data$percentile.X <- cut(pgs.percentiles, breaks = seq(0, 1, 1/n.percentiles), labels = FALSE);
        # replace column name with the number of percentiles
        colnames(pgs.percentile.data) <- gsub('percentile.X', paste0('percentile.', n.percentiles), colnames(pgs.percentile.data));
        }

    # go back to the original order
    pgs.percentile.data <- pgs.percentile.data[order(pgs),];
    return(pgs.percentile.data);
    }