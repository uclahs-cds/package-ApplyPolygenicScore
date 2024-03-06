# utility function for transforming long data to SNP by sample matrix
# wrapper for reshape2::dcast which handles extra columns introduced by missing sites
# and cleans up the matrix by moving the row ID column to rownames.
get.variant.by.sample.matrix <- function(long.data, variant.id, value.var) {
    # transform to SNP by sample matrix
    variant.by.sample.matrix <- reshape2::dcast(
        data = long.data,
        formula = variant.id ~ Indiv,
        value.var = value.var
        );
    # remove row ID column
    rownames(variant.by.sample.matrix) <- variant.by.sample.matrix$variant.id;
    variant.by.sample.matrix$variant.id <- NULL;
    # remove column named NA if it exists (caused by missing sites)
    variant.by.sample.matrix <- variant.by.sample.matrix[, !colnames(variant.by.sample.matrix) %in% 'NA'];

    return(variant.by.sample.matrix);

    }
