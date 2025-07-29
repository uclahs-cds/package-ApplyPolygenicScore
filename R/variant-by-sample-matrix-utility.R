# Handling CRAN warnings for data.table syntax:
if (getRversion() >= '2.15.1') utils::globalVariables(c('Indiv', 'N'));

type_stable_na <- function(x) {
  if (is.integer(x)) return(NA_integer_);
  if (is.numeric(x)) return(NA_real_);
  if (is.character(x)) return(NA_character_);
  if (is.logical(x)) return(NA);
  stop("Unhandled type in NA fallback");
}

custom.na.aggregation.function <- function(x) {
  if (length(x) > 1) {
    if (all(is.na(x))) {
      return(type_stable_na(x));
    } else {
      stop("Duplicate entries found");
    }
  } else {
    return(x[1]);
  }
}

custom.na.aggregation.function <- function(x) {
    if (length(x) > 1) {
        if (all(is.na(x))) {
            return(type_stable_na(x));
            } else {
                stop('Duplicate variant/effect-allele/sample combinations detected');
                }
        } else {
            return((x[1]));
            }
    }

# utility function for transforming long data to variant by sample matrix by user-provided variant ID (row ID) columns
# wrapper for data.table::dcast which handles extra columns introduced by missing sites
# and cleans up the matrix by moving the row ID column to rownames.
# Expects an "Indiv" column in the long data which contains sample identifiers.
get.variant.by.sample.matrix <- function(long.data, row.id.cols, value.col) {

    # transform to SNP by sample matrix
    setDT(long.data);

    # Get all unique row identifiers from the original data
    all.row.ids <- unique(long.data[, c(row.id.cols), with = FALSE]);

    # Filter the data to get only valid sample entries for pivoting
    filtered.data <- long.data[!is.na(Indiv) & !is.na(get(row.id.cols[1]))];

    # Perform validation to ensure no duplicates before casting
    check.cols <- c(row.id.cols, 'Indiv');
    duplicate.check <- filtered.data[, .N, by = check.cols];

    if (any(duplicate.check$N > 1)) {
        # Identify all combinations with duplicates
        duplicate.combos <- duplicate.check[N > 1];

        # Recreate the original key string format from the duplicate rows
        key.values.list <- as.list(duplicate.combos[, check.cols, with = FALSE]);
        duplicate.keys <- do.call(paste, c(key.values.list, sep = '_'));

        # Stop if duplicate variant entries (aka rows) are present.
        stop(
            paste('Duplicate variant/sample combinations detected:\n',
                paste(duplicate.keys, collapse = '\n'),
                '\nPlease ensure that each sample has only one genotype call for each variant id combination.\n'
                )
            );
        }

    # Construct the dcast formula dynamically based on user-provided columns
    formula.string <- paste(
        paste(row.id.cols, collapse = ' + '),
        '~',
        'Indiv'
        );

    # Pivot the data using dcast
    wide.data.dt <- data.table::dcast(
        data = filtered.data,
        formula = as.formula(formula.string),
        value.var = value.col
        );

    # Merge the pivoted data with all original row IDs to ensure all variants are present
    final.wide.dt <- merge(all.row.ids, wide.data.dt, by = row.id.cols, all.x = TRUE);

    # Convert the resulting data.table to a matrix
    # Create a character vector of all column names in the final data.table
    all.cols <- names(final.wide.dt);

    # Identify the columns that contain the data (not the row identifiers)
    data.cols.names <- setdiff(all.cols, row.id.cols);

    # Use the column names vector for robust subsetting
    data.cols.dt <- final.wide.dt[, data.cols.names, with = FALSE];

    # Get the row names vector directly from the merged data.table
    row.names.vector <- final.wide.dt[, do.call(paste, c(.SD, sep = ':')), .SDcols = row.id.cols];

    # Convert the data columns data.table to a matrix
    matrix.data <- as.matrix(data.cols.dt);

    # Assign the row names to the matrix
    rownames(matrix.data) <- row.names.vector;

    return(matrix.data);

    }

convert.long.vcf.to.wide.vcf <- function(long.vcf) {
    long.vcf <- as.data.frame(long.vcf);
    fixed.colnames <- c('CHROM', 'POS', 'ID', 'REF', 'ALT');
    fixed.data <- long.vcf[, fixed.colnames];
    fixed.data$variant.id <- paste0(fixed.data$CHROM, ':', fixed.data$POS);
    fixed.data <- unique(fixed.data); # remove duplicates
    # convert to data.table
    fixed.data <- data.table::as.data.table(fixed.data);

    allele.matrix <- get.variant.by.sample.matrix(
        long.data = long.vcf,
        row.id.cols = c('CHROM', 'POS'),
        value.col = 'gt_GT_alleles'
        );

    # sort allele matrix in order of corresponding fixed data variant id
    allele.matrix <- allele.matrix[match(fixed.data$variant.id, rownames(allele.matrix)), ];

    # save allele matrix row indices in fixed data
    fixed.data$allele.matrix.row.index <- seq_len(nrow(allele.matrix));

    # format output according to import.vcf standards
    output <- list(
        genotyped.alleles = allele.matrix,
        vcf.fixed.fields = fixed.data
        );

    return(output);

    }

custom.mean.dosage.aggregation.function <- function(x) {
    if (length(x) > 1) {
        if (all(is.na(x))) {
            return(as.numeric(NA));
            } else {
                return(mean(x, na.rm = TRUE));
                }
        } else {
            return(as.numeric(x[1]));
            }
    }

# utility function for transforming long data to SNP by sample matrix
# which handles multiallelic sites by averaging values or returning NA
get.combined.multiallelic.variant.by.sample.matrix <- function(long.data, variant.id, value.var) {
    setDT(long.data);

    # transform to SNP by sample matrix
    variant.by.sample.dt <- data.table::dcast(
        data = long.data,
        formula = variant.id ~ Indiv,
        value.var = value.var,
        fun.aggregate = custom.mean.dosage.aggregation.function
        );

    # Remove NA-named sample column if it exists
    if ('NA' %in% names(variant.by.sample.dt)) {
        variant.by.sample.dt[, 'NA' := NULL];
        }

    # remove row ID colum
    row.ids <- variant.by.sample.dt$variant.id;

    # Pull numeric matrix efficiently
    variant.by.sample.matrix <- as.matrix(variant.by.sample.dt[, -1, with = FALSE]);

    rownames(variant.by.sample.matrix) <- row.ids;

    return(variant.by.sample.matrix);
    }
