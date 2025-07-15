# Handling CRAN warnings for data.table syntax:
if (getRversion() >= '2.15.1') utils::globalVariables(c(
    'ID.vcf.unsplit', 'ID.vcf', 'POS.pgs', 'POS.vcf',
    'REF', 'REF.vcf', 'CHROM.pgs', 'CHROM.vcf',
    'merge.strategy', 'POS', 'CHROM', 'ID', 'ALT', '.'
    ));

#' @title Combine VCF with PGS
#' @description Match PGS SNPs to corresponding VCF information by genomic coordinates or rsID using a merge operation.
#' @param vcf.data A data.frame containing VCF data. Required columns: \code{CHROM, POS}.
#' @param pgs.weight.data A data.frame containing PGS data. Required columns: \code{CHROM, POS}.
#' @return A list containing a data.frame of merged VCF and PGS data and a data.frame of PGS SNPs missing from the VCF.
#'
#' A primary merge is first performed on chromosome and base pair coordinates. For SNPs that could not be matched in the first mergs, a second merge is attempted by rsID if available.
#' This action can account for short INDELs that can have coordinate mismatches between the PGS and VCF data.
#' The merge is a left outer join: all PGS SNPs are kept as rows even if they are missing from the VCF, and all VCF SNPs that are not a component of the PGS are dropped.
#' If no PGS SNPs are present in the VCF, the function will terminate with an error.
#' @examples
#' # Example VCF
#' vcf.path <- system.file(
#'     'extdata',
#'     'HG001_GIAB.vcf.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' vcf.import <- import.vcf(vcf.path);
#'
#' # Example pgs weight file
#' pgs.weight.path <- system.file(
#'     'extdata',
#'     'PGS000662_hmPOS_GRCh38.txt.gz',
#'     package = 'ApplyPolygenicScore',
#'     mustWork = TRUE
#'     );
#' pgs.import <- import.pgs.weight.file(pgs.weight.path);
#'
#' merge.data <- combine.vcf.with.pgs(
#'     vcf.data = vcf.import$dat,
#'     pgs.weight.data = pgs.import$pgs.weight.data
#'     );
#' @export
combine.vcf.with.pgs <- function(vcf.data, pgs.weight.data) {

    # check that inputs are data.frames
    if (!is.data.frame(vcf.data)) {
        stop('vcf.data must be a data.frame');
        }

    if (!is.data.frame(pgs.weight.data)) {
        stop('pgs.weight.data must be a data.frame');
        }

    data.table::setDT(vcf.data);
    data.table::setDT(pgs.weight.data);

    # check that inputs contain required columns for mergeing
    required.vcf.columns <- c('CHROM', 'POS');
    required.pgs.columns <- c('CHROM', 'POS');
    if (!all(required.vcf.columns %in% colnames(vcf.data))) {
        stop('vcf.data must contain columns named CHROM and POS');
        }
    if (!all(required.pgs.columns %in% colnames(pgs.weight.data))) {
        stop('pgs.weight.data must contain columns named CHROM and POS');
        }

    # Ensure POS columns are of a consistent type for the merge
    if (typeof(vcf.data$POS) != 'integer') {
      #warning('Converting VCF POS column from ', typeof(vcf.data$POS), ' to integer.')
      vcf.data[, POS := as.integer(POS)];
    }

    if (typeof(pgs.weight.data$POS) != 'integer') {
      # warning('Converting PGS POS column from ', typeof(pgs.weight.data$POS), ' to integer.')
      pgs.weight.data[, POS := as.integer(POS)];
    }

    # Add a check to ensure no NAs were introduced by the conversion
    if(any(is.na(vcf.data$POS)) | any(is.na(pgs.weight.data$POS))) {
        stop('Data type conversion of POS column introduced NA values. Please check your input files for non-numeric POS entries.')
    }

    # check for optional ID column
    rsid.available <- TRUE;
    if (!'ID' %in% colnames(vcf.data)) {
        rsid.available <- FALSE;
        warning('ID column not found in VCF data. Merging by rsID will not be possible.');
        }

    if (!'ID' %in% colnames(pgs.weight.data)) {
        rsid.available <- FALSE;
        warning('ID column not found in PGS weight data. Merging by rsID will not be possible.');
        }

    # match CHROM notation in pgs data to vcf data formatting
    chr.prefix <- grepl('^chr', vcf.data$CHROM[1]);
    numeric.sex.chr <- any(grepl('23$', vcf.data$CHROM));

    # pgs.weight.data$CHROM <- ApplyPolygenicScore::format.chromosome.notation(
    #     chromosome = pgs.weight.data$CHROM,
    #     chr.prefix = chr.prefix,
    #     numeric.sex.chr = numeric.sex.chr
    #     );

    pgs.weight.data[, CHROM := ApplyPolygenicScore::format.chromosome.notation(
        chromosome = CHROM,
        chr.prefix = chr.prefix,
        numeric.sex.chr = numeric.sex.chr
        )];

    # 'left outer join' on CHROM and POS columns
    # this merge keeps all pgs SNPs even if they are missing from the VCF
    # and drops all VCF SNPs that are not a component of the pgs

    # merged.vcf.with.pgs.data <- merge(
    #     x = pgs.weight.data,
    #     y = vcf.data,
    #     by = c('CHROM', 'POS'),
    #     suffixes = c('.pgs', '.vcf'),
    #     all.x = TRUE
    #     );

    # merged.vcf.with.pgs.data$merge.strategy <- 'genomic coordinate';

    # Left outer join on CHROM and POS
    # 'left outer join' on CHROM and POS columns using data.table::merge()
    merged.vcf.with.pgs.data <- data.table::merge.data.table(
        x = pgs.weight.data,
        y = vcf.data,
        by = c('CHROM', 'POS'),
        suffixes = c('.pgs', '.vcf'),
        all.x = TRUE
        );
    merged.vcf.with.pgs.data[, merge.strategy := 'genomic coordinate'];

    # check for pgs SNPs missing from the VCF data
    missing.snp.data <- NULL;
    missing.pgs.snp.index <- is.na(merged.vcf.with.pgs.data$REF);

    # attempt additional merge for missing SNPs, this time on rsID
    if (any(missing.pgs.snp.index) & rsid.available) {

        # Subset merged data corresponding to missing PGS SNPs (data.table does this by reference)
        missing.snp.merged.data <- merged.vcf.with.pgs.data[which(missing.pgs.snp.index)];

        # Expand the VCF ID column (HIGHLY OPTIMIZED)
        cols.to.keep <- setdiff(colnames(vcf.data), 'ID');
        # Expand the VCF ID column programmatically, preserving all other columns
        split.rsid.vcf.data <- vcf.data[, {
            # Store the original ID before unlisting
            original.id <- ID;
            # Unlist the split IDs
            split.ids <- unlist(strsplit(ID, ';', fixed = TRUE));
            # Return a list of the new ID and the preserved original ID
            list(ID = split.ids, ID.vcf.unsplit = original.id);
        }, by = cols.to.keep];

        # Drop NA-filled unmatched VCF columns from first merge
        cols.to.drop.from.missing.data <- setdiff(cols.to.keep, c('CHROM', 'POS'));
        missing.snp.merged.data[, (cols.to.drop.from.missing.data) := NULL]

        # Merge missing SNP data on rsID using data.table::merge()
        merged.vcf.with.missing.pgs.data <- data.table::merge.data.table(
            x = missing.snp.merged.data,
            y = split.rsid.vcf.data,
            by.x = 'ID.pgs',
            by.y = 'ID',
            suffixes = c('.pgs', '.vcf'),
            all.x = TRUE
            );

        # Update missing SNP index
        second.merge.missing.pgs.snp.index <- is.na(merged.vcf.with.missing.pgs.data$REF);

        if (all(missing.pgs.snp.index) & (sum(missing.pgs.snp.index) == sum(second.merge.missing.pgs.snp.index))) {
            stop('All PGS SNPs are missing from the VCF, terminating merge.');
        } else if (any(second.merge.missing.pgs.snp.index)) {
            warning(paste('PGS is missing', sum(second.merge.missing.pgs.snp.index)), ' SNPs from VCF');
            missing.snp.data <- merged.vcf.with.missing.pgs.data[which(second.merge.missing.pgs.snp.index)];
        } else {
            missing.snp.data <- NULL;
        }

        # Keep coordinates from VCF data for matched SNPs
        # Keep coordinates from PGS data for unmatched (missing) SNPs
        merged.vcf.with.missing.pgs.data[, `:=`(
            CHROM = fifelse(is.na(REF), CHROM.pgs, CHROM.vcf),
            POS = fifelse(is.na(REF), POS.pgs, POS.vcf)
            )];

        # Unconditionally update the ID.vcf column with the VCF's original unsplit ID
        merged.vcf.with.missing.pgs.data[, ID.vcf := ID.vcf.unsplit];

        # Add and rename columns to match original merge
        # add columns to match original merge
        merged.vcf.with.missing.pgs.data[, merge.strategy := 'rsID'];

        # subset columns to match original merge
        merged.vcf.with.missing.pgs.data <- merged.vcf.with.missing.pgs.data[, colnames(merged.vcf.with.pgs.data), with=FALSE];

        # combine merged data
        merged.vcf.with.pgs.data <- data.table::rbindlist(list(
            merged.vcf.with.pgs.data[which(!missing.pgs.snp.index)],
            merged.vcf.with.missing.pgs.data
            ));

        } else if (any(missing.pgs.snp.index)) {

        warning(paste('PGS is missing', sum(missing.pgs.snp.index)), ' SNPs from VCF');
        missing.snp.data <- merged.vcf.with.pgs.data[which(missing.pgs.snp.index)];

        } else {

        missing.snp.data <- NULL;

        }

    output <- list(
        merged.vcf.with.pgs.data = merged.vcf.with.pgs.data,
        missing.snp.data = missing.snp.data
        );

    return(output);
    }
