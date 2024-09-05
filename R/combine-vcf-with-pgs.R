#' @title Combine VCF with PGS
#' @description Match PGS SNPs to corresponding VCF information by genomic coordinates using a merge operation.
#' @param vcf.data A data.frame containing VCF data. Required columns: \code{CHROM, POS}.
#' @param pgs.weight.data A data.frame containing PGS data. Required columns: \code{CHROM, POS}.
#' @return A list containing a data.frame of merged VCF and PGS data and a data.frame of PGS SNPs missing from the VCF.
#'
#' Merge is performed on chromosome and base pair coordinates.
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

    # check that inputs contain required columns for mergeing
    required.vcf.columns <- c('CHROM', 'POS', 'ID');
    required.pgs.columns <- c('CHROM', 'POS', 'ID');
    if (!all(required.vcf.columns %in% colnames(vcf.data))) {
        stop('vcf.data must contain columns named CHROM, POS, and ID');
        }
    if (!all(required.pgs.columns %in% colnames(pgs.weight.data))) {
        stop('pgs.weight.data must contain columns named CHROM, POS, and ID');
        }

    # match CHROM notation in pgs data to vcf data formatting
    chr.prefix <- grepl('^chr', vcf.data$CHROM[1]);
    numeric.sex.chr <- any(grepl('23$', vcf.data$CHROM));

    pgs.weight.data$CHROM <- ApplyPolygenicScore::format.chromosome.notation(
        chromosome = pgs.weight.data$CHROM,
        chr.prefix = chr.prefix,
        numeric.sex.chr = numeric.sex.chr
        );

    # "left outer join" on CHROM and POS columns
    # this merge keeps all pgs SNPs even if they are missing from the VCF
    # and drops all VCF SNPs that are not a component of the pgs
    merged.vcf.with.pgs.data <- merge(
        x = pgs.weight.data,
        y = vcf.data,
        by = c('CHROM', 'POS'),
        suffixes = c('.pgs', '.vcf'),
        all.x = TRUE
        );

    merged.vcf.with.pgs.data$merge.strategy <- 'genomic coordinate';

    # check for pgs SNPs missing from the VCF data
    missing.snp.data <- NULL;
    missing.pgs.snp.index <- is.na(merged.vcf.with.pgs.data$REF);

    # attempt additional merge for missing SNPs, this time on rsID
    if (any(missing.pgs.snp.index)) {
        # extract merged rows corresponding to missing PGS SNPs
        missing.snp.merged.data <- merged.vcf.with.pgs.data[missing.pgs.snp.index, ];
        # rename ID column to match pre-merge pgs weight data
        missing.snp.merged.data$ID <- missing.snp.merged.data$ID.pgs;
        missing.snp.merged.data$ID.vcf <- NULL;
        missing.snp.merged.data$ID.pgs <- NULL;
        # extract only original pgs weight columns from merged missing SNP data (for neater second merge)
        missing.snp.pgs.weight.data <- subset(missing.snp.merged.data, select = colnames(pgs.weight.data));
        rm(missing.snp.merged.data);

        # Split VCF$ID column into separate rows for each rsID (multiple rsIDs are separated by ;)
        # most efficient way to do this is to use the data.table package
        if (any(grepl(';', vcf.data$ID))) {
            data.table::setDT(vcf.data);
            split.rsid.vcf.data <- merge(
                x = vcf.data,
                # split only entries with multiple rsIDs, save in new column, and merge back with the original data
                y = vcf.data[grepl(';', ID), strsplit(as.character(ID), ";", fixed=TRUE), by = .(Indiv, CHROM, POS)
                    ][,.(new.ID = V1, Indiv, CHROM, POS)],
                by = c('CHROM', 'POS', 'Indiv'),
                all = TRUE
                );
            # replace entries with multiple rsIDs with the new, split, rsID
            split.rsid.vcf.data <- split.rsid.vcf.data[!is.na(new.ID), ID := new.ID][, new.ID := NULL];
            } else {
            split.rsid.vcf.data <- vcf.data;
            }

        # merge missing SNP data on split rsID
        merged.vcf.with.missing.pgs.data <- merge(
            x = missing.snp.pgs.weight.data,
            y = split.rsid.vcf.data,
            suffixes = c('.pgs', '.vcf'),
            by = 'ID',
            all.x = TRUE
            );
        rm(split.rsid.vcf.data);

        # update missing SNP index
        second.merge.missing.pgs.snp.index <- is.na(merged.vcf.with.missing.pgs.data$REF);

        if (all(missing.pgs.snp.index) & (sum(missing.pgs.snp.index) == sum(second.merge.missing.pgs.snp.index))) {
            stop('All PGS SNPs are missing from the VCF, terminating merge.');
            } else if (any(second.merge.missing.pgs.snp.index)) {
                warning(paste('PGS is missing', sum(second.merge.missing.pgs.snp.index)), ' SNPs from VCF');
                missing.snp.data <- merged.vcf.with.missing.pgs.data[second.merge.missing.pgs.snp.index, ];
            } else {
            missing.snp.data <- NULL;
            }

        # keep coordinates from VCF data for matched SNPs with coordinate mismatch
        merged.vcf.with.missing.pgs.data[!is.na(merged.vcf.with.missing.pgs.data$REF), 'CHROM'] <- merged.vcf.with.missing.pgs.data[!is.na(merged.vcf.with.missing.pgs.data$REF), 'CHROM.vcf'];
        merged.vcf.with.missing.pgs.data[!is.na(merged.vcf.with.missing.pgs.data$REF), 'POS'] <- merged.vcf.with.missing.pgs.data[!is.na(merged.vcf.with.missing.pgs.data$REF), 'POS.vcf'];

        # keep coordinates from PGS data for unmatched (missing) SNPs
        merged.vcf.with.missing.pgs.data[is.na(merged.vcf.with.missing.pgs.data$REF), 'CHROM'] <- merged.vcf.with.missing.pgs.data[is.na(merged.vcf.with.missing.pgs.data$REF), 'CHROM.pgs'];
        merged.vcf.with.missing.pgs.data[is.na(merged.vcf.with.missing.pgs.data$REF), 'POS'] <- merged.vcf.with.missing.pgs.data[is.na(merged.vcf.with.missing.pgs.data$REF), 'POS.pgs'];

        # add columns to match original merge
        merged.vcf.with.missing.pgs.data$ID.pgs <- merged.vcf.with.missing.pgs.data$ID;
        merged.vcf.with.missing.pgs.data$ID.vcf <- merged.vcf.with.missing.pgs.data$ID;
        merged.vcf.with.missing.pgs.data$merge.strategy <- 'rsID';

        # subset columns to match original merge
        merged.vcf.with.missing.pgs.data <- subset(merged.vcf.with.missing.pgs.data, select = colnames(merged.vcf.with.pgs.data));

        # combine merged data
        merged.vcf.with.pgs.data <- rbind(merged.vcf.with.pgs.data[!missing.pgs.snp.index, ], merged.vcf.with.missing.pgs.data);
        rm(merged.vcf.with.missing.pgs.data);
        }

    output <- list(
        merged.vcf.with.pgs.data = merged.vcf.with.pgs.data,
        missing.snp.data = missing.snp.data
        );

    return(output);
    }
