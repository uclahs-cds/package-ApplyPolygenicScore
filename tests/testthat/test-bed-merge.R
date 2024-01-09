test_that(
    'merge.pgs.bed correctly checks input', {
        # check that pgs.bed.list is a named list
        expect_error(
            merge.pgs.bed(pgs.bed.list = 'not.a.named.list'),
            'pgs.bed.list must be a named list'
            );
        
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    ),
            'pgs.bed.list must be a named list'
            ))
        )

        # check that pgs.bed.list is a list of data.frames
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                test = 'not.a.data.frame'
                )),
            'all elements of pgs.bed.list must be data.frames'
            );

        # check that all elements of pgs.bed.list have the same column names
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4),
                    foo = c(1, 2, 3)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    )
                )),
            'all elements of pgs.bed.list must have the same column names'
            );

        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4),
                    foo = c(1, 2, 3)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4),
                    bar = c(1, 2, 3)
                    )
                )),
            'all elements of pgs.bed.list must have the same column names'
            );


        # check that pgs.bed.list has required columns
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    foo = c(1, 2, 3)
                    )
                )),
            'all elements of pgs.bed.list must have columns named chr, start, and end'
            );

        # check that provided coordinates are one bp in length
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(3, 4, 5)
                    )
                )),
            'all intervals specified in pgs.bed.list must represent one SNP and be one bp in length'
            );

        # check that annotation.column.index is whitin the bounds of the data
        expect_error(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    )
                ),
                add.annotation.data = TRUE,
                annotation.column.index = 5
                ),
            'annotation.column.index must be within the range of the number of columns in the data.frames in pgs.bed.list'
            );

        # check that correct input is accepted
        expect_silent(
            merge.pgs.bed(pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4),
                    foo = c(1, 2, 3)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4),
                    foo = c(1, 2, 3)
                    )
                ))
            );
        }
    );

test_that(
    'merge.pgs.bed correctly formats output', {

        simple.test.output <- merge.pgs.bed(pgs.bed.list = list(
            name1 = data.frame(
                chr = c('1', '2', '3'),
                start = c(1, 2, 3),
                end = c(2, 3, 4)
                ),
            name2 = data.frame(
                chr = c('1', '2', '3'),
                start = c(1, 2, 3),
                end = c(2, 3, 4)
                )
            ));

        # check that the output is a data.frame
        expect_equal(
            class(simple.test.output),
            'data.frame'
            );
        
        # check that the output has the correct number of columns
        expect_equal(
            ncol(simple.test.output),
            4
            );

        # check that the output has the correct column names
        expect_equal(
            colnames(simple.test.output),
            c('chr', 'start', 'end', 'annotation')
            );
        
    }
);

test_that(
    'merge.pgs.bed correctly merges overlapping intervals', {
        load('data/tiny.bed.test.data.Rda')

        simple.test.output <- merge.pgs.bed(
            pgs.bed.list = tiny.bed.test.data,
            add.annotation.data = TRUE,
            annotation.column.index = 4
            );
        
        simple.test.expected.output <- data.frame(
            chr = c('chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr3'),
            start = c(1, 2, 5, 3, 4, 4),
            end = c(2, 3, 6, 4, 5, 5),
            annotation = c(
                'no overlap with pgs2|pgs1',
                'overlap with pgs2|pgs1,overlap with pgs1|pgs2',
                'no overlap with pgs1|pgs2',
                'overlap with pgs2|pgs1,overlap with pgs1|pgs2',
                'no overlap with pgs2|pgs1',
                'no overlap with pgs1|pgs2'
                )
            );
        
        # check that the output matches expected output
        # comparison is peformed column by column to avoid errors due to rowname differences
        # caused by sorting
        expect_equal(
            dim(simple.test.output),
            dim(simple.test.expected.output)
            );

        expect_equal(
            simple.test.output$chr,
            simple.test.expected.output$chr
            );

        expect_equal(
            simple.test.output$start,
            simple.test.expected.output$start
            );

        expect_equal(
            simple.test.output$end,
            simple.test.expected.output$end
            );
        
        expect_equal(
            simple.test.output$annotation,
            simple.test.expected.output$annotation
            );

        }
    );

test_that(
    'merge.pgs.bed correctly adds slop', {

        small.slop.test.output <- merge.pgs.bed(
            pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    )
                ),
            slop = 1
            );
        
        expect_equal(
            small.slop.test.output$start,
            c(0, 1, 2)
            );
        
        expect_equal(
            small.slop.test.output$end,
            c(3, 4, 5)
            );

        large.slop.test.output <- merge.pgs.bed(
            pgs.bed.list = list(
                name1 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
                    )
                ),
            slop = 10
            );
        
        expect_equal(
            large.slop.test.output$start,
            c(0, 0, 0)
            );
        expect_equal(
            large.slop.test.output$end,
            c(12, 13, 14)
            );

        expect_warning(
            merge.pgs.bed(
                pgs.bed.list = list(
                    name1 = data.frame(
                        chr = c('1', '2', '3'),
                        start = c(1, 2, 3),
                        end = c(2, 3, 4)
                        ),
                    name2 = data.frame(
                        chr = c('1', '2', '3'),
                        start = c(1, 2, 3),
                        end = c(2, 3, 4)
                        )
                    ),
                slop = 10
                ),
            'Slop caused negative start coordinates.'
            );

        }
    );

test_that(
    'merge.pgs.bed works on real data', {
        pgs1 <- import.pgs.weight.file('data/PGS000662_hmPOS_GRCh38.txt');
        pgs1.bed <- convert.pgs.to.bed(pgs.weight.data = pgs1$pgs.weight.data);
        pgs1.bed <- subset(pgs1.bed, select = c('chr', 'start', 'end'));
        pgs2 <- import.pgs.weight.file('data/PGS003378_hmPOS_GRCh38.txt');
        pgs2.bed <- convert.pgs.to.bed(pgs.weight.data = pgs2$pgs.weight.data);
        pgs2.bed <- subset(pgs2.bed, select = c('chr', 'start', 'end'));
        pgs.bed.list <- list(
            PGS000662 = pgs1.bed,
            PGS003378 = pgs2.bed
            );
        
        expect_silent(
            merge.pgs.bed(pgs.bed.list = pgs.bed.list)
            );
        }
        
    );
