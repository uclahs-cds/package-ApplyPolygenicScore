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
                    end = c(2, 3, 4)
                    ),
                name2 = data.frame(
                    chr = c('1', '2', '3'),
                    start = c(1, 2, 3),
                    end = c(2, 3, 4)
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
)