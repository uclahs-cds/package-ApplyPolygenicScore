test_that(
    'flip.DNA.allele correctly checks input', {

        expect_error(
            flip.DNA.allele(1),
            'allele must be a character vector'
            );
        
        expect_error(
            flip.DNA.allele(c('A', 'T', 'C', 'G', 1)),
            'Invalid allele: 1'
            );

        expect_error(
            flip.DNA.allele(c('A', 'T', 'C', 'G', 'N')),
            'Invalid allele: N'
            );

        expect_error(
            flip.DNA.allele('N'),
            'Invalid allele: N'
            );
        }

    );

test_that(
    'flip.DNA.allele correctly flips SNP alleles', {

        flipped.alleles <- c('T', 'A', 'G', 'C', NA);
        names(flipped.alleles) <- c('A', 'T', 'C', 'G', NA);
        expect_equal(
            flip.DNA.allele(c('A', 'T', 'C', 'G', NA)),
            flipped.alleles
            );
        }
    );
