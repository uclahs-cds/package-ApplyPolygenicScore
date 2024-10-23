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
    'flip.DNA.allele correctly formats output', {

        expect_equal(
            class(flip.DNA.allele(c('A', 'T', 'C', 'G', NA))),
            'character'
            );
        }
    );

test_that(
    'flip.DNA.allele correctly flips SNP alleles', {

        flipped.alleles <- c('T', 'A', 'G', 'C', NA);
        expect_equal(
            flip.DNA.allele(c('A', 'T', 'C', 'G', NA)),
            flipped.alleles
            );
        }
    );

test_that(
    'flip.DNA.allele correctly handles INDEL alleles', {

        expect_equal(
            flip.DNA.allele(c('ATCG', 'TATA', 'AAAA', 'TTTA', 'AACCTTGGAACCTTGG', NA)),
            c('ATCG', 'TATA', 'AAAA', 'TTTA', 'AACCTTGGAACCTTGG', NA) # no flip expected
            );

        expect_equal(
            flip.DNA.allele(c('ATCG', 'TATA', 'AAAA', 'TTTA', 'AACCTTGGAACCTTGG', NA), return.indels.as.missing = TRUE),
            c(NA, NA, NA, NA, NA, NA) # missing expected
            );
        }
    );

# test data (single alleles):
# 1. no strand flips
# 2. effect allele switch
# 3. strand flip
# 4. palindromic (ambiguous) alleles

vcf.ref.allele <- c('A', 'A', 'A', 'A')
vcf.alt.allele <- c('G', 'G', 'G', 'T')
pgs.ref.allele <- c('A', 'T', 'C', 'T')
pgs.effect.allele <- c('G', 'C', 'T', 'A')

# test data (indels):
# 1. no strand flips insertion
# 2. no strand flips deletion
# 3. effect allele switch insertion
# 4. effect allele switch deletion
# 3. strand flip insertion
# 4. strand flip deletion
# strand flip and effect allele switch insertion

vcf.ref.indel.allele <- c('A', 'ATCG', 'A', 'ATCG', 'A', 'ATCG', 'A');
vcf.alt.indel.allele <- c('ATCG', 'A', 'ATCG', 'A', 'ATCG', 'A', 'ATCG');
pgs.ref.indel.allele <- c('A', 'ATCG', 'ATCG', 'A', 'G', 'GCGA', 'CGAT');
pgs.alt.indel.allele <- c('ATCG', 'A', 'A', 'ATCG', 'GCGA', 'G', 'A');


# test data (multi-allelic):
# 1. no strand flips
# 2. effect allele switch
# 3. strand flip
# 4. palindromic (ambiguous) alleles

vcf.ref.allele <- c('A', 'A', 'A', 'A')
vcf.alt.allele <- c('G,T', 'G,T', 'G,T', 'T,C')
pgs.ref.allele <- c('A', 'T', 'C', 'T')
pgs.effect.allele <- c('G', 'C', 'T', 'A')