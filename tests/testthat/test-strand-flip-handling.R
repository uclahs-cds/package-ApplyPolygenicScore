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
        expect_equal(
            flip.DNA.allele(c('A', 'T', 'C', 'G', NA)),
            flipped.alleles
            );
        }
    );

test_that(
    'flip.DNA.allele correctly flips INDEL alleles', {

        expect_equal(
            flip.DNA.allele(c('ATCG', 'TATA', 'AAAA', 'TTTA', 'AACCTTGGAACCTTGG', NA)),
            c('CGAT', 'TATA', 'TTTT', 'TAAA', 'CCAAGGTTCCAAGGTT', NA)
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

# test data (multi-allelic):
# 1. no strand flips
# 2. effect allele switch
# 3. strand flip
# 4. palindromic (ambiguous) alleles

vcf.ref.allele <- c('A', 'A', 'A', 'A')
vcf.alt.allele <- c('G,T', 'G,T', 'G,T', 'T,C')
pgs.ref.allele <- c('A', 'T', 'C', 'T')
pgs.effect.allele <- c('G', 'C', 'T', 'A')