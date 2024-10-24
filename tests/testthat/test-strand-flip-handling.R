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

        flipped.alleles <- c('T', 'A', 'G', 'C', NA, '*');
        expect_equal(
            flip.DNA.allele(c('A', 'T', 'C', 'G', NA, '*')),
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


test_that(
    'assess.strand.flip correctly checks input', {

        # input class
        expect_error(
            assess.strand.flip(1, 'A', 'A', 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.strand.flip('A', 1, 'A', 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.strand.flip('A', 'A', 1, 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.strand.flip('A', 'A', 'T', 1),
            'alleles must be character vectors'
            );

        # input length
        expect_error(
            assess.strand.flip('A', c('A', 'T'), 'A', 'TT'),
            'vcf.ref.allele and vcf.alt.allele must be the same length.'
            );

        expect_error(
            assess.strand.flip('A', 'A', 'TT', c('A', 'T')),
            'pgs.ref.allele and pgs.effect.allele must be the same length.'
            );

        expect_error(
            assess.strand.flip(c('A', 'A'), c('A', 'A'), 'A', 'TT'),
            'vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, and pgs.effect.allele must be the same length.'
            );

        # Accepted alleles
        expect_silent(
            assess.strand.flip(
                c('A', '*', 'C', 'G', 'T'),
                c('T', '*', 'G', 'C', 'A'),
                c('A', 'T', 'C', 'G', 'T'),
                c('T', 'A', 'G', 'C', 'A')
                )
            );
        }
    );

test_that(
    'assess.strand.flip correctly formats outputs', {
        test.output <- assess.strand.flip(c('A', 'A'), c('C', 'G'), c('A', 'A'), c('C', 'G'));

        expect_equal(
            class(test.output),
            'list'
            );

        expect_equal(
            length(test.output),
            3
            );

        expect_equal(
            names(test.output),
            c('match.status', 'new.pgs.effect.allele', 'new.pgs.other.allele')
            );
        }
    );

# test data (single alleles):
# 1. no strand flips
# 2. effect allele switch
# 3. strand flip
# 4. effect allele switch AND strand flip
# 5. palindromic (ambiguous) alleles
# 6. unresolved mismatch

vcf.ref.allele <- c('A', 'A', 'A', 'A', 'A', 'A')
vcf.alt.allele <- c('G', 'G', 'G', 'G', 'T', 'G')
pgs.ref.allele <- c('A', 'G', 'T', 'C', 'T', 'A')
pgs.effect.allele <- c('G', 'A', 'C', 'T', 'A', 'C')

test_that(
    'assess.strand.flip correctly handles SNP cases', {
        test.output <- assess.strand.flip(vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, pgs.effect.allele);

        expect_equal(
            test.output$match.status,
            c('default_match', 'effect_switch', 'strand_flip', 'effect_switch_with_strand_flip', 'ambiguous_flip', 'unresolved_mismatch')
            );

        expect_equal(
            test.output$new.pgs.effect.allele,
            c('G', 'A', 'G', 'A', 'A', 'C')
            );

        expect_equal(
            test.output$new.pgs.other.allele,
            c('A', 'G', 'A', 'G', 'T', 'A')
            );
        }
    );

vcf.ref.allele.indel.test <- c('A', 'A', 'A', 'A', 'A')
vcf.alt.allele.indel.test <- c('G', 'G', 'G', 'G', 'T')
pgs.ref.allele.indel.test <- c('A', 'G', 'T', 'C', 'T')
pgs.effect.allele.indel.test <- c('G', 'A', 'C', 'T', 'A')

test_that(
    'assess.strand.flip correctly handles INDEL cases', {
        # indel mismatch case
        expect_warning(
            assess.strand.flip('A', 'G', 'A', 'ATCG'),
            'Mismatch detected in INDEL PGS allele. Skipping strand flip assessment.'
            );
        # indel mismatch case
        expect_warning(
            assess.strand.flip('ATCG', 'A', 'G', 'A'),
            'Mismatch detected in INDEL VCF allele. Skipping strand flip assessment.'
            );
        # indel but everything matches
        expect_silent(
            assess.strand.flip('ATCG', 'A', 'ATCG', 'A')
            );

        # check indel handling when return.indels.as.missing == FALSE
        # VCF ALT allele is an INDEL
        test.leave.indels.vcf.alt <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'ATCG'),
            c('A', 'A'),
            c('G', 'G'),
            return.indels.as.missing = FALSE
            );

        expect_equal(
            test.leave.indels.vcf.alt$match.status,
            c('default_match', 'indel_mismatch')
            );

        expect_equal(
            test.leave.indels.vcf.alt$new.pgs.effect.allele,
            c('G', 'G')
            );

        expect_equal(
            test.leave.indels.vcf.alt$new.pgs.other.allele,
            c('A', 'A')
            );

        # VCF REF allele is an INDEL
        test.leave.indels.vcf.ref <- assess.strand.flip(
            c('A', 'ATCG'),
            c('G', 'A'),
            c('A', 'A'),
            c('G', 'G'),
            return.indels.as.missing = FALSE
            );

        expect_equal(
            test.leave.indels.vcf.ref$match.status,
            c('default_match', 'indel_mismatch')
            );

        expect_equal(
            test.leave.indels.vcf.ref$new.pgs.effect.allele,
            c('G', 'G')
            );

        expect_equal(
            test.leave.indels.vcf.ref$new.pgs.other.allele,
            c('A', 'A')
            );

        # PGS effect allele is an INDEL
        test.leave.indels.pgs.effect <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'G'),
            c('A', 'A'),
            c('ATCG', 'G'),
            return.indels.as.missing = FALSE
            );

        expect_equal(
            test.leave.indels.pgs.effect$match.status,
            c('indel_mismatch', 'default_match')
            );

        expect_equal(
            test.leave.indels.pgs.effect$new.pgs.effect.allele,
            c('ATCG', 'G')
            );

        expect_equal(
            test.leave.indels.pgs.effect$new.pgs.other.allele,
            c('A', 'A')
            );

        # PGS other allele is an INDEL
        test.leave.indels.pgs.other <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'G'),
            c('ATCG', 'A'),
            c('G', 'G'),
            return.indels.as.missing = FALSE
            );

        expect_equal(
            test.leave.indels.pgs.other$match.status,
            c('indel_mismatch', 'default_match')
            );

        expect_equal(
            test.leave.indels.pgs.other$new.pgs.effect.allele,
            c('G', 'G')
            );
        
        expect_equal(
            test.leave.indels.pgs.other$new.pgs.other.allele,
            c('ATCG', 'A')
            );

        # check indel handling when return.indels.as.missing == TRUE
        # VCF ALT allele is an INDEL
        test.missing.indels.vcf.alt <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'ATCG'),
            c('A', 'A'),
            c('G', 'G'),
            return.indels.as.missing = TRUE
            );

        expect_equal(
            test.missing.indels.vcf.alt$match.status,
            c('default_match', 'indel_mismatch')
            );

        expect_equal(
            test.missing.indels.vcf.alt$new.pgs.effect.allele,
            c('G', NA)
            );

        expect_equal(
            test.missing.indels.vcf.alt$new.pgs.other.allele,
            c('A', NA)
            );

        # VCF REF allele is an INDEL
        test.missing.indels.vcf.ref <- assess.strand.flip(
            c('A', 'ATCG'),
            c('G', 'A'),
            c('A', 'A'),
            c('G', 'G'),
            return.indels.as.missing = TRUE
            );

        expect_equal(
            test.missing.indels.vcf.ref$match.status,
            c('default_match', 'indel_mismatch')
            );

        expect_equal(
            test.missing.indels.vcf.ref$new.pgs.effect.allele,
            c('G', NA)
            );

        expect_equal(
            test.missing.indels.vcf.ref$new.pgs.other.allele,
            c('A', NA)
            );

        # PGS effect allele is an INDEL
        test.missing.indels.pgs.effect <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'G'),
            c('A', 'A'),
            c('ATCG', 'G'),
            return.indels.as.missing = TRUE
            );

        expect_equal(
            test.missing.indels.pgs.effect$match.status,
            c('indel_mismatch', 'default_match')
            );

        expect_equal(
            test.missing.indels.pgs.effect$new.pgs.effect.allele,
            c(NA, 'G')
            );

        expect_equal(
            test.missing.indels.pgs.effect$new.pgs.other.allele,
            c(NA, 'A')
            );

        # PGS other allele is an INDEL
        test.missing.indels.pgs.other <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'G'),
            c('ATCG', 'A'),
            c('G', 'G'),
            return.indels.as.missing = TRUE
            );

        expect_equal(
            test.missing.indels.pgs.other$match.status,
            c('indel_mismatch', 'default_match')
            );

        expect_equal(
            test.missing.indels.pgs.other$new.pgs.effect.allele,
            c(NA, 'G')
            );

        expect_equal(
            test.missing.indels.pgs.other$new.pgs.other.allele,
            c(NA, 'A')
            );
        }

    );

test_that(
    'assess.strand.flip correctly handles ambiguous and unresolved cases', {

        # palindromic case, return.ambiguous.as.missing == TRUE

        test.ambiguous.missing <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'T'),
            c('A', 'T'),
            c('G', 'A'),
            return.ambiguous.as.missing = TRUE
            );

        expect_equal(
            test.ambiguous.missing$match.status,
            c('default_match', 'ambiguous_flip')
            );

        expect_equal(
            test.ambiguous.missing$new.pgs.effect.allele,
            c('G', NA)
            );

        expect_equal(
            test.ambiguous.missing$new.pgs.other.allele,
            c('A', NA)
            );

        # unresolved case, return.ambiguous.as.missing == TRUE
        test.unresolved.missing <- assess.strand.flip(
            c('A', 'A'),
            c('G', 'T'),
            c('A', 'T'),
            c('G', 'T'),
            return.ambiguous.as.missing = TRUE
            );

        expect_equal(
            test.unresolved.missing$match.status,
            c('default_match', 'unresolved_mismatch')
            );

        expect_equal(
            test.unresolved.missing$new.pgs.effect.allele,
            c('G', NA)
            );

        expect_equal(
            test.unresolved.missing$new.pgs.other.allele,
            c('A', NA)
            );
        }
    );

test_that(
    
)

assess.strand.flip(
    vcf.ref.allele = vcf.ref.allele,
    vcf.alt.allele = vcf.alt.allele,
    pgs.ref.allele = pgs.ref.allele,
    pgs.effect.allele = pgs.effect.allele,
    return.indels.as.missing = FALSE,
    return.ambiguous.as.missing = FALSE
    );

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