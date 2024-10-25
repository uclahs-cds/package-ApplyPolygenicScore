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
    'assess.pgs.vcf.allele.match correctly checks input', {

        # input class
        expect_error(
            assess.pgs.vcf.allele.match(1, 'A', 'A', 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.pgs.vcf.allele.match('A', 1, 'A', 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.pgs.vcf.allele.match('A', 'A', 1, 'T'),
            'alleles must be character vectors'
            );

        expect_error(
            assess.pgs.vcf.allele.match('A', 'A', 'T', 1),
            'alleles must be character vectors'
            );

        # input length
        expect_error(
            assess.pgs.vcf.allele.match('A', c('A', 'T'), 'A', 'TT'),
            'vcf.ref.allele and vcf.alt.allele must be the same length.'
            );

        expect_error(
            assess.pgs.vcf.allele.match('A', 'A', 'TT', c('A', 'T')),
            'pgs.ref.allele and pgs.effect.allele must be the same length.'
            );

        expect_error(
            assess.pgs.vcf.allele.match(c('A', 'A'), c('A', 'A'), 'A', 'TT'),
            'vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, and pgs.effect.allele must be the same length.'
            );

        # Accepted alleles
        expect_silent(
            assess.pgs.vcf.allele.match(
                c('A', '*', 'C', 'G', 'T'),
                c('T', '*', 'G', 'C', 'A'),
                c('A', 'T', 'C', 'G', 'T'),
                c('T', 'A', 'G', 'C', 'A')
                )
            );
        }
    );

test_that(
    'assess.pgs.vcf.allele.match correctly formats outputs', {
        test.output <- assess.pgs.vcf.allele.match(c('A', 'A'), c('C', 'G'), c('A', 'A'), c('C', 'G'));

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
    'assess.pgs.vcf.allele.match correctly handles SNP cases', {
        test.output <- assess.pgs.vcf.allele.match(vcf.ref.allele, vcf.alt.allele, pgs.ref.allele, pgs.effect.allele);

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


test_that(
    'assess.pgs.vcf.allele.match correctly handles INDEL cases', {
        # indel mismatch case
        expect_warning(
            assess.pgs.vcf.allele.match('A', 'G', 'A', 'ATCG'),
            'Mismatch detected in INDEL PGS allele. Skipping strand flip assessment.'
            );
        # indel mismatch case
        expect_warning(
            assess.pgs.vcf.allele.match('ATCG', 'A', 'G', 'A'),
            'Mismatch detected in INDEL VCF allele. Skipping strand flip assessment.'
            );
        # indel but everything matches
        expect_silent(
            assess.pgs.vcf.allele.match('ATCG', 'A', 'ATCG', 'A')
            );

        # check indel handling when return.indels.as.missing == FALSE
        # VCF ALT allele is an INDEL
        test.leave.indels.vcf.alt <- assess.pgs.vcf.allele.match(
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
        test.leave.indels.vcf.ref <- assess.pgs.vcf.allele.match(
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
        test.leave.indels.pgs.effect <- assess.pgs.vcf.allele.match(
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
        test.leave.indels.pgs.other <- assess.pgs.vcf.allele.match(
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
        test.missing.indels.vcf.alt <- assess.pgs.vcf.allele.match(
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
        test.missing.indels.vcf.ref <- assess.pgs.vcf.allele.match(
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
        test.missing.indels.pgs.effect <- assess.pgs.vcf.allele.match(
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
        test.missing.indels.pgs.other <- assess.pgs.vcf.allele.match(
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
    'assess.pgs.vcf.allele.match correctly handles ambiguous and unresolved cases', {

        # palindromic case, return.ambiguous.as.missing == TRUE

        test.ambiguous.missing <- assess.pgs.vcf.allele.match(
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
        test.unresolved.missing <- assess.pgs.vcf.allele.match(
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
    'assess.pgs.vcf.allele.match correctly handles multi-allelic cases', {
        # case with extra alleles that don't change the result
        vcf.alt.multi.allele <- c('G,T', 'G,T', 'G,T', 'G,T', 'G,T', 'G,T')
        test.output.benign.multi <- assess.pgs.vcf.allele.match(vcf.ref.allele, vcf.alt.multi.allele, pgs.ref.allele, pgs.effect.allele);

        expect_equal(
            test.output.benign.multi$match.status,
            c('default_match', 'effect_switch', 'strand_flip', 'effect_switch_with_strand_flip', 'ambiguous_flip', 'unresolved_mismatch')
            );

        expect_equal(
            test.output.benign.multi$new.pgs.effect.allele,
            c('G', 'A', 'G', 'A', 'A', 'C')
            );
        
        expect_equal(
            test.output.benign.multi$new.pgs.other.allele,
            c('A', 'G', 'A', 'G', 'T', 'A')
            );

        # case with extra alleles that can change unresolved result
        vcf.alt.multi.allele <- c('G,T', 'G,A', 'G,C', 'G,T', 'A,T', 'G,C')
        test.output.unresolved.multi <- assess.pgs.vcf.allele.match(vcf.ref.allele, vcf.alt.multi.allele, pgs.ref.allele, pgs.effect.allele);

        expect_equal(
            test.output.unresolved.multi$match.status,
            c('default_match', 'effect_switch', 'strand_flip', 'effect_switch_with_strand_flip', 'ambiguous_flip', 'default_match')
            );

        expect_equal(
            test.output.unresolved.multi$new.pgs.effect.allele,
            c('G', 'A', 'G', 'A', 'A', 'C')
            );

        expect_equal(
            test.output.unresolved.multi$new.pgs.other.allele,
            c('A', 'G', 'A', 'G', 'T', 'A')
            );

        # case with extra indels that don't change the result
        vcf.alt.multi.allele <- c('G,TGCA', 'G,TGCA', 'G,TGCA', 'G,TGCA', 'GTCA,T', 'G,TGCA')
        test.output.benign.multi.indel <- assess.pgs.vcf.allele.match(vcf.ref.allele, vcf.alt.multi.allele, pgs.ref.allele, pgs.effect.allele);

        expect_equal(
            test.output.benign.multi.indel$match.status,
            c('default_match', 'effect_switch', 'strand_flip', 'effect_switch_with_strand_flip', 'ambiguous_flip', 'unresolved_mismatch')
            );

        expect_equal(
            test.output.benign.multi.indel$new.pgs.effect.allele,
            c('G', 'A', 'G', 'A', 'A', 'C')
            );

        expect_equal(
            test.output.benign.multi.indel$new.pgs.other.allele,
            c('A', 'G', 'A', 'G', 'T', 'A')
            );

        # case with only indel multi-alleles in a mismatch
        expect_warning(
            assess.pgs.vcf.allele.match(
                c('A', 'A'),
                c('TG,TGCA', 'G,TGCA'),
                c('A', 'A'),
                c('G', 'TGCA')
                ),
            'Mismatch detected in INDEL VCF allele. Skipping strand flip assessment.'
            );

        # case with only indel multi-alleles in a match
        test.output.indel.multi <- assess.pgs.vcf.allele.match(
            c('A', 'A'),
            c('TG,TGCA', 'G,TGCA'),
            c('A', 'A'),
            c('TGCA', 'G')
            );

        expect_equal(
            test.output.indel.multi$match.status,
            c('default_match', 'default_match')
            );

        expect_equal(
            test.output.indel.multi$new.pgs.effect.allele,
            c('TGCA', 'G')
            );

        expect_equal(
            test.output.indel.multi$new.pgs.other.allele,
            c('A', 'A')
            );
        }
    );
