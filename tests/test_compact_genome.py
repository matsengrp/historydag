import historydag.compact_genome as compact_genome
from historydag.compact_genome import CompactGenome
import historydag.parsimony_utils
from historydag.parsimony_utils import default_nt_transitions


def test_sequence_cg_convert():
    seqs = [
        "AAAA",
        "TAAT",
        "CTGA",
        "TGCA",
    ]
    for refseq in seqs:
        for seq in seqs:
            cg = compact_genome.compact_genome_from_sequence(seq, refseq)
            reseq = cg.to_sequence()
            if reseq != seq:
                print("\nUnmatched reconstructed sequence:")
                print("ref sequence:", refseq)
                print("sequence:", seq)
                print("cg:", cg)
                print("reconstructed sequence:", reseq)
                assert False


def test_cg_diff():
    refseq = "C" * 1000
    cgs = [
        CompactGenome({287: ("C", "G")}, refseq),
        CompactGenome({287: ("C", "G"), 318: ("C", "A"), 495: ("C", "T")}, refseq),
        CompactGenome(
            {287: ("C", "G"), 80: ("C", "T"), 257: ("C", "G"), 591: ("C", "A")}, refseq
        ),
        CompactGenome(
            {
                287: ("C", "G"),
                191: ("C", "G"),
                492: ("C", "G"),
                612: ("C", "G"),
                654: ("C", "G"),
            },
            refseq,
        ),
        CompactGenome({287: ("C", "G"), 318: ("C", "A"), 495: ("C", "T")}, refseq),
    ]
    # First try with apply_muts
    for parent_cg in cgs:
        for child_cg in cgs:
            assert (
                parent_cg.apply_muts(
                    str(p) + str(key) + str(c)
                    for p, c, key in compact_genome.cg_diff(parent_cg, child_cg)
                )
                == child_cg
            )
            assert all(
                par_nuc != child_nuc
                for par_nuc, child_nuc, idx in compact_genome.cg_diff(
                    parent_cg, child_cg
                )
            )
    # Now with apply_muts_raw
    for parent_cg in cgs:
        for child_cg in cgs:
            assert (
                parent_cg.apply_muts_raw(
                    (key, p, c)
                    for p, c, key in compact_genome.cg_diff(parent_cg, child_cg)
                )
                == child_cg
            )
            assert all(
                par_nuc != child_nuc
                for par_nuc, child_nuc, idx in compact_genome.cg_diff(
                    parent_cg, child_cg
                )
            )


def test_ambiguous_cg_distance():
    reference_seq = "AAAAAAAN"

    s1 = compact_genome.CompactGenome({1: ("A", "T"), 2: ("A", "G")}, reference_seq)
    s2 = compact_genome.CompactGenome({1: ("A", "T"), 2: ("A", "N")}, reference_seq)
    s3 = compact_genome.CompactGenome({1: ("A", "T")}, reference_seq)
    s4 = compact_genome.CompactGenome({8: ("N", "C")}, reference_seq)

    assert default_nt_transitions.min_weighted_cg_hamming_distance(s1, s2) == 0
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s1, s3) == 1
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s1, s4) == 3
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s4, s1) == 2
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s3, s2) == 0
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s4, s2) == 1
    assert default_nt_transitions.min_weighted_cg_hamming_distance(s4, s3) == 1


def test_reconcile_cgs():
    # ambiguitymap = historydag.parsimony_utils.standard_nt_ambiguity_map_gap_as_char
    ambiguitymap = historydag.parsimony_utils.standard_nt_ambiguity_map
    reference_seq = "AAAAAA-N"

    s1 = compact_genome.CompactGenome({1: ("A", "T"), 2: ("A", "G")}, reference_seq)
    s2 = compact_genome.CompactGenome({1: ("A", "C"), 2: ("A", "G")}, reference_seq)
    s3 = compact_genome.CompactGenome({1: ("A", "T")}, reference_seq)
    s4 = compact_genome.CompactGenome({8: ("N", "C")}, reference_seq)
    s5 = compact_genome.CompactGenome({1: ("A", "G"), 2: ("A", "T")}, reference_seq)

    a1 = compact_genome.CompactGenome({1: ("A", "Y"), 2: ("A", "G")}, reference_seq)
    a2 = compact_genome.CompactGenome({1: ("A", "Y"), 2: ("A", "R")}, reference_seq)
    a3 = compact_genome.CompactGenome({1: ("A", "B"), 2: ("A", "K")}, reference_seq)
    a4 = compact_genome.CompactGenome({1: ("A", "R"), 2: ("A", "W")}, reference_seq)
    a5 = compact_genome.CompactGenome({1: ("A", "D"), 2: ("A", "W")}, reference_seq)

    seqs = [s1, s2, s3, s4, s5]

    for s in seqs:
        assert s == compact_genome.reconcile_cgs((s,), ambiguitymap=ambiguitymap)[0]

    testcases = [
        ((s1, s2), a1),
        ((s1, s2, s3), a2),
        ((s1, s2, s5), a3),
        ((s4, s5), a4),
        ((s3, s4, s5), a5),
        ((s1, s1), s1),
    ]

    for cgs, answer in testcases:
        assert compact_genome.reconcile_cgs(cgs, ambiguitymap=ambiguitymap)[0] == answer


def test_ambiguous_cg_diff():
    ref = "AAAAAA"
    c1 = compact_genome.CompactGenome({}, ref)
    c2 = compact_genome.CompactGenome(
        {1: ("A", "C"), 2: ("A", "N"), 4: ("A", "K")}, ref
    )
    muts = set(compact_genome.ambiguous_cg_diff(c1, c2))
    assert muts == {("A", "C", 1), ("A", "G", 4)} or muts == {
        ("A", "C", 1),
        ("A", "T", 4),
    }
