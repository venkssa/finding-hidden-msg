
from ch1.part2 import *
from ch1.part2 import _neighbors
import unittest


class Part2Tests(unittest.TestCase):

    def test_skew_i(self):
        self.assertListEqual(skew_i("CATGGGCATCGGCCATACGCC"),
                             [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2])

    def test_min_skew(self):
        self.assertListEqual(min_skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"), [11, 24])

    def test_hamming_count(self):
        self.assertEqual(hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC"), 3)

    def test_approximate_pattern_match(self):
        self.assertListEqual(approximate_pattern_matching(
            "ATTCTGGA", "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT", 3),
            [6, 7, 26, 27])

    def test_count_2(self):
        self.assertEqual(count_2("AACAAGCTGATAAACATTTAAAGAG", "AAAAA"), 11)

    def test_approximate_pattern_count(self):
        self.assertEqual(
            approximate_pattern_count("GAGG", "TTTAGAGCCTTCAGAGG", 2),
            4)

    def test_neighbor(self):
        self.assertSetEqual(_neighbors("ACG", 1),
                            {"CCG", "TCG", "GCG", "AAG", "ATG", "AGG", "ACA", "ACC", "ACT", "ACG"})

        self.assertSetEqual(_neighbors("AAT", 0), {"AAT"})

    def test_frequent_words_with_mismatch(self):
        self.assertSetEqual(
            frequent_words_with_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1),
            {"GATG", "ATGC", "ATGT"})

    def test_frequent_words_with_mismatch_and_revese_complement(self):
        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1),
            {"ATGT", "ACAT"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("AAAAAAAAAA", 2, 1),
            {"AT", "TA"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("AGTCAGTC", 4, 2),
            {"AATT", "GGCC"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("AATTAATTGGTAGGTAGGTA", 4, 0),
            {"AATT"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("ATA", 3, 1),
            {"AAA", "AAT", "ACA", "AGA", "ATA", "ATC", "ATG", "ATT", "CAT", "CTA", "GAT", "GTA", "TAA", "TAC", "TAG",
             "TAT", "TCT", "TGT", "TTA", "TTT"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("AAT", 3, 0),
            {"AAT", "ATT"}
        )

        self.assertSetEqual(
            frequent_words_with_mismatch_and_reverse_complement("TAGCG", 2, 1),
            {"CA", "CC", "GG", "TG"}
        )
