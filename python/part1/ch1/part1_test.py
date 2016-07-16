
import unittest

from ch1.part1 import *
from ch1.part1 import _pattern_count, _computing_frequencies, number_to_pattern, pattern_to_number


class Part1Tests(unittest.TestCase):

    def test_valid_pattern_count(self):
        self.assertEqual(_pattern_count("ACAACTATGCATACTATCGGGAACTATCCT", "ACTAT"), 3)
        self.assertEqual(_pattern_count("GCGCG", "GCG"), 2)

    def test_no_valid_pattern_count(self):
        self.assertEqual(_pattern_count("ACAT", "GCT"), 0)

    def test_frequent_words(self):
        self.assertSetEqual(frequent_words("ACAACTATGCATACTATCGGGAACTATCCT", 5), {"ACTAT"})
        self.assertSetEqual(frequent_words("CGATATATCCATAG", 3), {"ATA"})
        self.assertSetEqual(frequent_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4), {"CATG", "GCAT"})

    def test_pattern_to_number(self):
        self.assertEqual(pattern_to_number("AAAA"), 0)
        self.assertEqual(pattern_to_number("AAAC"), 1)
        self.assertEqual(pattern_to_number("AAGA"), 8)

    def test_number_to_pattern(self):
        self.assertEqual(number_to_pattern(8, 4), "AAGA")
        self.assertEqual(number_to_pattern(0, 4), "AAAA")
        self.assertEqual(number_to_pattern(1, 4), "AAAC")

        self.assertEqual(number_to_pattern(8, 2), "GA")
        self.assertEqual(number_to_pattern(0, 2), "AA")
        self.assertEqual(number_to_pattern(1, 2), "AC")

    def test_computing_frequencies(self):
        self.assertListEqual(_computing_frequencies("ACGCGGCTCTGAAA", 2),
                             [2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0])

    def test_faster_frequent_words(self):
        self.assertSetEqual(faster_frequent_words("ACAACTATGCATACTATCGGGAACTATCCT", 5), {"ACTAT"})
        self.assertSetEqual(faster_frequent_words("CGATATATCCATAG", 3), {"ATA"})
        self.assertSetEqual(faster_frequent_words("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4), {"CATG", "GCAT"})

    def test_faster_frequent_words_by_sorting(self):
        self.assertSetEqual(faster_frequent_words_by_sorting("ACAACTATGCATACTATCGGGAACTATCCT", 5), {"ACTAT"})
        self.assertSetEqual(faster_frequent_words_by_sorting("CGATATATCCATAG", 3), {"ATA"})
        self.assertSetEqual(faster_frequent_words_by_sorting("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4), {"CATG", "GCAT"})

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AAAACCCGGT"), "ACCGGGTTTT")

    def test_pattern_matching(self):
        self.assertListEqual(pattern_matching("ATAT", "GATATATGCATATACTT"), [1, 3, 9])

    def test_find_clumps_slow(self):
        self.assertSetEqual(
            find_clumps_slow("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", k=5, el=50, t=4),
            {"CGACA", "GAAGA"})

    def test_find_clumps(self):
        self.assertSetEqual(
            find_clumps("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", k=5, el=50, t=4),
            {"CGACA", "GAAGA"})
