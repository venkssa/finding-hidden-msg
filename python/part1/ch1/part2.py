
from ch1.part1 import reverse_complement
import itertools


_skew_value = {'A': 0, 'C': -1, 'G': 1, 'T': 0}


def skew_i(genome):
    return list(itertools.accumulate(itertools.chain([0], (_skew_value[x] for x in genome))))


def min_skew(genome):
    """
    Find a position in a genome where the skew diagram attains a minimum.

    :param genome:
    :return: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
    """
    skew = skew_i(genome)
    smallest_skew = min(skew)
    return [idx for idx, x in enumerate(skew) if x == smallest_skew]


def hamming_distance(genome1, genome2):
    """
    Compute the Hamming distance between two strings.

    :param genome1:
    :param genome2:
    :return: The Hamming distance between these strings.
    """
    return sum((1 for first, second in zip(genome1, genome2) if first != second))


def approximate_pattern_matching(pattern, genome, d):
    """
    Find all approximate occurrences of a pattern in a string.

    :param pattern:
    :param genome:
    :param  d:
    :return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.
    """
    plen = len(pattern)
    return [idx for idx in range(0, len(genome) - plen + 1)
            if hamming_distance(pattern, genome[idx: idx + plen]) <= d]


def count_2(genome, pattern):
    """
    Total number of occurrences of Pattern in Text with at most d mismatches

    :param genome:
    :param pattern:
    """
    return len(approximate_pattern_matching(pattern, genome, 2))


def approximate_pattern_count(pattern, genome, d):
    return len(approximate_pattern_matching(pattern, genome, d))


def _neighbors(pattern, d):
    """
     Finds the d-neighborhood of a string.
    :param pattern:
    :param d:
    :return: The collection of strings Neighbors(Pattern, d).
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set(_skew_value.keys())

    suffix_pattern = pattern[1:]
    neighbors = set()
    for text in _neighbors(suffix_pattern, d):
        if hamming_distance(text, suffix_pattern) < d:
            v = [x + text for x in _skew_value.keys()]
            neighbors.update(v)
        else:
            neighbors.add(pattern[0] + text)

    return neighbors


def frequent_words_with_mismatch(genome, k, d):
    """
    Find the most frequent k-mers with mismatches in a string.

    :param genome:
    :param k:
    :param d:
    :return:
    """
    kmers = {genome[idx:idx + k] for idx in range(0, len(genome) - k + 1)}
    patterns = set(itertools.chain(*(_neighbors(kmer, d) for kmer in kmers)))
    pattern_to_count = [(pattern, approximate_pattern_count(pattern, genome, d)) for pattern in patterns]
    max_count = max(pattern_to_count, key=lambda pc: pc[1])[1]
    return {pattern for pattern, count in pattern_to_count if count == max_count}


def _fast_approx_pattern_count(pattern, kmers, d):
    count = sum((1 for kmer in kmers if hamming_distance(pattern, kmer) <= d))
    return count


def frequent_words_with_mismatch_and_reverse_complement(genome, k, d):
    """
    Find the most frequent k-mers (with mismatches and reverse complements) in a string.
    :param genome:
    :param k:
    :param d:
    :return: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Pattern)
      over all possible k-mers.
    """
    kmers_set = {genome[idx:idx + k] for idx in range(0, len(genome) - k + 1)}
    kmers = [genome[idx:idx + k] for idx in range(0, len(genome) - k + 1)]
    patterns = set(itertools.chain(*(_neighbors(kmer, d) for kmer in kmers_set)))

    pattern_to_count = [
        (pattern, rc_pattern,
         _fast_approx_pattern_count(pattern, kmers, d) + _fast_approx_pattern_count(rc_pattern, kmers, d))
        for pattern, rc_pattern in zip(patterns, (reverse_complement(p) for p in patterns))]

    max_count = max(pattern_to_count, key=lambda pc: pc[2])[2]

    return set(itertools.chain(*([pattern, rc_pattern]
                                 for pattern, rc_pattern, count in pattern_to_count if count == max_count)))
