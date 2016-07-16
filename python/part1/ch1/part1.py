
import math


def _pattern_count(text, pattern):
    counts = (1 for idx in range(0, len(text) - len(pattern) + 1) if text[idx:idx + len(pattern)] == pattern)
    return sum(counts)


def frequent_words(text, k):
    kmers = set(text[idx:idx + k] for idx in range(0, len(text) - k + 1))
    kmer_count = [(kmer, _pattern_count(text, kmer)) for kmer in kmers]
    max_count = max(kmer_count, key=lambda ck: ck[1])[1]
    return set((kmer for kmer, count in kmer_count if count == max_count))


_nbase = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_rbase = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def pattern_to_number(pattern):
    value = 0

    for ch in pattern:
        value = value * 4 + _nbase[ch]

    return value


def number_to_pattern(number, k):
    pattern = ['A' for _ in range(k)]

    idx = k - 1
    while number > 0:
        pattern[idx] = _rbase[number % 4]
        number = math.floor(number / 4)
        idx -= 1

    return "".join(pattern)


def _computing_frequencies(text, k):
    freq = [0 for _ in range(0, 4 ** k)]

    for kmer in (text[idx:idx + k] for idx in range(0, len(text) - k + 1)):
        num = pattern_to_number(kmer)
        freq[num] += 1

    return freq


def faster_frequent_words(text, k):
    freq = _computing_frequencies(text, k)
    max_count = max(freq)
    return set((number_to_pattern(num, k) for num, count in enumerate(freq) if count == max_count))


def faster_frequent_words_by_sorting(text, k):
    encoded_kmers = sorted([pattern_to_number(text[idx:idx + k]) for idx in range(0, len(text) - k + 1)])
    counts = [1 for _ in range(0, len(encoded_kmers))]

    for first, second in zip(enumerate(encoded_kmers), enumerate(encoded_kmers[1:], 1)):
        if first[1] == second[1]:
            counts[second[0]] += counts[first[0]]

    max_count = max(counts)
    return set((number_to_pattern(kmer, k) for kmer, count in zip(encoded_kmers, counts) if count == max_count))


_complement = str.maketrans('ATGC', 'TACG')


def reverse_complement(pattern):
    return pattern.translate(_complement)[::-1]


def pattern_matching(pattern, genome):
    """
     Find all occurrences of a pattern in a string.

    :return: All starting positions in Genome where Pattern appears as a substring.
    """
    occurrences = []
    idx = genome.find(pattern)
    while idx != -1:
        occurrences.append(idx)
        idx = genome.find(pattern, idx + 1)
    return occurrences


def find_clumps_slow(genome, k, el, t):
    """
    Find patterns forming clumps in a string.
    :param genome:
    :param k: length of k-mer
    :param el: length of window
    :param t: min num of appearances for selection
    :return: All distinct k-mers forming (el, t)-clumps in Genome.
    """
    genome_slices = (genome[idx:idx + el] for idx in range(0, len(genome) - el))
    frequencies = (_computing_frequencies(genome_slice, k) for genome_slice in genome_slices)
    return set((number_to_pattern(number, k) for frequency in frequencies
                for (number, freq_count) in enumerate(frequency) if freq_count >= t))


class _FreqArray(object):

    def __init__(self, t):
        self._freq_count = {}
        self._clumps = set()
        self._t = t

    def dec_count(self, key):
        count = self._freq_count.get(key)
        if count is None or count == 0:
            raise KeyError("Cannot decrement {} as it does not exist in freq_count".format(key))
        if count == self._t:
            self._clumps.remove(key)
        if count == 1:
            self._freq_count.pop(key)
        else:
            self._freq_count[key] = count - 1

    def inc_count(self, key):
        self._freq_count[key] = self._freq_count.get(key, 0) + 1
        if self._freq_count[key] == self._t:
            self._clumps.add(key)

    def count(self, key):
        return self._freq_count.get(key, 0)

    def kmers_in_clumps(self):
        return self._clumps


def _genome_to_numeric_genome_convertor(genome, k):
    numeric_genome = []
    value = 0
    for ch in genome[0:k]:
        value = value * 4 + _nbase[ch]
    numeric_genome.append(value)

    p = 4 ** (k - 1)
    for prev, next in zip(genome[0:len(genome) - k], genome[k:]):
        value -= (p * _nbase[prev])
        value = value * 4 + _nbase[next]
        numeric_genome.append(value)

    return numeric_genome


def find_clumps(genome, k, el, t):
    """
    Find patterns forming clumps in a string.
    :param genome:
    :param k: length of k-mer
    :param el: length of window
    :param t: min num of appearances for selection
    :return: All distinct k-mers forming (el, t)-clumps in Genome.
    """
    numeric_genome = _genome_to_numeric_genome_convertor(genome, k)
    
    nel = el - k + 1

    freq_count = _FreqArray(t)

    for number in numeric_genome[0:el]:
        freq_count.inc_count(number)

    clumps = set(freq_count.kmers_in_clumps())

    for prev, next in zip(numeric_genome[0: len(numeric_genome) - nel], numeric_genome[nel:]):
        freq_count.dec_count(prev)
        freq_count.inc_count(next)
        clumps.update(freq_count.kmers_in_clumps())

    return set((number_to_pattern(kmer, k) for kmer in clumps))
