"""Tests for the Altschul and Erickson dinucleotide shuffle implementation"""
import random
import re

import pytest

from tfomics import dinucleotide_shuffle


def get_all_dinucleotide_pairs(sequence):
    """Get all the dinucleotide pairs in a string"""
    captures = re.finditer(r'(?=(\["ACTG"]{2}))', sequence)
    return [capture.group(1) for capture in captures]


def test_preserves_pair_frequency():
    """Sample 100 shuffles of a short snippet of dinucleotides and assert
    that all sequences have the same frequency of dinucleotide pairs
    """
    sequence = ("AGCAGAAGCAGGATACAGGGCAGCTCTGAGGCAAGGTAGGC"
                + "AGGTGCTGTGGTGCTCCCAGGTAGCCTAGTGGGATGCAGGAG")

    samples = [dinucleotide_shuffle(sequence) for _ in range(100)]

    sequence_dinucleotides = get_all_dinucleotide_pairs(sequence).sort()
    for sample in samples:
        sample_dinucleotides = get_all_dinucleotide_pairs(sample).sort()
        assert sample_dinucleotides == sequence_dinucleotides


def test_error_invalid_nucleotide():
    """Check that we validate the nucleotides by passing in a sequence
    containing an X
    """
    sequence = ("AGCAGAAGCAGGATACAGGGCAGCTCTGAGGCAAGGTAGGC"
                + "AGGTGCTGTGGTGCTCCCAGGTAGCCTXGTGGGATGCAGGAG")

    with pytest.raises(AssertionError):
        dinucleotide_shuffle(sequence)


def test_can_handle_mixed_case():
    """Pass in a string with mixed case. The result will all be uppercase,
    but the shuffle will run successfully.
    """
    sequence = ("TGCTTACTGGCtaattATTGGttaaggTATTTACTGATTGTCACTT" +
                "ATTATTggttaaggtATTtactGATTGtcactTACAGGGGTTAGCA")

    shuffled = dinucleotide_shuffle(sequence)
    assert shuffled == shuffled.upper()


def test_deterministic_once_seeded():
    """Redo the same set of shuffles twice, re-seeding the RNG with the
    same seed each time. This should produce the same results twice.
    """
    sequence = ("AGCAGAAGCAGGATACAGGGCAGCTCTGAGGCAAGGTAGGC"
                + "AGGTGCTGTGGTGCTCCCAGGTAGCCTCGTGGGATGCAGGAG")

    rng = random.Random(42)
    first_set = [dinucleotide_shuffle(sequence, rng=rng) for _ in range(10)]

    rng = random.Random(42)
    second_set = [dinucleotide_shuffle(sequence, rng=rng) for _ in range(10)]
    assert first_set == second_set
