"""Tests for the statistics module"""
import random

import pandas as pd
import pytest
import tfomics.statistics


def test_estimate_p():
    """Check a few values for a binomial distribution"""
    coin_flips = pd.DataFrame({"heads": [1, 9, 15, 12],
                               "tails": [9, 1, 15, 33]})
    probs = tfomics.statistics.estimate_binomial_probability(coin_flips,
                                                             positives_column="heads",
                                                             negatives_column="tails")
    assert probs["ar"].round(5).tolist() == [0.1, 0.9, 0.5, 0.26667]
    assert probs["ar_sterr"].round(5).tolist() == [0.09487, 0.09487, 0.09129, 0.06592]


def test_add_one_if_no_reads():
    """Check that we add a read if there's a zero in either column before estimating p"""
    coin_flips = pd.DataFrame({"heads": [1, 0, 9],
                               "tails": [9, 9, 0]})
    probs = tfomics.statistics.estimate_binomial_probability(coin_flips,
                                                             positives_column="heads",
                                                             negatives_column="tails")
    assert probs["ar"].round(5).tolist() == [0.1, 0.1, 0.9]
    assert probs["ar_sterr"].round(5).tolist() == [0.09487, 0.09487, 0.09487]


def test_swap_positive_and_negative():
    """If we swap success for failure in a binomial, we should have p->1-p"""
    coin_flips = pd.DataFrame({"heads": [random.randint(0, 100) for _ in range(100)],
                               "tails": [random.randint(0, 100) for _ in range(100)]})

    heads = tfomics.statistics.estimate_binomial_probability(coin_flips,
                                                             positives_column="heads",
                                                             negatives_column="tails")
    tails = tfomics.statistics.estimate_binomial_probability(coin_flips,
                                                             positives_column="tails",
                                                             negatives_column="heads")
    assert all(heads["ar"].round(5) == (1-tails["ar"]).round(5))
    assert all(heads["ar_sterr"].round(5) == tails["ar_sterr"].round(5))


def test_rescale_p_and_variance():
    """Check a few sample values are rescaled correclty by calculate_effect_size"""
    ratios = pd.DataFrame({"ar": [0.1, 0.9, 0.5, 0.26667],
                           "ar_sterr": [0.09487, 0.09487, 0.09129, 0.06592]})
    effect_sizes = tfomics.statistics.calculate_effect_size(ratios)

    # Not the greatest test, but this is at least a regression check. These are
    # manually calculated.
    assert effect_sizes["es"].round(4).tolist() == [0.8, -0.8, 0., 0.4667]
    assert effect_sizes["es_sterr"].round(4).tolist() == [0.1897, 0.1897, 0.1826, 0.1318]


def test_group_reffects():
    """"Check that grouping of effects and standard deviations behaves as expected
    at a few sample points
    """
    data = pd.DataFrame({
        "group": ["A", "A", "B", "B", "B"],
        "effect": [1., 0., 1., 1., 0.],
        "sterr": [1., 1., 1., 1., 1.]
    })
    grouped_results = tfomics.statistics.group_statistics(data,
                                                          param_column="effect",
                                                          sterr_column="sterr",
                                                          group_columns=("group",))
    assert grouped_results["effect"].round(4).tolist() == [0.5, 0.6667]
    assert grouped_results["sterr"].round(4).tolist() == [0.7071, 0.5774]


def test_allele_seq_missing_columns():
    """We throw an error if we call an allele-seq specific method with a dataframe that doesn't
    have allele-seq's column structure
    """
    bad_data_frame = pd.DataFrame({"Foo": [1, 2, 3]})
    with pytest.raises(KeyError):
        tfomics.statistics.allele_seq_effect_size(bad_data_frame)


def test_full_allele_seq():
    """Run the full allele seq pipeline to estimate effect sizes and pool by SNP and position"""
    allele_seq_data = pd.DataFrame({'chrm': ["chr1", "chr1", "chr2", "chr2", "chr2"],
                                    'snppos': [42, 42, 55, 55, 150],
                                    'mat_all': ["C", "A", "G", "G", "T"],
                                    'pat_all': ["A", "C", "C", "C", "G"],
                                    'cA': [10, 10, 0, 0, 0],
                                    'cC': [2, 2, 5, 5, 0],
                                    'cG': [0, 0, 10, 10, 12],
                                    'cT': [0, 0, 0, 0, 8],
                                    'ref': ["C", "C", "G", "G", "T"]})
    results = tfomics.statistics.allele_seq_effect_size(allele_seq_data)
    # I'm assuming these come out in order of chromosome and SNP
    assert results["es"].round(4).tolist() == [0.6667, -0.3333, 0.2]
    assert results["es_sterr"].round(4).tolist() == [0.1521, 0.1721, 0.2191]
