"""Unit tests for the mendelian randomisation module"""
import pytest
import tfomics.mendelian_randomisation as mr
import pandas as pd


@pytest.fixture
def snps():
    return pd.DataFrame(
        {
            "snp": {6: "rs3814316", 372: "rs7622475", 4285: "rs8105903"},
            "REF": {6: "G", 372: "T", 4285: "C"},
            "ALT": {6: "A", 372: "C", 4285: "A"},
            "es": {
                6: -0.3676815447570283,
                372: 0.8662327873562312,
                4285: -0.8788795836308416,
            },
            "es_sterr": {
                6: 0.0795660627120642,
                372: 0.0352281703029674,
                4285: 0.0489189501992082,
            },
        }
    ).reset_index(drop=True)


@pytest.fixture
def gwas():
    return pd.DataFrame(
        {
            "rsid": {
                0: "rs7622475",
                1: "rs7622475",
                2: "rs7622475",
                3: "rs7622475",
                4: "rs8105903",
                5: "rs8105903",
                6: "rs8105903",
                7: "rs8105903",
                8: "rs3814316",
                9: "rs3814316",
                10: "rs3814316",
                11: "rs3814316",
            },
            "allele": {
                0: "C",
                1: "C",
                2: "C",
                3: "C",
                4: "A",
                5: "A",
                6: "A",
                7: "A",
                8: "A",
                9: "A",
                10: "A",
                11: "A",
            },
            "iscore": {
                0: 0.994035,
                1: 0.994035,
                2: 0.994035,
                3: 0.994035,
                4: 0.983091,
                5: 0.983091,
                6: 0.983091,
                7: 0.983091,
                8: 0.991887,
                9: 0.991887,
                10: 0.991887,
                11: 0.991887,
            },
            "beta": {
                0: -0.47226,
                1: -1.1963,
                2: 14.333,
                3: -0.45951,
                4: -0.49113,
                5: 6.7073,
                6: -0.29062,
                7: -0.24232,
                8: -9.6309,
                9: 0.52655,
                10: 0.49062,
                11: 0.8951600000000001,
            },
            "NSE": {
                0: 0.067366,
                1: 0.13735999999999998,
                2: 1.6988,
                3: 0.06737,
                4: 0.11014000000000003,
                5: 1.3594,
                6: 0.053982,
                7: 0.053976,
                8: 1.4164,
                9: 0.056216,
                10: 0.05621,
                11: 0.11478,
            },
            "trait": {
                0: "23108-0.0",
                1: "23106-0.0",
                2: "23105-0.0",
                3: "23107-0.0",
                4: "23106-0.0",
                5: "23105-0.0",
                6: "23108-0.0",
                7: "23107-0.0",
                8: "23105-0.0",
                9: "23108-0.0",
                10: "23107-0.0",
                11: "23106-0.0",
            },
            "MAF": {
                0: 0.200793,
                1: 0.200793,
                2: 0.200793,
                3: 0.200793,
                4: 0.450407,
                5: 0.450407,
                6: 0.450407,
                7: 0.450407,
                8: 0.369195,
                9: 0.369195,
                10: 0.369195,
                11: 0.369195,
            },
            "HWE": {
                0: 0.4474,
                1: 0.4474,
                2: 0.4474,
                3: 0.4474,
                4: 0.7424,
                5: 0.7424,
                6: 0.7424,
                7: 0.7424,
                8: 0.2104,
                9: 0.2104,
                10: 0.2104,
                11: 0.2104,
            },
        }
    )


@pytest.fixture
def mr_result():
    return pd.DataFrame(
        {
            "snp": {
                0: "rs3814316",
                1: "rs3814316",
                2: "rs3814316",
                3: "rs3814316",
                4: "rs7622475",
                5: "rs7622475",
                6: "rs7622475",
                7: "rs7622475",
                8: "rs8105903",
                9: "rs8105903",
                10: "rs8105903",
                11: "rs8105903",
            },
            "REF": {
                0: "G",
                1: "G",
                2: "G",
                3: "G",
                4: "T",
                5: "T",
                6: "T",
                7: "T",
                8: "C",
                9: "C",
                10: "C",
                11: "C",
            },
            "ALT": {
                0: "A",
                1: "A",
                2: "A",
                3: "A",
                4: "C",
                5: "C",
                6: "C",
                7: "C",
                8: "A",
                9: "A",
                10: "A",
                11: "A",
            },
            "es": {
                0: -0.3676815447570283,
                1: -0.3676815447570283,
                2: -0.3676815447570283,
                3: -0.3676815447570283,
                4: 0.8662327873562312,
                5: 0.8662327873562312,
                6: 0.8662327873562312,
                7: 0.8662327873562312,
                8: -0.8788795836308416,
                9: -0.8788795836308416,
                10: -0.8788795836308416,
                11: -0.8788795836308416,
            },
            "es_sterr": {
                0: 0.0795660627120642,
                1: 0.0795660627120642,
                2: 0.0795660627120642,
                3: 0.0795660627120642,
                4: 0.0352281703029674,
                5: 0.0352281703029674,
                6: 0.0352281703029674,
                7: 0.0352281703029674,
                8: 0.0489189501992082,
                9: 0.0489189501992082,
                10: 0.0489189501992082,
                11: 0.0489189501992082,
            },
            "rsid": {
                0: "rs3814316",
                1: "rs3814316",
                2: "rs3814316",
                3: "rs3814316",
                4: "rs7622475",
                5: "rs7622475",
                6: "rs7622475",
                7: "rs7622475",
                8: "rs8105903",
                9: "rs8105903",
                10: "rs8105903",
                11: "rs8105903",
            },
            "allele": {
                0: "A",
                1: "A",
                2: "A",
                3: "A",
                4: "C",
                5: "C",
                6: "C",
                7: "C",
                8: "A",
                9: "A",
                10: "A",
                11: "A",
            },
            "iscore": {
                0: 0.991887,
                1: 0.991887,
                2: 0.991887,
                3: 0.991887,
                4: 0.994035,
                5: 0.994035,
                6: 0.994035,
                7: 0.994035,
                8: 0.983091,
                9: 0.983091,
                10: 0.983091,
                11: 0.983091,
            },
            "beta": {
                0: -9.6309,
                1: 0.52655,
                2: 0.49062,
                3: 0.8951600000000001,
                4: -0.47226,
                5: -1.1963,
                6: 14.333,
                7: -0.45951,
                8: -0.49113,
                9: 6.7073,
                10: -0.29062,
                11: -0.24232,
            },
            "NSE": {
                0: 1.4164,
                1: 0.056216,
                2: 0.05621,
                3: 0.11478,
                4: 0.067366,
                5: 0.13735999999999998,
                6: 1.6988,
                7: 0.06737,
                8: 0.11014000000000003,
                9: 1.3594,
                10: 0.053982,
                11: 0.053976,
            },
            "trait": {
                0: "23105-0.0",
                1: "23108-0.0",
                2: "23107-0.0",
                3: "23106-0.0",
                4: "23108-0.0",
                5: "23106-0.0",
                6: "23105-0.0",
                7: "23107-0.0",
                8: "23106-0.0",
                9: "23105-0.0",
                10: "23108-0.0",
                11: "23107-0.0",
            },
            "MAF": {
                0: 0.369195,
                1: 0.369195,
                2: 0.369195,
                3: 0.369195,
                4: 0.200793,
                5: 0.200793,
                6: 0.200793,
                7: 0.200793,
                8: 0.450407,
                9: 0.450407,
                10: 0.450407,
                11: 0.450407,
            },
            "HWE": {
                0: 0.2104,
                1: 0.2104,
                2: 0.2104,
                3: 0.2104,
                4: 0.4474,
                5: 0.4474,
                6: 0.4474,
                7: 0.4474,
                8: 0.7424,
                9: 0.7424,
                10: 0.7424,
                11: 0.7424,
            },
            "MR total causal effect": {
                0: 26.193590995610894,
                1: -1.4320816682489605,
                2: -1.3343612345955846,
                3: -2.4346068296453134,
                4: -0.5451883222307388,
                5: -1.381037542634635,
                6: 16.546360527110444,
                7: -0.5304694150430838,
                8: 0.5588137546340942,
                9: -7.631648436172215,
                10: 0.3306710104692453,
                11: 0.2757146764052974,
            },
            "MR se": {
                0: 6.85340574323195,
                1: 0.34556522834616554,
                2: 0.3267273497771018,
                3: 0.6123885844616807,
                4: 0.08086778177952504,
                5: 0.16822431251270717,
                6: 2.0733699442809868,
                7: 0.08071016598670662,
                8: 0.1291209346499008,
                9: 1.6040108128959152,
                10: 0.06411975471957629,
                11: 0.06330293028608873,
            },
            "z score": {
                0: 3.821981650725737,
                1: -4.144171782278948,
                2: -4.08402062302377,
                3: -3.975591464993508,
                4: -6.74172470461871,
                5: -8.209500291643728,
                6: 7.980418821421891,
                7: -6.572522910316587,
                8: 4.327832323621765,
                9: -4.757853485035974,
                10: 5.157084769201226,
                11: 4.355480467637809,
            },
            "p value": {
                0: 0.00013238354251153526,
                1: 3.410440300095233e-05,
                2: 4.4263099938862287e-05,
                3: 7.020455467429789e-05,
                4: 1.565174356732314e-11,
                5: 2.2211092625853213e-16,
                6: 1.4583758922178308e-15,
                7: 4.94697726749336e-11,
                8: 1.5058403614300026e-05,
                9: 1.9566242805186023e-06,
                10: 2.5082426838562196e-07,
                11: 1.3277535680564031e-05,
            },
            "effect_allele": {
                0: "alt",
                1: "alt",
                2: "alt",
                3: "alt",
                4: "alt",
                5: "alt",
                6: "alt",
                7: "alt",
                8: "alt",
                9: "alt",
                10: "alt",
                11: "alt",
            },
            "q values": {
                0: 0.00013238354251153526,
                1: 4.54725373346031e-05,
                2: 5.311571992663474e-05,
                3: 7.658678691741587e-05,
                4: 6.260697426929255e-11,
                5: 2.6653311151023855e-15,
                6: 8.750255353306986e-15,
                7: 1.4840931802480083e-10,
                8: 2.258760542145004e-05,
                9: 3.913248561037205e-06,
                10: 6.019782441254927e-07,
                11: 2.258760542145004e-05,
            },
        },
    )


def test_causal_effect_no_gwas_effect():
    """Test that the causal effect calculation with gwas_effect = 0
    returns expected values for a few different cases"""

    # Setting gwas_effect to zero means standard error is equal
    # to gwas_error/exposure_effect^2 and effect is zero
    assert (0.0, 2.0) == mr.causal_effect(1.0, 42.3, 0.0, 2.0)

    # Changing exposure_error does not change that.
    # FIXME: that is mathematically the case, but why?
    assert (0.0, 2.0) == mr.causal_effect(1.0, 22.1, 0.0, 2.0)

    # Error is inversely proportional to exposure_effect
    assert (0.0, 1.0) == mr.causal_effect(2.0, 22.1, 0.0, 2.0)

    # Error scales linearly with gwas_error
    assert (0.0, 4.0) == mr.causal_effect(1.0, 22.1, 0.0, 4.0)


def test_causal_effect_full():
    """Test causal effect with gwas_effect != 0"""
    assert (3.0, 5.0) == mr.causal_effect(1.0, 1.0, 3.0, 4.0)
    # FIXME: add some more here?


def test_singularity_if_exposure_is_zero():
    """If the exposure is zero, we hit a division by zero"""
    with pytest.raises(ZeroDivisionError):
        mr.causal_effect(0.0, 23.0, 3.12, 2.17)


def test_naive_mr(snps, gwas, mr_result):
    """Run a full end-to-end test of the naive MR analysis"""
    gwas_codes = list(gwas.trait)
    columns_BVs = {
        "rsid": "snp",
        "es": "es",
        "es_sterr": "es_sterr",
        "alt": "ALT",
        "ref": "REF",
    }

    columns_GWAS = {
        "rsid": "rsid",
        "trait": "trait",
        "MAF": "MAF",
        "HWE": "HWE",
        "iscore": "iscore",
        "allele": "allele",
        "beta": "beta",
        "NSE": "NSE",
    }

    mr_snps = mr.naive_mr_analysis(
        snps,
        gwas,
        gwas_codes,
        permute=False,
        columns_bvs=columns_BVs,
        columns_gwas=columns_GWAS,
    )

    pd.testing.assert_frame_equal(mr_snps, mr_result)


@pytest.fixture
def filter_frame():
    return pd.DataFrame(
        {
            "id": ["high", "low"],
            "MAF": [1.0, 0.05],
            "HWE": [1.0, 0.05],
            "trait": ["foo", "bar"],
            "iscore": [1.0, 0.05],
        }
    )


def test_filter_gwas_min_maf(filter_frame):

    # if min_MAF < 0.05, both columns should survive
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # Should be the same for min_MAF == 0.05
    filtered = mr.filter_gwas(filter_frame, 0.05, 0.04, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # 1. > min_maf > 0.05 results in one entry, that entry being the low one
    filtered = mr.filter_gwas(filter_frame, 0.25, 0.04, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame[filter_frame.id == "high"], filtered)


def test_filter_gwas_min_hwe(filter_frame):
    # if min_HWE < 0.05, both columns should survive
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # Should be the same for min_HWE == 0.05
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.05, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # 1. > min_HWE > 0.05 results in one entry, that entry being the low one
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.25, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame[filter_frame.id == "high"], filtered)


def test_filter_gwas_min_iscore(filter_frame):
    # if min_iscore < 0.05, both columns should survive
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 0.04, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # Should be the same for min_iscore == 0.05
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 0.05, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame, filtered)

    # 1. > min_iscore > 0.05 results in one entry, that entry being the low one
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 0.25, list(filter_frame.trait))
    pd.testing.assert_frame_equal(filter_frame[filter_frame.id == "high"], filtered)


def test_high_filter_empty_result(filter_frame):
    """Test that a high threshold returns an empty dataframe"""
    # MAF
    filtered = mr.filter_gwas(filter_frame, 1.2, 0.04, 0.04, list(filter_frame.trait))
    assert filtered.empty

    # HWE
    filtered = mr.filter_gwas(filter_frame, 0.04, 1.2, 0.04, list(filter_frame.trait))
    assert filtered.empty

    # iscore
    filtered = mr.filter_gwas(filter_frame, 0.04, 0.04, 1.2, list(filter_frame.trait))
    assert filtered.empty


def test_filter_traits(filter_frame):
    """Check filtering based on a list of traits"""
    traits = ["foo"]
    filtered = mr.filter_gwas(filter_frame, 0.01, 0.01, 0.01, traits)
    assert list(filtered.trait.unique()) == traits

    traits = ["bar"]
    filtered = mr.filter_gwas(filter_frame, 0.01, 0.01, 0.01, traits)
    assert list(filtered.trait.unique()) == traits

    traits = ["baz"]
    filtered = mr.filter_gwas(filter_frame, 0.01, 0.01, 0.01, traits)
    assert filtered.empty
