"""Unit tests for the mendelian randomisation module"""
import pytest
import tfomics.mendelian_randomisation as mr


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
