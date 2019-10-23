"""Statistical methods used for estimating the effect size at
potential allele-specific binding sites.
"""

from scipy.optimize import curve_fit


def _fit_multiple_binomial_samples(ref_counts, alt_counts):
    """Fit multiple ASB samples by doing a least-squares fit of
    a curve with points (0, Ref/Total) and (1, Alt/Total).

    Returns an estimated effect-size (Ref-Alt)/Total and
    an estimated variance in the form of sqrt(Cov[0,0]),
    wher Cov is the covariance matrix returned by the least-squares
    fit.
    """
    x_coords = []
    y_coords = []
    sigma = []
    for ref, alt in zip(ref_counts, alt_counts):
        probability, variance = _binomial_probability_and_variance(ref, alt)

        x_coords.append(0.)
        x_coords.append(1.)

        y_coords.append(1.-probability)
        y_coords.append(probability)

        sigma.append(variance)
        sigma.append(variance)

    params, cov = curve_fit(lambda x, a, b: (a*x)+b,
                            x_coords,
                            y_coords,
                            sigma=sigma,
                            absolute_sigma=True)

    return params[0], cov[0][0]**0.5


def _binomial_effect_size(ref_count, alt_count):
    """Treat ref vs. alt as a binomial distribution and
    rescale/translate the AR to a range of [-1, 1]
    """
    prob, var = _binomial_probability_and_variance(ref_count, alt_count)

    return 2.*(prob-0.5), 2.*var


def _binomial_probability_and_variance(ref, alt):
    """Return an estimate for the probability and
    variance of a single binomial experiment where
    ref and alt are treated as the two possible outcomes.
    """
    total = alt+ref

    # If a SNP with 0 reads has been passed to us, throw an error.
    if total == 0:
        raise ValueError("Encountered SNP with 0 reads")

    p_estimate = alt/total
    var_estimate = (p_estimate*(1.-p_estimate)/total)**0.5

    return p_estimate, var_estimate


def asb_effect_size(ref_counts, alt_counts):
    """Estimate the effect size and variance at an allele-specific binding site
    by treating each sample as an independent sample following a binomial distribution.

    If we have multiple samples, we perform a least-squares fit to a line with coordinates
    (0, Ref/Total) and (1, Alt/Total) and return the slope of the line and the
    square root of the covariance of the slope.
    """

    if len(ref_counts) == 1:
        return _binomial_effect_size(ref_counts[0], alt_counts[0])
    else:
        return _fit_multiple_binomial_samples(ref_counts, alt_counts)
