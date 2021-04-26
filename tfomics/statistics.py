"""Statistical methods used for estimating the effect size at
potential allele-specific binding sites.
"""

import pandas as pd
from scipy.optimize import curve_fit


def calculate_effect_size(data_frame, param_column="ar", sterr_column="ar_sterr"):
    """Calculate the ASB effect size by translating and rescaling the estimated allelic ratio
    and associated standard deviations.
    """
    return data_frame.apply(
        lambda row: pd.Series(
            (1 - 2.0 * row[param_column], 2 * row[sterr_column]),
            index=["es", "es_sterr"],
        ),
        axis=1,
    )


def estimate_binomial_probability(
    data_frame, positives_column="ref_count", negatives_column="alt_count"
):
    """Estimate the allelic ratio from the read counts to a reference and alternate allele,
    treating them as positives and negatives in a binomial experiment.

    Returns a new data frame with the columns
    ar - allelic ratio
    ar_sterr - standard deviation for the allelic ratio
    """
    return data_frame.apply(
        lambda row: pd.Series(
            _binomial_probability_and_variance(
                row[positives_column], row[negatives_column]
            ),
            index=["ar", "ar_sterr"],
        ),
        axis=1,
    )


def _binomial_probability_and_variance(positives, negatives):
    """Return an estimate for the probability and
    variance of a single binomial experiment with the number of
    positive and negative outcomes as given.
    """
    # Require a minimum of 1 for both positive and negative outcomes.
    # This is equivalent to setting a bound on your estimator for
    # p by assuming that if you ran one more trial it would come up
    # with the opposite outcome.

    positives = max(positives, 1)
    negatives = max(negatives, 1)

    total = positives + negatives

    p_estimate = positives / total
    var_estimate = (p_estimate * (1.0 - p_estimate) / total) ** 0.5

    return p_estimate, var_estimate


def _pool_summary_statistics(points, sterr):
    """Pool summary statistics from multiple samples by performing
    a weighted least-squares fit of a constant to the data points.

    points - points to be averaged
    sterr - list of standard deviations for each point

    Returns a tuple with the pooled average of the input points, and
    the estimated standard deviation of the pooled average.
    """
    if len(points) == 1:
        return (points.iat[0], sterr.iat[0])

    params, cov = curve_fit(
        lambda x, a: a, range(len(points)), points, sigma=sterr, absolute_sigma=True
    )

    return (params[0], cov[0][0] ** 0.5)


def group_statistics(
    data_frame,
    param_column="ar",
    sterr_column="ar_sterr",
    group_columns=("CHROM", "POS"),
):
    """Combine effect sizes and variances across multiple samples.
    Inputs:
    data_frame - data frame containing summary statistics
    param_column - column containing the effect size to be pooled
    stedv_column - column containing the standard deviations of the effect size
    group_columns - column names to group by (e.g. chromosome and position for SNPs)

    Returns a new data frame indexed by `group_columns` with the pooled statistics.
    """
    return data_frame.groupby(list(group_columns)).apply(
        lambda group: pd.Series(
            _pool_summary_statistics(group[param_column], group[sterr_column]),
            index=[param_column, sterr_column],
        )
    )


def get_ref_and_alt_counts(row):
    """Get the read counts for the reference and alternate allele from
    an allele-seq dataframe.
    """
    ref = row["ref"]

    # This won't quite work if neither the maternal or paternal
    # alleles are the reference allele, though this is rare.
    if row["mat_all"] == row["ref"]:
        alt = row["pat_all"]
    else:
        alt = row["mat_all"]

    return pd.Series((row[f"c{ref}"], row[f"c{alt}"]), index=["ref_count", "alt_count"])


def allele_seq_effect_size(data_frame):
    """Calculate the ASB effect-size of a data frame containing allele-seq data.
    This method is designed to work on an alle-seq dataframe and presumes the
    following columns are present:
    chrm - chromosome
    snppos - position of heterozygous SNP
    mat_all/pat_all - maternal and paternal alleles
    cA/cC/cG/cT - read counts for each nucleotide at SNP site
    ref - nucleotide considered as the reference allele
    """
    # Check that all required columns are present.
    required_cols = {
        "chrm",
        "snppos",
        "mat_all",
        "pat_all",
        "cA",
        "cC",
        "cG",
        "cT",
        "ref",
    }
    missing_cols = required_cols - set(data_frame.columns)

    if missing_cols:
        raise KeyError("data frame missing columns: {}".format(missing_cols))

    # Compute reference and alternate counts

    probabilities = data_frame.apply(get_ref_and_alt_counts, axis=1).pipe(
        estimate_binomial_probability
    )
    return group_statistics(
        data_frame.join(probabilities), group_columns=("chrm", "snppos")
    ).pipe(calculate_effect_size)
