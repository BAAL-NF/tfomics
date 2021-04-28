"""MR analysis function for combinging AlleleSeq or BaalChIP data with GeneAtlas GWAS results"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests


def filter_gwas(
    effect,
    min_MAF,
    min_HWE,
    min_iscore,
    trait_list,
    columns={"MAF": "MAF", "HWE": "HWE", "iscore": "iscore", "trait": "trait"},
):
    query = f"""({columns['MAF']} >= {min_MAF} ) and \
                ({columns['HWE']} >= {min_HWE}) and \
                ({columns['iscore']} >= {min_iscore})"""

    if trait_list is not None:
        query += f" and ({columns['trait']} in {trait_list})"

    return effect.dropna(axis=0, how="any").query(query)


def _calculate_causal_effect(exposure_effect, exposure_error, gwas_effect, gwas_error):
    """Calulcate causal effect from BV effect size and GWAS beta"""

    causal_effect = float(gwas_effect) / float(exposure_effect)

    errors_squared = gwas_error ** 2 / exposure_effect ** 2
    errors_squared += (gwas_effect ** 2) * (exposure_error ** 2) / exposure_effect ** 4
    standard_error = errors_squared ** 0.5

    return causal_effect, standard_error


def _fit_effects(row, columns_bvs, columns_gwas):
    # FIXME: The swap sign here is odd. Why is ref the one with sign -> -1
    return_index = [
        "MR total causal effect",
        "MR se",
        "z score",
        "p value",
        "effect_allele",
    ]

    if row[columns_gwas["allele"]] == row[columns_bvs["alt"]]:
        allele_swap_sign = 1
        effect_allele = "alt"
    elif row[columns_gwas["allele"]] == row[columns_bvs["ref"]]:
        allele_swap_sign = -1
        effect_allele = "ref"
    else:
        return pd.Series(np.nan(len(return_index)), index=return_index)

    causal_es, causal_se = _calculate_causal_effect(
        allele_swap_sign * row[columns_bvs["es"]],
        row[columns_bvs["es_sterr"]],
        row[columns_gwas["beta"]],
        row[columns_gwas["NSE"]],
    )
    z = causal_es / causal_se
    p = norm.sf(abs(z)) * 2.0

    return pd.Series(
        data=[causal_es, causal_se, z, p, effect_allele], index=return_index
    )


def naive_effect_on_trait(
    exposure,
    effect,
    trait_list=None,
    permute=False,
    columns_bvs={},
    columns_gwas={},
    min_MAF=1.0e-3,
    min_HWE=1.0e-50,
    min_iscore=0.9,
):
    """
    MR analysis of a set of SNPs

    Inputs:
    - exposure: dataframe of variants used as exposure variables for the MR analysis, along with their effect on
    - effect: dataframe containing the estimated effect of SNPs on phenotypic trait, as determined by GWAS analysis
    - trait_list: list of traits to select from the effect-SNPs
    - permute: whether to permute the RSIDs for a permutation test
    - min_MAF: required minimum Minor Allele Frequency for inclusion in the analysis
    - min_HWE: required minimum HWE for inclusion in the analysis
    - min_iscore: required minimum iscore for inclusion in the analysis
    """

    # FIXME: the MR here isn't just going by what allele, it's adjusting for the strength of the effect.
    #        it shouldn't matter, since it's just a rescaling, but it's not straight MR.
    # FIXME: ditch the column dicts
    # FIXME: what's HWE and iscore?
    # FIXME: maybe the filtering should be separate from the MR analysis.

    effect = filter_gwas(effect, min_MAF, min_HWE, min_iscore, trait_list, columns_gwas)
    if permute:
        effect[columns_gwas["rsid"]] = np.random.permutation(
            effect[columns_gwas["rsid"]]
        )

    candidates = exposure.merge(effect, left_on="snp", right_on="rsid", how="left")

    out = candidates.join(
        candidates.apply(
            _fit_effects, axis=1, columns_bvs=columns_bvs, columns_gwas=columns_gwas
        )
    )

    # Multiple testing correction with the Benjamini-Hotchberg method
    out["q values"] = multipletests(out["p value"], method="fdr_bh")[1]

    return out
