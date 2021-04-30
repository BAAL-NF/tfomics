"""MR analysis function for combinging AlleleSeq or BaalChIP data with GeneAtlas GWAS results.
Intended use: 
- Filter your GWAS results using filter_effect_snps
- Calculate naive_effect_on_trait using your list of binding variants and 
"""
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests


def filter_effect_snps(
    effect: pd.DataFrame,
    min_MAF: float = 1.0e-3,
    min_HWE: float = 1.0e-50,
    min_iscore: float = 0.9,
    trait_list: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Apply a set of filters to a GWAS dataframe

    Input:
    - effect: dataframe of GWAS results
    - min_MAF: required minimum Minor Allele Frequency for inclusion in the analysis
    - min_HWE: required minimum HWE for inclusion in the analysis
    - min_iscore: required minimum iscore for inclusion in the analysis
    - trait_list: list of traits to select from the effect-SNPs

    Output:
    - Filtered set of SNPs from the GWAS dataframe
    """
    # FIXME: what's HWE and iscore?

    query = f"""(MAF >= {min_MAF} ) and \
                (HWE >= {min_HWE}) and \
                (iscore >= {min_iscore})"""

    if trait_list is not None:
        query += f" and (trait in {trait_list})"

    return effect.dropna(axis=0, how="any").query(query)


def _calculate_causal_effect(
    exposure_effect: float, exposure_error: float, gwas_effect: float, gwas_error: float
) -> Tuple[float, float]:
    """Calulcate causal effect from BV effect size and GWAS beta"""

    causal_effect = gwas_effect / exposure_effect

    errors_squared = gwas_error ** 2 / exposure_effect ** 2
    errors_squared += (gwas_effect ** 2) * (exposure_error ** 2) / exposure_effect ** 4
    standard_error = errors_squared ** 0.5

    return causal_effect, standard_error


def _fit_effects(row: pd.Series) -> pd.Series:
    # FIXME: The swap sign here is odd. Why is ref the one with sign -> -1
    return_index = [
        "MR total causal effect",
        "MR se",
        "z score",
        "p value",
        "effect_allele",
    ]

    # Determine the effect allele, and orient the effect accordingly.
    if row.allele == row.alt:
        allele_swap_sign = 1
        effect_allele = "alt"
    elif row.allele == row.ref:
        allele_swap_sign = -1
        effect_allele = "ref"
    else:
        return pd.Series(np.full(len(return_index), np.nan), index=return_index)

    causal_es, causal_se = _calculate_causal_effect(
        allele_swap_sign * row.es,
        row.es_sterr,
        row.beta,
        row.NSE,
    )
    # Calculate the z-score and associated p-value
    z = causal_es / causal_se
    p = norm.sf(abs(z)) * 2.0

    return pd.Series(
        data=[causal_es, causal_se, z, p, effect_allele], index=return_index
    )


def naive_effect_on_trait(
    exposure: pd.DataFrame,
    effect: pd.DataFrame,
    permute: bool = False,
) -> pd.DataFrame:
    """
    MR analysis of a set of SNPs. This function takes two sets of SNPs for effect and exposure measurements and performs
    a mendelian randomisation analysis for all snp-trait pairs. It uses the Benjamini-Hochberg method to determine
    q-values for each SNP based on the distribution of p-values in the analysis.

    If permute is set, the analysis is treated as a permutation test, and the effect SNPs are randomly reshuffled.

    Inputs:
    - exposure: dataframe of variants used as exposure variables for the MR analysis, along with their effect on
    - effect: dataframe containing the estimated effect of SNPs on phenotypic trait, as determined by GWAS analysis
    - permute: whether to permute the RSIDs for a permutation test

    Output:
    - Set of SNPs after MR analysis, including effect size, p-values and q-values.
    """

    # FIXME: the MR here isn't just going by what allele, it's adjusting for the strength of the effect.
    #        it shouldn't matter, since it's just a rescaling, but it's not straight MR.
    if permute:
        effect[["rsid"]] = np.random.permutation(effect[["rsid"]])

    candidates = exposure.merge(effect, on="rsid", how="left")

    out = candidates.join(candidates.apply(_fit_effects, axis=1)).dropna()

    # Multiple testing correction with the Benjamini-Hotchberg method
    out["q values"] = multipletests(out["p value"], method="fdr_bh")[1]

    return out
