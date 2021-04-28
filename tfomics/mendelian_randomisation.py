"""MR analysis function for combinging AlleleSeq or BaalChIP data with GeneAtlas GWAS results"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests


def filter_gwas(
    gwas,
    min_MAF,
    min_HWE,
    min_iscore,
    traits,
    columns={"MAF": "MAF", "HWE": "HWE", "iscore": "iscore", "trait": "trait"},
):
    query = f"""({columns['MAF']} >= {min_MAF} ) and \
                ({columns['HWE']} >= {min_HWE}) and \
                ({columns['iscore']} >= {min_iscore}) and\
                ({columns['trait']} in {traits})"""

    return gwas.dropna(axis=0, how="any").query(query)


def mr_fit(row, columns_bvs, columns_gwas):
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

    causal_es, causal_se = causal_effect(
        allele_swap_sign * row[columns_bvs["es"]],
        row[columns_bvs["es_sterr"]],
        row[columns_gwas["beta"]],
        row[columns_gwas["NSE"]],
    )
    causal_es = causal_es * allele_swap_sign
    z = causal_es / causal_se
    p = norm.sf(abs(z)) * 2.0

    return pd.Series(
        data=[causal_es, causal_se, z, p, effect_allele], index=return_index
    )


def naive_mr_analysis(
    df_BVs,
    df_GWAS,
    GWAScodeslist,
    permute=False,
    columns_bvs={},
    columns_gwas={},
    min_MAF=1.0e-3,
    min_HWE=1.0e-50,
    min_iscore=0.9,
):
    """Perform naive MR analysis"""

    # define output format
    out = []
    out_cols = (
        list(df_BVs.columns)
        + list(df_GWAS.columns)
        + ["MR total causal effect", "MR se", "z score", "p value", "effect_allele"]
    )

    # Loop over BVs in the dataframe

    df_GWAS = filter_gwas(
        df_GWAS, min_MAF, min_HWE, min_iscore, GWAScodeslist, columns_gwas
    )
    if permute:
        df_GWAS[columns_gwas["rsid"]] = np.random.permutation(
            df_GWAS[columns_gwas["rsid"]]
        )

    candidates = df_BVs.merge(df_GWAS, left_on="snp", right_on="rsid", how="left")

    out = candidates.join(
        candidates.apply(
            mr_fit, axis=1, columns_bvs=columns_bvs, columns_gwas=columns_gwas
        )
    )
    out["q values"] = list(multipletests(list(out["p value"]), method="fdr_bh")[1])

    return out


def causal_effect(exposure_effect, exposure_error, gwas_effect, gwas_error):
    """Calulcate causal effect from BV effect size and GWAS beta"""

    causal_effect = float(gwas_effect) / float(exposure_effect)

    errors_squared = gwas_error ** 2 / exposure_effect ** 2
    errors_squared += (gwas_effect ** 2) * (exposure_error ** 2) / exposure_effect ** 4
    standard_error = errors_squared ** 0.5

    return causal_effect, standard_error
