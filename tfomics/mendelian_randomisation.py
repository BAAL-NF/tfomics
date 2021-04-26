"""MR analysis function for combinging AlleleSeq or BaalChIP data with GeneAtlas GWAS results"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.sandbox.stats.multicomp import multipletests


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
        + ["MR total causal effect", "MR se", "z score", "p value"]
    )

    # Loop over BVs in the dataframe  - - - - - - - - - - - - - - - - -

    for _, BV in df_BVs.iterrows():

        # This can be replaced by a join or a merge.
        # select BV matching rsids from GWAS dataframe
        rsid_GWAS = df_GWAS.loc[
            df_GWAS[columns_gwas["rsid"]] == BV[columns_bvs["rsid"]]
        ]

        if not rsid_GWAS.empty:

            # Filtering
            rsid_GWAS = rsid_GWAS.dropna(axis=0, how="any")
            rsid_GWAS = rsid_GWAS.loc[
                (rsid_GWAS[columns_gwas["MAF"]] >= min_MAF)
                & (rsid_GWAS[columns_gwas["HWE"]] >= min_HWE)
                & (rsid_GWAS[columns_gwas["iscore"]] >= min_iscore)
            ]

            # Permute all the rsids in the GWAS results
            if permute:
                rsid_GWAS[columns_gwas["rsid"]] = np.random.permutation(
                    rsid_GWAS[columns_gwas["rsid"]]
                )

            for _, trait_GWAS in rsid_GWAS.iterrows():

                # Only apply calculation if trait is in the traits of interest for study
                if trait_GWAS[columns_gwas["trait"]] in GWAScodeslist:

                    # determine effect allele
                    if trait_GWAS[columns_gwas["allele"]] == BV[columns_bvs["alt"]]:
                        allele_swap_sign = 1
                    elif trait_GWAS[columns_gwas["allele"]] == BV[columns_bvs["ref"]]:
                        allele_swap_sign = -1
                    else:
                        allele_swap_sign = 0

                    # calculate causal effect
                    if allele_swap_sign != 0:
                        causal_es, causal_se = causal_effect(
                            BV[columns_bvs["es"]],
                            BV[columns_bvs["es_sterr"]],
                            trait_GWAS[columns_gwas["beta"]],
                            trait_GWAS[columns_gwas["NSE"]],
                        )
                        z = causal_es / causal_se
                        p = norm.sf(abs(z)) * 2.0
                        out.append(
                            list(BV) + list(trait_GWAS) + [causal_es, causal_se, z, p]
                        )

                    # output nan if the GWAS allele matches neither allele at the BV
                    else:
                        out.append(
                            list(BV) + list(trait_GWAS) + list(np.full(4, np.nan))
                        )

        # output nan if rsid in BV does not exist in GWAS results
        else:
            out.append(
                list(BV)
                + [np.nan for i in range(len(list(df_GWAS.columns)))]
                + list(np.full(4, np.nan))
            )

    # Combine all results into dataframe  - - - - - - - - - - - - - - - - -

    out = pd.DataFrame(data=out, columns=out_cols)

    # Calculate p values and q values  - - - - - - - - - - - - - - - - -

    pvalues = list(out["p value"])
    qvalues = get_q_values(pvalues)
    out["q values"] = qvalues

    return out


def causal_effect(exposure_effect, exposure_error, gwas_effect, gwas_error):
    """Calulcate causal effect from BV effect size and GWAS beta"""

    causal_effect = float(gwas_effect) / float(exposure_effect)

    errors_squared = gwas_error ** 2 / exposure_effect ** 2
    errors_squared += (gwas_effect ** 2) * (exposure_error ** 2) / exposure_effect ** 4
    standard_error = errors_squared ** 0.5

    return causal_effect, standard_error


def get_q_values(list_of_p):
    """Get q values from a list of p values as input
    requires modules loaded:
    from statsmodels.sandbox.stats.multicomp import multipletests"""

    # first record all positions of the float p values
    # get the p-values with None values removed:
    pos_p_values = []
    p_values_no_none = []
    for index, p in enumerate(list_of_p):
        if (type(p) == float or isinstance(p, np.float64)) and not np.isnan(p):
            pos_p_values.append(index)
            p_values_no_none.append(p)

    # convert to q values
    q_values = list(multipletests(p_values_no_none, method="fdr_bh")[1])

    # recover the None values
    out = [np.nan for i in range(len(list_of_p))]

    for pos, q in zip(pos_p_values, q_values):
        out[pos] = q

    return out
