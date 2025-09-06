"""
Component for ProjectHER2 - Statistical Tests (plain functions)

"""

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu

def add_her2_group_by_median(df, signal_col):
    """
    Add a 'her2_group' column using a median cut on the given signal
    """
    
    out = df.copy()
    median_val = out[signal_col].median()
    out["her2_group"] = np.where(out[signal_col] >= median_val, "High", "Low")
    
    return out, float(median_val)

def survival_chi2_fisher(df_with_group):
    """
    Run Chi-square and Fisher's exact tests for her2_group vs vital_status.
    Assumes vital_status is coded 0=alive, 1=deceased
    """
    
    # Build contingency table (rows: group, cols: 0/1)
    ct = pd.crosstab(
        df_with_group["her2_group"],
        df_with_group["vital_status"],
        rownames=["HER2 Group"],
        colnames=["Deceased"]
    )

    # Ensure both columns 0 and 1 exist (ordering alive, deceased)
    for col in [0, 1]:
        if col not in ct.columns:
            ct[col] = 0
    ct = ct[[0, 1]]

    # Chi-square
    chi2, p, dof, expected = chi2_contingency(ct)

    # Fisher's exact on [[High-deceased, High-alive], [Low-deceased, Low-alive]]
    table = np.array([
        [ct.loc["High", 1], ct.loc["High", 0]],
        [ct.loc["Low",  1], ct.loc["Low",  0]],
    ])
    oddsr, fisher_p = fisher_exact(table)

    return {"chi2_p": float(p), "fisher_p": float(fisher_p), "odds_ratio": float(oddsr), "table": ct}

def frac_below(df, drugs, thresh=50):
    """
    Compute fraction of measurements with viability < thresh for each drug
    """
    
    d = df.copy()
    d["drug_name"] = d["drug_name"].astype(str).str.lower()
    out = {}
    
    for drug in drugs:
        name = str(drug).lower()
        sub = d[d["drug_name"] == name]["viability"].dropna()
        out[drug] = float((sub < thresh).mean()) if len(sub) else float("nan")
        
    return out

def mannwhitney_targeted_vs_comparators(df, targeted, comparators, alternative="two-sided"):
    """
    Mann–Whitney U test comparing viability for targeted vs comparator drugs
    """
    
    d = df.copy()
    d["drug_name"] = d["drug_name"].astype(str).str.lower()
    tset = set([s.lower() for s in targeted])
    cset = set([s.lower() for s in comparators])

    targ_vals = d.loc[d["drug_name"].isin(tset), "viability"].dropna().values
    comp_vals = d.loc[d["drug_name"].isin(cset), "viability"].dropna().values

    u_stat, p_val = mannwhitneyu(targ_vals, comp_vals, alternative=alternative)

    return {
        "u_stat": float(u_stat),
        "p_value": float(p_val),
        "n_targeted": int(len(targ_vals)),
        "n_comp": int(len(comp_vals))
    }

def mwu_status_vs_signal(df, signal_col, alternative="greater"):
    """
    Mann–Whitney U test for HER2+ vs HER2− pathway activity
    """
    
    pos_values = df.loc[df["her2_final_status"] == "Positive", signal_col].dropna().values
    neg_values = df.loc[df["her2_final_status"] == "Negative", signal_col].dropna().values

    stat, p_value = mannwhitneyu(pos_values, neg_values, alternative=alternative)

    return {
        "u_stat": float(stat),
        "p_value": float(p_value),
        "median_pos": float(np.median(pos_values)) if len(pos_values) else float("nan"),
        "median_neg": float(np.median(neg_values)) if len(neg_values) else float("nan")
    }