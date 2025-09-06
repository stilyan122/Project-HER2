"""
Component ProjectHER2 EDA helpers (matplotlib-only, plain functions, no type hints).

"""

import matplotlib.pyplot as plt
import numpy as np

def barplot_counts(df, column, title=None):
    """
    Display a bar plot for the value counts of a given column
    """
    
    # Count including NaN (shown as 'nan')
    counts = df[column].value_counts(dropna=False)

    # Basic bar chart
    plt.figure()
    plt.bar(counts.index.astype(str), counts.values)
    plt.title(title or f"Distribution of {column}")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def boxplot_by_group(df, signal, group, title=None):
    """
    Draw a boxplot of a numeric signal grouped by a categorical column
    """
    
    # Build list of values per group label
    grouped = df.groupby(group)[signal].apply(list)
    labels = grouped.index.tolist()
    data = [grouped[label] for label in labels]

    # compute n per group
    counts = df.groupby(group)[signal].apply(lambda s: s.notna().sum()).reindex(labels)
    lab_with_n = [f"{lab}\n(n={int(counts.loc[lab])})" for lab in labels]

    plt.figure()
    plt.boxplot(data, labels=lab_with_n, patch_artist=True, showfliers=False)
    plt.title(title or f"{signal} by {group}")
    plt.ylabel(signal)
    plt.xlabel(group)
    plt.tight_layout()
    plt.show()

def hist_density(df, column, bins=40):
    """
    Plot a histogram with density normalization for one numeric column
    """
    
    data = df[column].dropna()

    plt.figure()
    plt.hist(data, bins=bins, density=True, alpha=0.6)
    plt.title(f"Histogram & Density of {column}")
    plt.xlabel(column)
    plt.ylabel("Density")
    plt.tight_layout()
    plt.show()

def bar_top_drugs_by_count(drug_df, top_n=20):
    """
    Plot a bar chart of the top-N drugs by number of rows/measurements
    """
    
    freq = (
        drug_df["drug_name"]
        .value_counts()
        .rename_axis("drug")
        .reset_index(name="count")
    )

    plt.figure()
    plt.bar(freq["drug"].iloc[:top_n], freq["count"].iloc[:top_n])
    plt.xticks(rotation=45, ha="right")
    plt.title(f"Top {top_n} Drugs by Measurement Count")
    plt.ylabel("Count")
    plt.xlabel("Drug Name")
    plt.tight_layout()
    plt.show()

def median_response_logdose(df, drugs):
    """
    Plot median viability vs log10(dose) curves for a list of drugs
    """
    
    plt.figure()

    # Normalize drug names for matching
    d = df.copy()
    d["drug_name"] = d["drug_name"].astype(str).str.lower()

    # Track global x-range for nicer ticks
    all_log10 = []

    for drug in drugs:
        name = str(drug).lower()

        med = (
            d[d["drug_name"] == name]
            .groupby("dose")["viability"]
            .median()
            .reset_index()
        )
        if med.empty:
            continue

        med["log10_dose"] = np.log10(med["dose"])
        all_log10.extend(med["log10_dose"].tolist())

        plt.plot(
            med["log10_dose"],
            med["viability"],
            marker="o",
            label=drug.title()
        )

    plt.xlabel("log\u2081\u2080(Dose)")
    plt.ylabel("Median Viability")
    plt.title("Median Dose–Response Curves (log₁₀ dose)")

    # Nice integer ticks across the plotted range
    if len(all_log10) > 0:
        lo, hi = np.floor(min(all_log10)), np.ceil(max(all_log10))
        plt.xticks(np.arange(lo, hi + 1))

    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_ecdf(df, drug_names):
    """
    Draw ECDF curves of viability for a set of drugs
    """
    
    # Normalize names for filtering
    d = df.copy()
    d["drug_name"] = d["drug_name"].astype(str).str.lower()
    targets = [str(x).lower() for x in drug_names]

    plt.figure()
    for drug in targets:
        vals = np.sort(d.loc[d["drug_name"] == drug, "viability"].dropna().values)
        if len(vals) == 0:
            continue
        y = np.arange(1, len(vals) + 1) / len(vals)
        plt.step(vals, y, where="post", label=drug.title())

    plt.xlabel("Viability")
    plt.ylabel("ECDF")
    plt.title("ECDF of Viability for Selected Drugs")
    plt.legend()
    plt.tight_layout()
    plt.show()

def violin_viability(df, drugs):
    """
    Plot violin distributions of viability for specific drugs
    """
    
    # Normalize names
    d = df.copy()
    d["drug_name"] = d["drug_name"].astype(str).str.lower()
    targets = [str(x).lower() for x in drugs]

    # Collect viability arrays in the same order
    data = [d.loc[d["drug_name"] == drug, "viability"].dropna().values for drug in targets]
    labels = [drug.title() for drug in targets]

    plt.figure()
    plt.violinplot(data, showmedians=True)
    plt.xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right")
    plt.ylabel("Viability")
    plt.title("Viability Distributions (Violin Plot)")
    plt.tight_layout()
    plt.show()