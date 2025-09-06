"""
Component for Utilies in Dataset logic
===================================================

"""

import re
import pandas as pd
from pathlib import Path

def to_snake(col):
    """
    Convert a column name to snake_case. Returns the snake_cased column name.
    """
    
    s = re.sub(r"[.\-\s]+", "_", col)            
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s) 
    s = s.lower()
    s = re.sub(r"_+", "_", s).strip("_")
    
    return s

def load_data(data_dir):
    """
    Load mutations.csv and drug-sensitivity.csv from data_dir.
    """
    
    data_dir = Path(data_dir)
    mut_path = data_dir / "mutations.csv"
    drug_path = data_dir / "drug-sensitivity.csv"
    
    if not mut_path.exists():
        raise FileNotFoundError("Missing file: {}".format(mut_path))
        
    if not drug_path.exists():
        raise FileNotFoundError("Missing file: {}".format(drug_path))
        
    return pd.read_csv(mut_path), pd.read_csv(drug_path)

def resolve_signal(columns, preference=None):
    """
    Choose the HER2 pathway signal column from available columns.
    """
    
    if preference is None:
        preference = ["pp_her2", "pp_her2_py1248"]
        
    cols = set(columns)
    
    for cand in preference:
        if cand in cols:
            return cand

def clean_mutations(mutations_df, signal_preference=None):
    """
    Clean tumor/mutations dataframe and select relevant columns.
    """

    # Normalize columns to snake_case
    df = mutations_df.copy()
    df.columns = [to_snake(c) for c in df.columns]

    # Resolve the HER2 status column (support common variants)
    her2_col_candidates = ["her2_final_status", "her2_status", "her2_status_final"]
    her2_col = next((c for c in her2_col_candidates if c in df.columns), None)
    if her2_col is None:
        raise ValueError(
            "Could not find a HER2 status column among: {}. "
            "Saw columns like: {}".format(
                ", ".join(her2_col_candidates), ", ".join(list(df.columns)[:12])
            )
        )

    # Standardize HER2 labels; keep Positive/Negative
    df[her2_col] = df[her2_col].astype(str).str.strip().str.capitalize()
    df = df[df[her2_col].isin(["Positive", "Negative"])].copy()

    # Choose pathway signal
    signal_col = resolve_signal(df.columns, signal_preference)

    # Final column set (only those that exist)
    wanted = [her2_col, signal_col, "er_status", "pr_status", "vital_status", "histological_type"]
    cols_final = [c for c in wanted if c in df.columns]
    out = df[cols_final].dropna(subset=[signal_col]).copy()

    # Downstream expects this exact name
    if her2_col != "her2_final_status":
        out = out.rename(columns={her2_col: "her2_final_status"})

    # Add H indicator (1 if HER2 Positive else 0)
    out["H"] = (
        out["her2_final_status"].astype(str).str.strip().str.lower().isin(
            ["positive", "pos", "her2+", "1", "true", "yes", "high"]
        )
    ).astype(int)

    return out, signal_col

def clean_drugs(drug_df):
    """
    Clean drug sensitivity dataframe.
    """
    
   # Work on a copy and snake_case all headers
    df = drug_df.copy()
    df.columns = [to_snake(c) for c in df.columns]

    # Sanity check for essential columns
    if "drug_name" not in df.columns or "viability" not in df.columns:
        raise ValueError("Drug dataframe must contain 'drug_name' and 'viability'.")

    # Normalize drug names
    df["drug_name"] = (
        df["drug_name"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    # Drop rows missing essential fields
    df = df.dropna(subset=["drug_name", "viability"])

    # Keep only the important columns (keep those that exist)
    keep_cols = [c for c in ["cosmic_id", "drug_name", "dose", "viability"] if c in df.columns]
    df = df[keep_cols]
    
    return df

def encode_vital_status(df, col="vital_status"):
    """
    Encode vital_status to {0,1} where 0=alive and 1=deceased.
    Raises if unmapped values exist. If column missing, returns df unchanged
    """
    
    out = df.copy()
    
    if col not in out.columns:
        return out

    # Create mapping with common words
    mapping = {"dead":1,"deceased":1,"1":1,1:1,True:1,"alive":0,"0":0,0:0,False:0}
    
    def map_one(x):
        """ 
        Little Helper Function for Core mapping
        """
        
        return mapping.get(str(x).strip().lower(), None)
        
    out[col] = out[col].map(map_one)
    
    if out[col].isna().any():
        bad_rows = out.index[out[col].isna()].tolist()[:5]
        raise ValueError("Unmapped vital_status values at rows: {}".format(bad_rows))
        
    out[col] = out[col].astype(int)
    return out

def validate_drugs_numeric(df, clip_viability=(0, 200), min_dose=0, keep_zero_dose=True):
    """
    Coerce 'dose' and 'viability' to numeric; clip viability; drop invalid doses.
    Returns a cleaned copy of the input df
    """
    
    out = df.copy()
    
    if "dose" in out.columns:
        out["dose"] = pd.to_numeric(out["dose"], errors="coerce")
        
    if "viability" in out.columns:
        out["viability"] = pd.to_numeric(out["viability"], errors="coerce").clip(clip_viability[0], clip_viability[1])
   
    out = out.dropna(subset=[c for c in ["drug_name","viability"] if c in out.columns])
    
    if "dose" in out.columns:
        if keep_zero_dose:
            out = out[out["dose"].fillna(-1) >= min_dose]
        else:
            out = out[out["dose"].fillna(-1) > min_dose]
    return out

def missingness_report(df, top_n=25):
    """
    Return a compact table of: column, null_pct, n_unique, example (top_n rows)
    """
    
    rows = []
    
    for c in df.columns:
        s = df[c]
        rows.append({
            "column": c,
            "null_pct": round(float(s.isna().mean()), 3),
            "n_unique": int(s.nunique(dropna=True)),
            "example": (s.dropna().iloc[0] if s.notna().any() else None),
        })
        
    rep = pd.DataFrame(rows).sort_values("null_pct", ascending=False)
    return rep.head(top_n)