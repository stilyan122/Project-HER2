# Project HER2

Hypothesis-driven statistical analysis of **ERBB2/HER2** using curated datasets—transparent tests (no predictors), clear effect sizes, and reproducible workflows.

---

## 📖 Table of Contents

- [Introduction](#introduction)  
- [Features](#features)  
- [Repository Structure](#repository-structure)  
- [Installation](#installation)  
- [Key Findings](#key-findings)  
- [Contributing](#contributing)  
- [License](#license)  

---

## Introduction

**Project HER2** focuses on **testing biological hypotheses** about ERBB2/HER2 status and pathway activity. We emphasize a clean chain: **assumptions → test choice → effect size → uncertainty**, not prediction. The repo provides a reproducible pipeline for **EDA**, **group comparisons**, **categorical associations**, **effect sizes**, and **(optional) bootstrap** confidence intervals.

---

## Features

- ✅ **Group Comparison (HER2+ vs HER2−)**: **Mann–Whitney U** on pathway/protein signals (e.g., `pp_HER2`, `pp_HER2.pY1248`).  
- ✅ **Categorical Association**: **χ²** test of independence and **Fisher’s exact** (2×2) for median-split pathway level vs `vital.status`.  
- ✅ **Effect Sizes**: **Cramér’s V** (ϕ for 2×2).  
- ✅ **Visualization Helpers**: Quick EDA (counts, histograms, box/violin) with Matplotlib.  
- ✅ **(Optional) Bootstrap**: CIs for selected estimates.  

---

## Repository Structure

```text
Project HER2/
├─ data/
│  ├─ drug-sensitivity.csv
│  └─ mutations.csv
├─ images/
│  ├─ BreastCancer.jpg
│  ├─ ER.jpg
│  ├─ HER2.jpg
│  ├─ HumanCell.jpg
│  ├─ Pathway.jpg
│  ├─ Phosphorylation.jpg
│  └─ PR.jpg
├─ notebooks_code/
│  ├─ eda.py
│  ├─ ProjectHER2.ipynb
│  ├─ stats.py
│  └─ utils.py
└─ .gitignore
```

---

## Installation

1. **Clone the repo**  
   ```bash
   git clone https://github.com/<your-username>/project-her2.git
   cd project-her2
   ```
2. **Create & activate a virtual environment**  
   ```bash
   python3 -m venv venv
   source venv/bin/activate        # Windows: venv\Scripts\activate
   ```
3. **Install dependencies**  
   ```bash
   pip install -r requirements.txt
   ```
   *Requirements include:*  
   - Python 3.10+  
   - pandas, NumPy, SciPy, Matplotlib  
   - (Optional) statsmodels, Jupyter

---

## Key Findings

- **HER2 pathway signal difference (Mann–Whitney U, two-sided)**  
  - `pp_HER2`: **n⁺ = 86, n⁻ = 457**, medians (pos/neg) **1.904 / −0.131**; **U = 33,344**, **p = 1.09×10⁻²⁴**.  
  - `pp_HER2.pY1248`: **n⁺ = 86, n⁻ = 457**, medians (pos/neg) **1.258 / −0.065**; **U = 33,441**, **p = 5.12×10⁻²⁵**.  
  *Interpretation:* HER2-positive samples show markedly higher pathway readouts.

- **Vital status association (2×2, median-split of `pp_HER2`)**  
  - **χ² = 0.056** (df = 1), **p = 0.813**; **Cramér’s V = 0.009**.  
  - Fisher’s exact: **OR = 1.054**, **p = 0.826**.  
  - Contingency (rows High/Low, cols 0/1): `[[307, 46], [304, 48]]`.  
  *Interpretation:* No evidence of association between median-split HER2 level and vital status in this dataset.

---

## Contributing

1. Fork the repo  
2. Create a feature branch (`git checkout -b feature/XYZ`)  
3. Commit your changes (`git commit -m "Add XYZ"`)  
4. Push to your branch (`git push origin feature/XYZ`)  
5. Open a Pull Request

Please follow the existing code style, add tests where appropriate, and update this README if you add new functionality.

---

## License

This project is licensed under the [MIT License](LICENSE). Feel free to use, modify, and distribute under the terms of MIT.
