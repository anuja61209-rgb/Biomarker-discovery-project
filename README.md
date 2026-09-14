# Biomarker-discovery-project
Biomarker discovery in RIF - Recurrent Implantation failure 
# 🧬 Identification of X-Chromosomal and XCI Escape Biomarkers in Recurrent Implantation Failure (RIF) Using Consensus Machine Learning

![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)
![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Research%20--%20In%20Development-orange.svg)

---

## 📌 Project Overview

**Recurrent Implantation Failure (RIF)** remains a significant challenge in reproductive medicine and assisted reproductive technology (ART)[cite: 3]. Despite the transfer of high-quality embryos, current clinical diagnostics lack the molecular resolution needed to pinpoint subtle endometrial receptivity anomalies[cite: 2, 3].

This project presents a computational and machine learning-driven discovery pipeline focused on **X-chromosome genetics, X-Chromosome Inactivation (XCI) escape mechanics, and microRNA (miRNA) regulatory networks** in endometrial tissue[cite: 2, 3]. By investigating sex-linked dosage compensation alongside autosomal effectors, this study identifies novel candidate biomarkers driving maternal-fetal immunotolerance and structural tissue remodeling during the **Window of Implantation (WOI)**[cite: 2, 3].

---

## 🔬 Key Findings & Candidate Biomarkers

Our consensus machine learning models and multi-cohort transcriptomic validation identified a hybrid regulatory architecture:

### 🌟 Core Diagnostic Panel (Tier 1 & High Priority)
* **`SH2D1A`** *(Xq25)*: Key X-linked immune regulator controlling uterine Natural Killer (uNK) cell activation dynamics and maternal-fetal immune tolerance[cite: 3, 7].
* **`CHIC1`** *(Xq13.2)*: Located adjacent to the X-inactivation center (Xic); serves as a high-novelty tracker for sex-biased dosage discrepancies and epigenetic stability[cite: 2, 7].
* **`IL2RG`** *(Xq21.1)*: XCI escape candidate encoding the common gamma chain ($\gamma_c$) for interleukin receptors, maintaining local Th1/Th2/Th17 immunotolerance balance[cite: 7].
* **`FGF2`** *(4q27)*: Premier autosomal effector target of **miR-3607-5p**, driving endometrial tissue remodeling, stromal cell proliferation, and structural receptivity[cite: 7].

### 🧬 Regulatory Axis
* **`miR-3607-5p`**: Acts as a master post-transcriptional regulator bridging X-linked dosage compensation with downstream autosomal effectors (e.g., `DCAF6`, `NLRP7`, `PEG3`)[cite: 2, 3, 7].

---

## 💻 Machine Learning & Bioinformatics Workflow

The computational workflow integrates feature selection, supervised classification, and downstream functional annotation across independent public datasets (GEO accession **GSE108966**, **GSE111974**, **GSE58144**, **GSE243550**)[cite: 2, 3].


```

┌──────────────────────────────────────────────────────────┐
│ 1. Data Acquisition & Preprocessing                      │
│    • GEO miRNA Datasets (Cohort 1: n=91 | Cohort 2: n=111)│
│    • Standardization (StandardScaler) & Variance Filter  │
└────────────────────────────┬─────────────────────────────┘
│
▼
┌──────────────────────────────────────────────────────────┐
│ 2. Feature Selection & Consensus Modeling                │
│    • Supervised ANOVA (SelectKBest)                      │
│    • Parallel SVM & Logistic Regression Training         │
│    • SHAP (SHapley Additive exPlanations) Feature Mining │
└────────────────────────────┬─────────────────────────────┘
│
▼
┌──────────────────────────────────────────────────────────┐
│ 3. Target Prediction & Network Construction              │
│    • miRDB & TargetScan miRNA-Target Identification      │
│    • STRING Protein-Protein Interaction (PPI) Mapping    │
│    • GeneMANIA Functional Co-expression Analysis         │
└────────────────────────────┬─────────────────────────────┘
│
▼
┌──────────────────────────────────────────────────────────┐
│ 4. Differential Expression & External GEO Validation    │
│    • DESeq2 Transcriptomic Validation                    │
│    • Statistical Benchmarking & AUC Scoring             │
└──────────────────────────────────────────────────────────┘

```

---

## 📊 Summary of Candidate Biomarkers

| Gene Symbol | Chromosome & Context | Diagnostic / Biological Mechanism in RIF | Priority Tier |
| :--- | :--- | :--- | :--- |
| **`SH2D1A`** | Chromosome X (Xq25)[cite: 7] | Regulates uNK cell activation dynamics and immune tolerance[cite: 7]. | **Tier 1 (Core)**[cite: 7] |
| **`CHIC1`** | Chromosome X (Xq13.2)[cite: 7] | Xic-adjacent gene tracking sex-biased epigenetic dosage[cite: 7]. | **Tier 1 (Core)**[cite: 7] |
| **`IL2RG`** | Chromosome X (Xq21.1)[cite: 7] | Interleukin receptor common gamma chain maintaining immune balance[cite: 7]. | **Tier 1 (Core)**[cite: 7] |
| **`FGF2`** | Chromosome 4 (4q27)[cite: 7] | Autosomal angiogenic factor driving structural endometrial remodeling[cite: 7]. | **Tier 1 (Core)**[cite: 7] |
| **`NLRP7`** | Chromosome 19 (19q13.42)[cite: 7] | Activates maternal-effect inflammasome pathways during early implantation[cite: 7]. | **Tier 2 (Secondary)**[cite: 7] |
| **`DCAF6`** | Chromosome 1 (1q24.2)[cite: 7] | Stabilizes X-linked Androgen Receptor (AR) signaling[cite: 7]. | **Tier 2 (Secondary)**[cite: 7] |

---

## 🛠️ Tech Stack & Dependencies

* **Language**: Python 3.9+, R 4.2+
* **Machine Learning**: `scikit-learn`, `shap`, `numpy`, `pandas`
* **Differential Expression**: `DESeq2`, `limma`
* **Target Prediction & Networks**: TargetScan, miRDB, STRING-db, GeneMANIA[cite: 2, 3, 7]
* **Visualization**: `matplotlib`, `seaborn`, `ggplot2`

---



## 📝 Author & Acknowledgments

* **Author**: Anuja Dravid


* **Focus**: Biomarker Discovery, Reproductive Immunology, Machine Learning in Bioinformatics



```

```
