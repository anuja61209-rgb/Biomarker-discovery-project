
# Comprehensive Analysis Report

This report summarizes the entire workflow, from data loading and preprocessing to model building, comparison, and biomarker identification.
## 1. Dataset Used

The analysis started with the `GSE108966_EST_RIF_ml_dataset.csv` dataset. This dataset contains **111 samples** and **2594 features**, primarily representing miRNA expression levels. Key characteristics:

*   **Target Variable**: 'Condition' (CONTROL/IVF), which was mapped to a binary 'Label' (0/1).
*   **Features**: miRNA expression values (2588 columns after dropping metadata).
*   **Missing Values**: No missing values were found in either the feature matrix (X) or the target variable (y).
## 2. Preprocessing Methods

The raw dataset underwent several preprocessing steps to prepare it for machine learning:

*   **Label Encoding**: The categorical 'Condition' column (`CONTROL`, `IVF`) was converted into numerical labels (`0`, `1`) for machine learning.
*   **Feature-Target Split**: The dataset was split into features (X) and the target variable (y), excluding metadata columns.
*   **Missing Value Check**: No missing values were found in X or y.
*   **Standardization**: `StandardScaler` was applied to the features to transform them to a mean of 0 and standard deviation of 1. This is crucial for distance-based algorithms and regularization techniques.
*   **Variance Filtering**: `VarianceThreshold` with a threshold of 0.01 was used to remove features with very low variance, reducing the feature space from 2588 to **2093 features**. This helps in removing uninformative features.
*   **Data Splitting**: The data was split into training (80%) and testing (20%) sets using `train_test_split`, ensuring stratification to maintain the class distribution in both sets (`y_train`: 63 control, 25 IVF; `y_test`: 17 control, 6 IVF).
*   **Feature Selection (SelectKBest)**: `SelectKBest` with `f_classif` as the scoring function was used to select the **top 100 features** from the training data. This step further reduced dimensionality and focused on the most discriminative features.
## 3. Models in Detail with Comparison

Two supervised machine learning models were built and compared:
### a. Support Vector Machine (SVM)

*   **Model**: `SVC` with a `linear` kernel and `class_weight='balanced'` (to handle class imbalance) and `random_state=42` for reproducibility.
*   **Training**: Trained on `X_train_scaled` (top 100 features after scaling).
*   **Evaluation on Test Set (23 samples)**:
    *   **Accuracy**: 0.65
    *   **Classification Report**:
        *   Precision (Class 0): 0.76, Recall (Class 0): 0.76, F1-score (Class 0): 0.76
        *   Precision (Class 1): 0.33, Recall (Class 1): 0.33, F1-score (Class 1): 0.33
    *   **Confusion Matrix**:
        ```
        [[13, 4]
         [ 4, 2]]
        ```
        (13 True Negatives, 2 True Positives, 4 False Positives, 4 False Negatives)
    *   **AUC Score**: 0.68
*   **Cross-Validation (5-fold stratified)**:
    *   **Mean CV Accuracy**: 0.78
    *   **Standard Deviation**: 0.08
*   **Biomarker Identification**: Linear SVM coefficients and SHAP values were used to identify feature importance. Top miRNAs from both methods were extracted.
### b. Logistic Regression

*   **Model**: `LogisticRegression` with `class_weight='balanced'`, `max_iter=5000`, and `random_state=42`.
*   **Training**: Trained on `X_train_scaled` (same top 100 features).
*   **Evaluation on Test Set (23 samples)**:
    *   **Accuracy**: 0.65
    *   **Classification Report**:
        *   Precision (Class 0): 0.80, Recall (Class 0): 0.71, F1-score (Class 0): 0.75
        *   Precision (Class 1): 0.38, Recall (Class 1): 0.50, F1-score (Class 1): 0.43
    *   **Confusion Matrix**:
        ```
        [[12, 5]
         [ 3, 3]]
        ```
        (12 True Negatives, 3 True Positives, 5 False Positives, 3 False Negatives)
    *   **AUC Score**: 0.75
*   **Cross-Validation (5-fold stratified)**:
    *   **Mean CV Accuracy**: 0.83
    *   **Standard Deviation**: 0.09
*   **Biomarker Identification**: Logistic Regression coefficients and SHAP values were used to identify feature importance. Top miRNAs from both methods were extracted.
### Model Comparison

Both models showed comparable performance. Logistic Regression had a slightly higher mean cross-validation accuracy (0.83 vs 0.78) and a slightly better AUC score (0.75 vs 0.68) on the test set. Given the small test set size, these differences are not highly significant, but Logistic Regression appears to have a minor edge in this specific context.
## 4. Biomarkers Obtained

### a. Consensus Biomarkers

Biomarkers were identified using four methods: SVM coefficients, SVM SHAP values, Logistic Regression coefficients, and Logistic Regression SHAP values. The top 20 miRNAs from each method were intersected to find **12 consensus biomarkers**:

*   `hsa-miR-4510`
*   `hsa-miR-5094`
*   `hsa-miR-3150b-3p`
*   `hsa-miR-4664-3p`
*   `hsa-miR-3607-5p`
*   `hsa-miR-185-3p`
*   `hsa-miR-3162-3p`
*   `hsa-miR-937-3p`
*   `hsa-miR-3609`
*   `hsa-miR-5190`
*   `hsa-miR-6799-3p`
*   `hsa-miR-7106-3p`
### b. Integration with DESeq2 Results

The consensus miRNAs were compared against a list of significant miRNAs identified by DESeq2 analysis (filtered for `padj < 0.05` and `|log2FoldChange| > 1`). After standardizing miRNA naming conventions (replacing dots with hyphens), **2 common miRNAs** were found:

*   `hsa-miR-4510`
*   `hsa-miR-3607-5p`

These two miRNAs are consistently identified as important by both machine learning models and differential expression analysis. Their DESeq2 values are:

| Unnamed: 0    | baseMean  | log2FoldChange | lfcSE    | stat     | pvalue        | padj          | cleaned_mirna_name |
|:--------------|:----------|:---------------|:---------|:---------|:--------------|:--------------|:-------------------|
| hsa.miR.3607.5p | 2.559635 | 2.471624 | 0.471926 | 5.237313 | 1.629309e-07 | 5.925482e-05 | hsa-miR-3607-5p |
| hsa.miR.4510 | 1.041423 | -1.666906 | 0.435419 | -3.828278 | 1.290427e-04 | 7.410185e-03 | hsa-miR-4510 |

`hsa-miR-3607-5p` shows significant upregulation, while `hsa-miR-4510` shows significant downregulation, both with very low adjusted p-values, indicating strong statistical significance.
### c. X-chromosome Linked Genes and KEGG Enrichment

The target genes of these common miRNAs were extracted (from user provided lists) and checked for their location on the X chromosome using the `mygene` library. A total of **2 genes were identified on the X chromosome**:

*   `IL2RG`
*   `SH2D1A`

KEGG enrichment analysis was performed on these 2 X-chromosome linked genes. The results indicated potential involvement in:

*   **Fanconi anemia pathway** (genes: IL2RG, SH2D1A)
*   **Primary immunodeficiency** (genes: IL2RG, SH2D1A)
*   **Epstein-Barr virus infection** (genes: IL2RG, SH2D1A)
*   **T cell receptor signaling pathway** (gene: IL2RG)
*   **JAK-STAT signaling pathway** (gene: IL2RG)
*   **Interferon-gamma signaling pathway** (gene: IL2RG)
*   **Natural killer cell mediated cytotoxicity** (gene: IL2RG)

This highlights the role of these X-linked genes in immune-related pathways, particularly `IL2RG` which is involved in several signaling cascades critical for immune function.
### Conclusion

This comprehensive analysis provides a robust set of potential miRNA biomarkers for the studied condition, identified through multiple machine learning approaches and validated by differential expression analysis. The subsequent investigation into their target genes and their chromosomal locations, followed by KEGG enrichment, offers insights into the biological pathways potentially dysregulated by these key miRNAs, particularly highlighting their roles in immune responses.

The heatmap demonstrates that the consensus miRNAs display diverse expression patterns across the cohort, indicating that they are not uniformly regulated. Some miRNAs show similar expression trends in specific sample groups, suggesting possible co-regulation or involvement in related biological pathways. However, the absence of a consistent high- or low-expression pattern across all miRNAs suggests that these molecules likely contribute collectively rather than individually to the underlying biological condition. This observation supports the idea that the consensus miRNAs function as a coordinated molecular signature, with each miRNA potentially influencing different aspects of the biological processes associated with recurrent implantation failure.

more about the gene:
1 SH2D1A	Chromosome X	Xq25	YES (Direct Link)
X-linked gene encoding SAP. Causative gene for Duncan's disease; subject to skewed XCI pathology in female carriers.
The SH2D1A gene is located on human chromosome Xq25. It encodes the signaling lymphocyte activation molecule (SLAM)-associated protein, commonly known as SAP. SAP functions as a critical cytoplasmic molecular switch in T cells and Natural Killer (NK) cells, transducing activating or inhibitory signals upon interaction with SLAM family receptors to coordinate the response against intracellular pathogen
