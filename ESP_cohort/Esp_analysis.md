# Comprehensive Analysis Report

This report summarizes the entire workflow, from data loading and preprocessing to model building, comparison, and biomarker identification.

## 1. Dataset Used

The analysis started with the `GSE108966_RIF_ml_dataset.csv` dataset. This dataset contains **91 samples** and **2594 features**, primarily representing miRNA expression levels. Key characteristics:

*   **Target Variable**: 'Condition' (CONTROL/IVF), which was mapped to a binary 'Label' (0/1).
*   **Features**: miRNA expression values (2588 columns after dropping metadata).
*   **Missing Values**: No missing values were found in either the feature matrix (X) or the target variable (y).

## 2. Preprocessing Methods

The raw dataset underwent several preprocessing steps to prepare it for machine learning:

*   **Label Encoding**: The categorical 'Condition' column (`CONTROL`, `IVF`) was converted into numerical labels (`0`, `1`) for machine learning.
*   **Feature-Target Split**: The dataset was split into features (X) and the target variable (y), excluding metadata columns.
*   **Standardization**: `StandardScaler` was applied to the features to transform them to a mean of 0 and standard deviation of 1. This is crucial for distance-based algorithms and regularization techniques.
*   **Variance Filtering**: `VarianceThreshold` with a threshold of 0.01 was used to remove features with very low variance, reducing the feature space from 2588 to **2236 features**. This helps in removing uninformative features.
*   **Data Splitting**: The data was split into training (80%) and testing (20%) sets using `train_test_split`, ensuring stratification to maintain the class distribution in both sets (`y_train`: 49 control, 23 IVF; `y_test`: 13 control, 6 IVF).
*   **Feature Selection (SelectKBest)**: `SelectKBest` with `f_classif` as the scoring function was used to select the **top 100 features** from the training data. This step further reduced dimensionality and focused on the most discriminative features.

## 3. Models in Detail with Comparison

Two supervised machine learning models were built and compared:

### a. Support Vector Machine (SVM)

*   **Model**: `SVC` with a `linear` kernel and `class_weight='balanced'` (to handle class imbalance) and `random_state=42` for reproducibility.
*   **Training**: Trained on `X_train_scaled` (top 100 features after scaling).
*   **Evaluation on Test Set (19 samples)**:
    *   **Accuracy**: 0.79
    *   **Classification Report**:
        *   Precision (Class 0): 0.85, Recall (Class 0): 0.85, F1-score (Class 0): 0.85
        *   Precision (Class 1): 0.67, Recall (Class 1): 0.67, F1-score (Class 1): 0.67
    *   **Confusion Matrix**:
        ```
        [[11, 2]
         [ 2, 4]]
        ```
        (11 True Negatives, 4 True Positives, 2 False Positives, 2 False Negatives)
    *   **AUC Score**: 0.73
*   **Cross-Validation (5-fold stratified)**:
    *   **Mean CV Accuracy**: 0.74
    *   **Standard Deviation**: 0.13
*   **Biomarker Identification**: Linear SVM coefficients and SHAP values were used to identify feature importance. Top miRNAs from both methods were extracted.

### b. Logistic Regression

*   **Model**: `LogisticRegression` with `class_weight='balanced'`, `max_iter=5000`, and `random_state=42`.
*   **Training**: Trained on `X_train_scaled` (same top 100 features).
*   **Evaluation on Test Set (19 samples)**:
    *   **Accuracy**: 0.74
    *   **Classification Report**:
        *   Precision (Class 0): 0.83, Recall (Class 0): 0.77, F1-score (Class 0): 0.80
        *   Precision (Class 1): 0.57, Recall (Class 1): 0.67, F1-score (Class 1): 0.62
    *   **Confusion Matrix**:
        ```
        [[10, 3]
         [ 2, 4]]
        ```
        (10 True Negatives, 4 True Positives, 3 False Positives, 2 False Negatives)
    *   **AUC Score**: 0.76
*   **Cross-Validation (5-fold stratified)**:
    *   **Mean CV Accuracy**: 0.82
    *   **Standard Deviation**: 0.12
*   **Biomarker Identification**: Logistic Regression coefficients and SHAP values were used to identify feature importance. Top miRNAs from both methods were extracted.

### Model Comparison

Both models showed comparable performance. Logistic Regression had a slightly higher mean cross-validation accuracy (0.82 vs 0.74) and a slightly better AUC score (0.76 vs 0.73) on the test set. Given the small test set size, these differences are not highly significant, but Logistic Regression appears to have a minor edge in this specific context.

## 4. Biomarkers Obtained

### a. Consensus Biomarkers

Biomarkers were identified using four methods: SVM coefficients, SVM SHAP values, Logistic Regression coefficients, and Logistic Regression SHAP values. The top 20 miRNAs from each method were intersected to find **10 consensus biomarkers**:

*   `hsa-miR-548f-3p`
*   `hsa-miR-3146`
*   `hsa-miR-5095`
*   `hsa-miR-5701`
*   `hsa-miR-1273c`
*   `hsa-miR-3607-3p`
*   `hsa-miR-7152-5p`
*   `hsa-miR-548am-3p`
*   `hsa-miR-4478`
*   `hsa-miR-29b-3p`

### b. Integration with DESeq2 Results

The consensus miRNAs were compared against a list of significant miRNAs identified by DESeq2 analysis (filtered for `padj < 0.05` and `|log2FoldChange| > 1`). After standardizing miRNA naming conventions (replacing dots with hyphens), **2 common miRNAs** were found:

*   `hsa-miR-5701`
*   `hsa-miR-3607-3p`

These two miRNAs are consistently identified as important by both machine learning models and differential expression analysis. Their DESeq2 values are:

| Unnamed: 0    | baseMean  | log2FoldChange | lfcSE    | stat     | pvalue        | padj          | cleaned_mirna_name |
|:--------------|:----------|:---------------|:---------|:---------|:--------------|:--------------|:-------------------|
| hsa.miR.5701  | 83.148010 | 1.796157       | 0.220865 | 8.132365 | 4.209939e-16  | 2.805924e-13  | hsa-miR-5701       |
| hsa.miR.3607.3p | 43.617679 | 1.513542       | 0.225559 | 6.710170 | 1.943972e-11  | 3.701878e-09  | hsa-miR-3607-3p    |

Both miRNAs show significant upregulation (`log2FoldChange > 0`) and very low adjusted p-values, indicating strong statistical significance.

### c. X-chromosome Linked Genes and KEGG Enrichment

The target genes of these common miRNAs were extracted (from miRDB and TargetScan) and checked for their location on the X chromosome using the `mygene` library. A total of **6 genes were identified on the X chromosome**:

*   `CHIC1`
*   `ZNF711`
*   `MOSPD1`
*   `PABPC5`
*   `PAGE5`
*   `P2RY10`

KEGG enrichment analysis was performed on these 6 X-chromosome linked genes. The results indicated potential involvement in:

*   **RNA degradation** (gene: `PABPC5`): Highly significant enrichment, suggesting `PABPC5` plays a central role.
*   **mRNA surveillance pathway** (gene: `PABPC5`): Significant enrichment, reinforcing the role of `PABPC5` in mRNA quality control.
*   **RNA transport** (gene: `PABPC5`): Some enrichment, indicating possible involvement in RNA movement.
*   **Neuroactive ligand-receptor interaction** (gene: `P2RY10`): Suggests involvement in cell signaling pathways.

This highlights `PABPC5` as a key X-linked target gene involved in crucial RNA metabolic processes. The other X-linked genes did not show significant enrichment in these specific pathways.

### Conclusion

This comprehensive analysis provides a robust set of potential miRNA biomarkers for the studied condition, identified through multiple machine learning approaches and validated by differential expression analysis. The subsequent investigation into their target genes and their chromosomal locations, followed by KEGG enrichment, offers insights into the biological pathways potentially dysregulated by these key miRNAs.
futhur  analysis

1.   MOSPD1	Chromosome X	Xq26.3	YES (Direct Link)
X-linked motile sperm domain-containing protein 1; regulates membrane contact sites. Susceptible to epigenetic reactivation in breast cancer.
2. PIH1D3	Chromosome X	Xq22.3	YES (Direct Link)
X-linked gene encoding DNAAF6. Mutated in primary ciliary dyskinesia; shows random XCI mosaicism in female carriers.
3. P2RY10	Chromosome X	Xq21.1	YES (Direct Link)
X-linked G-protein coupled purinergic receptor. Highly validated escapee from X-chromosome inactivation in lymphocytes.
4. PAGE5	Chromosome X	Xq21.1	YES (Direct Link)
X-linked cancer-testis antigen subject to XCI. Targeted by hsa-miR-3607-3p; reactivated upon loss of XCI in tumors.
5. FAM122B	Chromosome X	Xq26.3	YES (Direct Link)
X-linked regulator of cell growth; adjacent to MOSPD1 and shares regulatory enhancer domains.
6. GAGE12B	Chromosome X	Xp11.23	YES (Direct Link)
X-linked cancer-testis antigen of the repetitive GAGE family; reactivated upon chromatin relaxation of the Xi in cancer cells.
7. ZNF711	Chromosome X	Xq21.1	YES (Direct Link)
X-linked zinc finger transcription factor associated with intellectual disability. Upregulated by genomic "loss of XCI" in ovarian cancers.
8. CHIC1	Chromosome X	Xq22.1	YES (Direct Link)
X-linked cysteine-rich hydrophobic domain protein; maps near the Xic center and coordinates the establishment of XCI.
9. RP6-24A23.6	Chromosome X	Xq22.3	YES (Direct Link)
Uncharacterized X-linked protein-coding transcript ENSG00000260548. Shows strong co-expression with autosomal TCF24.

Some functional linkage :
1. TCF24	Chromosome 20	20q11.21	Functional Linkage
Autosomal transcription factor with high expression co-correlation with the X-linked gene RP6-24A23.6.
2. DCAF6	Chromosome 1	1q24.2	Functional Linkage
Autosomal substrate receptor of CRL4 E3 ligase ; directly binds and stabilizes the X-linked Androgen Receptor against degradation.



