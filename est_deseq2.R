{
  "cells": [
    {
      "cell_type": "code",
      "source": [
        "# STEP 1\n",
        "install.packages(\"BiocManager\")\n",
        "\n",
        "# STEP 2\n",
        "options(repos = BiocManager::repositories())\n",
        "\n",
        "# STEP 3\n",
        "BiocManager::install(c(\n",
        "  \"DelayedArray\",\n",
        "  \"SummarizedExperiment\",\n",
        "  \"BiocParallel\",\n",
        "  \"DESeq2\"\n",
        "), ask = FALSE, update = FALSE)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "collapsed": true,
        "id": "FL3wWAyD6MKq",
        "outputId": "911dcddc-9f48-46ae-f0e1-0b364f1a1fb1"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "Installing package into ‘/usr/local/lib/R/site-library’\n",
            "(as ‘lib’ is unspecified)\n",
            "\n",
            "'getOption(\"repos\")' replaces Bioconductor standard repositories, see\n",
            "'help(\"repositories\", package = \"BiocManager\")' for details.\n",
            "Replacement repositories:\n",
            "    BioCsoft: https://bioconductor.org/packages/3.23/bioc\n",
            "    BioCann: https://bioconductor.org/packages/3.23/data/annotation\n",
            "    BioCexp: https://bioconductor.org/packages/3.23/data/experiment\n",
            "    BioCworkflows: https://bioconductor.org/packages/3.23/workflows\n",
            "    BioCbooks: https://bioconductor.org/packages/3.23/books\n",
            "    CRAN: https://cran.rstudio.com\n",
            "\n",
            "'getOption(\"repos\")' replaces Bioconductor standard repositories, see\n",
            "'help(\"repositories\", package = \"BiocManager\")' for details.\n",
            "Replacement repositories:\n",
            "    BioCsoft: https://bioconductor.org/packages/3.23/bioc\n",
            "    BioCann: https://bioconductor.org/packages/3.23/data/annotation\n",
            "    BioCexp: https://bioconductor.org/packages/3.23/data/experiment\n",
            "    BioCworkflows: https://bioconductor.org/packages/3.23/workflows\n",
            "    BioCbooks: https://bioconductor.org/packages/3.23/books\n",
            "    CRAN: https://cran.rstudio.com\n",
            "\n",
            "Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.0 (2026-04-24)\n",
            "\n",
            "Warning message:\n",
            "“package(s) not installed when version(s) same as or greater than current; use\n",
            "  `force = TRUE` to re-install: 'DelayedArray' 'BiocParallel'”\n",
            "Installing package(s) 'SummarizedExperiment', 'DESeq2'\n",
            "\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "N9GGIc07YwLq",
        "collapsed": true
      },
      "outputs": [],
      "source": [
        "library(DESeq2)\n",
        "\n",
        "# load data\n",
        "countsS <- read.csv(\"/content/gse108966_esp_counts.csv\", row.names=1)\n",
        "labelsS <- read.csv(\"/content/gse108966_esp_labels.csv\")\n",
        "\n",
        "# transpose counts\n",
        "countsS <- t(countsS)\n",
        "\n"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "collapsed": true,
        "id": "rClHhaLoFzwY"
      },
      "outputs": [],
      "source": [
        "library(DESeq2)\n",
        "\n",
        "# load data\n",
        "counts <- read.csv(\"/content/gse108966_esT_counts.csv\", row.names=1)\n",
        "labels <- read.csv(\"/content/gse108966_esT_labels.csv\")\n",
        "\n",
        "# transpose counts\n",
        "counts <- t(counts)\n",
        "\n"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "# create metadata\n",
        "coldata <- data.frame(\n",
        "    condition = as.factor(labels$condition)\n",
        ")\n",
        "\n",
        "# build DESeq dataset\n",
        "dds <- DESeqDataSetFromMatrix(\n",
        "    countData = counts,\n",
        "    colData = coldata,\n",
        "    design = ~ condition\n",
        ")\n",
        "\n",
        "# run DESeq2\n",
        "dds <- DESeq(dds)\n",
        "\n",
        "# get results\n",
        "res <- results(dds)\n",
        "\n",
        "# sort by adjusted p-value\n",
        "resOrdered <- res[order(res$padj), ]\n",
        "\n",
        "# save results\n",
        "write.csv(as.data.frame(resOrdered),\n",
        "          \"est_DESeq2_results.csv\")"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "collapsed": true,
        "id": "qbiUljaL3WGN",
        "outputId": "ade1cbf2-708b-4d01-be79-d8945d5c5000"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "estimating size factors\n",
            "\n",
            "estimating dispersions\n",
            "\n",
            "gene-wise dispersion estimates\n",
            "\n",
            "mean-dispersion relationship\n",
            "\n",
            "final dispersion estimates\n",
            "\n",
            "fitting model and testing\n",
            "\n",
            "-- replacing outliers and refitting for 7 genes\n",
            "-- DESeq argument 'minReplicatesForReplace' = 7 \n",
            "-- original counts are preserved in counts(dds)\n",
            "\n",
            "estimating dispersions\n",
            "\n",
            "fitting model and testing\n",
            "\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "print(resOrdered)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "collapsed": true,
        "id": "7c3B7U-A7ZGC",
        "outputId": "dddfe210-936f-469b-e994-70e317d63220"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "log2 fold change (MLE): condition 1 vs 0 \n",
            "Wald test p-value: condition 1 vs 0 \n",
            "DataFrame with 2588 rows and 6 columns\n",
            "                 baseMean log2FoldChange     lfcSE       stat      pvalue\n",
            "                <numeric>      <numeric> <numeric>  <numeric>   <numeric>\n",
            "hsa.miR.501.3p  609.60588      -1.514447  0.254244   -5.95667 2.57435e-09\n",
            "hsa.miR.934       3.77960       2.051007  0.373419    5.49251 3.96253e-08\n",
            "hsa.miR.3607.5p   2.55964       2.471624  0.471926    5.23731 1.62931e-07\n",
            "hsa.miR.1827      3.74032      -2.551829  0.504135   -5.06180 4.15318e-07\n",
            "hsa.miR.425.3p   79.39674      -0.752515  0.148374   -5.07175 3.94183e-07\n",
            "...                   ...            ...       ...        ...         ...\n",
            "hsa.miR.936     0.0000000             NA        NA         NA          NA\n",
            "hsa.miR.937.5p  0.0368619       0.062960  2.082696  0.0302301    0.975884\n",
            "hsa.miR.938     0.0000000             NA        NA         NA          NA\n",
            "hsa.miR.9500    0.0000000             NA        NA         NA          NA\n",
            "hsa.miR.96.3p   0.1342945      -0.625466  0.834186 -0.7497921    0.453380\n",
            "                       padj\n",
            "                  <numeric>\n",
            "hsa.miR.501.3p  2.80862e-06\n",
            "hsa.miR.934     2.16156e-05\n",
            "hsa.miR.3607.5p 5.92525e-05\n",
            "hsa.miR.1827    7.55187e-05\n",
            "hsa.miR.425.3p  7.55187e-05\n",
            "...                     ...\n",
            "hsa.miR.936              NA\n",
            "hsa.miR.937.5p           NA\n",
            "hsa.miR.938              NA\n",
            "hsa.miR.9500             NA\n",
            "hsa.miR.96.3p            NA\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "significant <- subset(resOrdered,\n",
        "                      padj < 0.05 &\n",
        "                      abs(log2FoldChange) > 1)"
      ],
      "metadata": {
        "id": "KyAtd7K97zdK"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "head(significant)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/",
          "height": 486
        },
        "id": "YR_6Acte71Bd",
        "outputId": "4a4a12fb-7113-44db-d147-67ae36a64d22"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "display_data",
          "data": {
            "text/plain": [
              "log2 fold change (MLE): condition 1 vs 0 \n",
              "Wald test p-value: condition 1 vs 0 \n",
              "DataFrame with 6 rows and 6 columns\n",
              "                 baseMean log2FoldChange     lfcSE      stat      pvalue\n",
              "                <numeric>      <numeric> <numeric> <numeric>   <numeric>\n",
              "hsa.miR.501.3p  609.60588       -1.51445  0.254244  -5.95667 2.57435e-09\n",
              "hsa.miR.934       3.77960        2.05101  0.373419   5.49251 3.96253e-08\n",
              "hsa.miR.3607.5p   2.55964        2.47162  0.471926   5.23731 1.62931e-07\n",
              "hsa.miR.1827      3.74032       -2.55183  0.504135  -5.06180 4.15318e-07\n",
              "hsa.miR.96.5p   371.12504       -1.19449  0.234643  -5.09068 3.56776e-07\n",
              "hsa.miR.335.3p  298.55079       -1.24556  0.252125  -4.94023 7.80303e-07\n",
              "                       padj\n",
              "                  <numeric>\n",
              "hsa.miR.501.3p  2.80862e-06\n",
              "hsa.miR.934     2.16156e-05\n",
              "hsa.miR.3607.5p 5.92525e-05\n",
              "hsa.miR.1827    7.55187e-05\n",
              "hsa.miR.96.5p   7.55187e-05\n",
              "hsa.miR.335.3p  9.45900e-05"
            ]
          },
          "metadata": {}
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "write.csv(as.data.frame(significant),\n",
        "          \"significant_deseq2_results.csv\")"
      ],
      "metadata": {
        "id": "-UwY__gg8HiG"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\n",
        "    install.packages(\"BiocManager\")\n",
        "}\n",
        "BiocManager::install(\"edgeR\", update = FALSE, ask = FALSE)\n",
        "\n",
        "library(edgeR)\n",
        "\n",
        "# Load data for edgeR\n",
        "edger_counts <- read.csv(\"/content/gse108966_esT_counts.csv\", row.names=1)\n",
        "edger_labels <- read.csv(\"/content/gse108966_esT_labels.csv\")\n",
        "\n",
        "# Transpose counts if necessary (edgeR expects genes as rows, samples as columns)\n",
        "# Assuming the loaded counts are already in the correct orientation (genes as rows, samples as columns).\n",
        "# If not, you might need to transpose: edger_counts <- t(edger_counts)\n",
        "# For now, we will proceed assuming the files are structured correctly for edgeR.\n"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "collapsed": true,
        "id": "ap8oeZjgIwxH",
        "outputId": "23812767-e20f-4ea1-c1db-e5f5b598bd16"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "'getOption(\"repos\")' replaces Bioconductor standard repositories, see\n",
            "'help(\"repositories\", package = \"BiocManager\")' for details.\n",
            "Replacement repositories:\n",
            "    BioCsoft: https://bioconductor.org/packages/3.23/bioc\n",
            "    BioCann: https://bioconductor.org/packages/3.23/data/annotation\n",
            "    BioCexp: https://bioconductor.org/packages/3.23/data/experiment\n",
            "    BioCworkflows: https://bioconductor.org/packages/3.23/workflows\n",
            "    BioCbooks: https://bioconductor.org/packages/3.23/books\n",
            "    CRAN: https://cran.rstudio.com\n",
            "\n",
            "Bioconductor version 3.23 (BiocManager 1.30.27), R 4.6.0 (2026-04-24)\n",
            "\n",
            "Installing package(s) 'edgeR'\n",
            "\n",
            "also installing the dependencies ‘statmod’, ‘limma’\n",
            "\n",
            "\n",
            "Loading required package: limma\n",
            "\n",
            "\n",
            "Attaching package: ‘limma’\n",
            "\n",
            "\n",
            "The following object is masked from ‘package:DESeq2’:\n",
            "\n",
            "    plotMA\n",
            "\n",
            "\n",
            "The following object is masked from ‘package:BiocGenerics’:\n",
            "\n",
            "    plotMA\n",
            "\n",
            "\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "library(edgeR)\n",
        "\n",
        "# Convert to matrix\n",
        "counts <- as.matrix(edger_counts)\n",
        "\n",
        "# Replace NA\n",
        "counts[is.na(counts)] <- 0\n",
        "\n",
        "# IMPORTANT:\n",
        "# transpose because edgeR expects:\n",
        "# rows = genes/miRNAs\n",
        "# cols = samples\n",
        "\n",
        "counts <- t(counts)\n",
        "\n",
        "# Ensure numeric\n",
        "mode(counts) <- \"numeric\"\n",
        "\n",
        "# Ensure integer counts\n",
        "counts <- round(counts)\n",
        "\n",
        "# Check dimensions\n",
        "dim(counts)\n",
        "\n",
        "# Create DGEList\n",
        "dge <- DGEList(counts = counts)\n",
        "\n",
        "# Labels\n",
        "condition <- factor(edger_labels$condition)\n",
        "\n",
        "# Check matching\n",
        "length(condition)\n",
        "ncol(dge)\n",
        "\n",
        "# Must match\n",
        "stopifnot(length(condition) == ncol(dge))\n",
        "\n",
        "# Design matrix\n",
        "design <- model.matrix(~condition)\n",
        "\n",
        "# Filter low-expression miRNAs\n",
        "keep <- filterByExpr(dge, design)\n",
        "\n",
        "dge <- dge[keep,,keep.lib.sizes=FALSE]\n",
        "\n",
        "# Normalize\n",
        "dge <- calcNormFactors(dge)\n",
        "\n",
        "# Estimate dispersion\n",
        "dge <- estimateDisp(dge, design)\n",
        "\n",
        "# Fit model\n",
        "fit <- glmFit(dge, design)\n",
        "\n",
        "# Differential expression\n",
        "lrt <- glmLRT(fit)\n",
        "\n",
        "# Results\n",
        "results <- topTags(lrt, n=Inf)\n",
        "\n",
        "# View\n",
        "head(results$table)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/",
          "height": 497
        },
        "id": "c8Myw67IKa9D",
        "outputId": "96001ca4-2cc7-47ab-a5a4-d48b08f3d4ba"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "display_data",
          "data": {
            "text/html": [
              "<style>\n",
              ".list-inline {list-style: none; margin:0; padding: 0}\n",
              ".list-inline>li {display: inline-block}\n",
              ".list-inline>li:not(:last-child)::after {content: \"\\00b7\"; padding: 0 .5ex}\n",
              "</style>\n",
              "<ol class=list-inline><li>2588</li><li>111</li></ol>\n"
            ],
            "text/markdown": "1. 2588\n2. 111\n\n\n",
            "text/latex": "\\begin{enumerate*}\n\\item 2588\n\\item 111\n\\end{enumerate*}\n",
            "text/plain": [
              "[1] 2588  111"
            ]
          },
          "metadata": {}
        },
        {
          "output_type": "display_data",
          "data": {
            "text/html": [
              "111"
            ],
            "text/markdown": "111",
            "text/latex": "111",
            "text/plain": [
              "[1] 111"
            ]
          },
          "metadata": {}
        },
        {
          "output_type": "display_data",
          "data": {
            "text/html": [
              "111"
            ],
            "text/markdown": "111",
            "text/latex": "111",
            "text/plain": [
              "[1] 111"
            ]
          },
          "metadata": {}
        },
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "calcNormFactors has been renamed to normLibSizes\n",
            "\n"
          ]
        },
        {
          "output_type": "display_data",
          "data": {
            "text/html": [
              "<table class=\"dataframe\">\n",
              "<caption>A data.frame: 6 × 5</caption>\n",
              "<thead>\n",
              "\t<tr><th></th><th scope=col>logFC</th><th scope=col>logCPM</th><th scope=col>LR</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>\n",
              "\t<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>\n",
              "</thead>\n",
              "<tbody>\n",
              "\t<tr><th scope=row>hsa.miR.934</th><td> 1.9681948</td><td>-0.4481669</td><td>44.58712</td><td>2.432873e-11</td><td>1.377006e-08</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.345.5p</th><td> 0.7852197</td><td> 7.8954247</td><td>29.22196</td><td>6.454338e-08</td><td>1.826578e-05</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.501.3p</th><td>-1.5197357</td><td> 6.5272911</td><td>26.21658</td><td>3.051916e-07</td><td>5.757947e-05</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.96.5p</th><td>-1.2962472</td><td> 5.7759800</td><td>24.55046</td><td>7.238902e-07</td><td>1.024305e-04</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.877.5p</th><td>-1.1754730</td><td> 1.2631808</td><td>23.76602</td><td>1.087861e-06</td><td>1.231459e-04</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.30d.5p</th><td> 0.7047399</td><td>11.4124232</td><td>23.03757</td><td>1.588660e-06</td><td>1.497741e-04</td></tr>\n",
              "</tbody>\n",
              "</table>\n"
            ],
            "text/markdown": "\nA data.frame: 6 × 5\n\n| <!--/--> | logFC &lt;dbl&gt; | logCPM &lt;dbl&gt; | LR &lt;dbl&gt; | PValue &lt;dbl&gt; | FDR &lt;dbl&gt; |\n|---|---|---|---|---|---|\n| hsa.miR.934 |  1.9681948 | -0.4481669 | 44.58712 | 2.432873e-11 | 1.377006e-08 |\n| hsa.miR.345.5p |  0.7852197 |  7.8954247 | 29.22196 | 6.454338e-08 | 1.826578e-05 |\n| hsa.miR.501.3p | -1.5197357 |  6.5272911 | 26.21658 | 3.051916e-07 | 5.757947e-05 |\n| hsa.miR.96.5p | -1.2962472 |  5.7759800 | 24.55046 | 7.238902e-07 | 1.024305e-04 |\n| hsa.miR.877.5p | -1.1754730 |  1.2631808 | 23.76602 | 1.087861e-06 | 1.231459e-04 |\n| hsa.miR.30d.5p |  0.7047399 | 11.4124232 | 23.03757 | 1.588660e-06 | 1.497741e-04 |\n\n",
            "text/latex": "A data.frame: 6 × 5\n\\begin{tabular}{r|lllll}\n  & logFC & logCPM & LR & PValue & FDR\\\\\n  & <dbl> & <dbl> & <dbl> & <dbl> & <dbl>\\\\\n\\hline\n\thsa.miR.934 &  1.9681948 & -0.4481669 & 44.58712 & 2.432873e-11 & 1.377006e-08\\\\\n\thsa.miR.345.5p &  0.7852197 &  7.8954247 & 29.22196 & 6.454338e-08 & 1.826578e-05\\\\\n\thsa.miR.501.3p & -1.5197357 &  6.5272911 & 26.21658 & 3.051916e-07 & 5.757947e-05\\\\\n\thsa.miR.96.5p & -1.2962472 &  5.7759800 & 24.55046 & 7.238902e-07 & 1.024305e-04\\\\\n\thsa.miR.877.5p & -1.1754730 &  1.2631808 & 23.76602 & 1.087861e-06 & 1.231459e-04\\\\\n\thsa.miR.30d.5p &  0.7047399 & 11.4124232 & 23.03757 & 1.588660e-06 & 1.497741e-04\\\\\n\\end{tabular}\n",
            "text/plain": [
              "               logFC      logCPM     LR       PValue       FDR         \n",
              "hsa.miR.934     1.9681948 -0.4481669 44.58712 2.432873e-11 1.377006e-08\n",
              "hsa.miR.345.5p  0.7852197  7.8954247 29.22196 6.454338e-08 1.826578e-05\n",
              "hsa.miR.501.3p -1.5197357  6.5272911 26.21658 3.051916e-07 5.757947e-05\n",
              "hsa.miR.96.5p  -1.2962472  5.7759800 24.55046 7.238902e-07 1.024305e-04\n",
              "hsa.miR.877.5p -1.1754730  1.2631808 23.76602 1.087861e-06 1.231459e-04\n",
              "hsa.miR.30d.5p  0.7047399 11.4124232 23.03757 1.588660e-06 1.497741e-04"
            ]
          },
          "metadata": {}
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "edgeR_results <- results$table\n",
        "\n",
        "sig_mirnas <- edgeR_results[\n",
        "    edgeR_results$FDR < 0.05 &\n",
        "    abs(edgeR_results$logFC) > 1,\n",
        "]\n",
        "\n",
        "head(sig_mirnas)"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/",
          "height": 410
        },
        "collapsed": true,
        "id": "3WZNVquIKkPK",
        "outputId": "f9a03e9c-b6b5-4074-e04a-9cf24862ade7"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "display_data",
          "data": {
            "text/html": [
              "<table class=\"dataframe\">\n",
              "<caption>A data.frame: 6 × 5</caption>\n",
              "<thead>\n",
              "\t<tr><th></th><th scope=col>logFC</th><th scope=col>logCPM</th><th scope=col>LR</th><th scope=col>PValue</th><th scope=col>FDR</th></tr>\n",
              "\t<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>\n",
              "</thead>\n",
              "<tbody>\n",
              "\t<tr><th scope=row>hsa.miR.934</th><td> 1.968195</td><td>-0.4481669</td><td>44.58712</td><td>2.432873e-11</td><td>1.377006e-08</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.501.3p</th><td>-1.519736</td><td> 6.5272911</td><td>26.21658</td><td>3.051916e-07</td><td>5.757947e-05</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.96.5p</th><td>-1.296247</td><td> 5.7759800</td><td>24.55046</td><td>7.238902e-07</td><td>1.024305e-04</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.877.5p</th><td>-1.175473</td><td> 1.2631808</td><td>23.76602</td><td>1.087861e-06</td><td>1.231459e-04</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.335.3p</th><td>-1.185642</td><td> 5.3224232</td><td>22.74242</td><td>1.852330e-06</td><td>1.497741e-04</td></tr>\n",
              "\t<tr><th scope=row>hsa.miR.885.5p</th><td> 1.556297</td><td> 0.4081597</td><td>16.85320</td><td>4.038523e-05</td><td>2.285804e-03</td></tr>\n",
              "</tbody>\n",
              "</table>\n"
            ],
            "text/markdown": "\nA data.frame: 6 × 5\n\n| <!--/--> | logFC &lt;dbl&gt; | logCPM &lt;dbl&gt; | LR &lt;dbl&gt; | PValue &lt;dbl&gt; | FDR &lt;dbl&gt; |\n|---|---|---|---|---|---|\n| hsa.miR.934 |  1.968195 | -0.4481669 | 44.58712 | 2.432873e-11 | 1.377006e-08 |\n| hsa.miR.501.3p | -1.519736 |  6.5272911 | 26.21658 | 3.051916e-07 | 5.757947e-05 |\n| hsa.miR.96.5p | -1.296247 |  5.7759800 | 24.55046 | 7.238902e-07 | 1.024305e-04 |\n| hsa.miR.877.5p | -1.175473 |  1.2631808 | 23.76602 | 1.087861e-06 | 1.231459e-04 |\n| hsa.miR.335.3p | -1.185642 |  5.3224232 | 22.74242 | 1.852330e-06 | 1.497741e-04 |\n| hsa.miR.885.5p |  1.556297 |  0.4081597 | 16.85320 | 4.038523e-05 | 2.285804e-03 |\n\n",
            "text/latex": "A data.frame: 6 × 5\n\\begin{tabular}{r|lllll}\n  & logFC & logCPM & LR & PValue & FDR\\\\\n  & <dbl> & <dbl> & <dbl> & <dbl> & <dbl>\\\\\n\\hline\n\thsa.miR.934 &  1.968195 & -0.4481669 & 44.58712 & 2.432873e-11 & 1.377006e-08\\\\\n\thsa.miR.501.3p & -1.519736 &  6.5272911 & 26.21658 & 3.051916e-07 & 5.757947e-05\\\\\n\thsa.miR.96.5p & -1.296247 &  5.7759800 & 24.55046 & 7.238902e-07 & 1.024305e-04\\\\\n\thsa.miR.877.5p & -1.175473 &  1.2631808 & 23.76602 & 1.087861e-06 & 1.231459e-04\\\\\n\thsa.miR.335.3p & -1.185642 &  5.3224232 & 22.74242 & 1.852330e-06 & 1.497741e-04\\\\\n\thsa.miR.885.5p &  1.556297 &  0.4081597 & 16.85320 & 4.038523e-05 & 2.285804e-03\\\\\n\\end{tabular}\n",
            "text/plain": [
              "               logFC     logCPM     LR       PValue       FDR         \n",
              "hsa.miR.934     1.968195 -0.4481669 44.58712 2.432873e-11 1.377006e-08\n",
              "hsa.miR.501.3p -1.519736  6.5272911 26.21658 3.051916e-07 5.757947e-05\n",
              "hsa.miR.96.5p  -1.296247  5.7759800 24.55046 7.238902e-07 1.024305e-04\n",
              "hsa.miR.877.5p -1.175473  1.2631808 23.76602 1.087861e-06 1.231459e-04\n",
              "hsa.miR.335.3p -1.185642  5.3224232 22.74242 1.852330e-06 1.497741e-04\n",
              "hsa.miR.885.5p  1.556297  0.4081597 16.85320 4.038523e-05 2.285804e-03"
            ]
          },
          "metadata": {}
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "# Get significant miRNA names from DESeq2 results\n",
        "deseq2_sig_mirnas <- rownames(significant)\n",
        "\n",
        "# Get significant miRNA names from edgeR results\n",
        "edger_sig_mirnas <- rownames(sig_mirnas)\n",
        "\n",
        "# Find the common miRNAs\n",
        "common_mirnas <- intersect(deseq2_sig_mirnas, edger_sig_mirnas)\n",
        "\n",
        "# Print the common miRNAs and their count\n",
        "print(paste(\"Number of common significant miRNAs:\", length(common_mirnas)))\n",
        "if (length(common_mirnas) > 0) {\n",
        "  print(\"Common significant miRNAs:\")\n",
        "  print(common_mirnas)\n",
        "} else {\n",
        "  print(\"No common significant miRNAs found between DESeq2 and edgeR at the specified thresholds.\")\n",
        "}"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "hXI-KotLKt7T",
        "outputId": "e07efe27-d5dc-4fea-b601-637346aa4658"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "[1] \"Number of common significant miRNAs: 11\"\n",
            "[1] \"Common significant miRNAs:\"\n",
            " [1] \"hsa.miR.501.3p\"  \"hsa.miR.934\"     \"hsa.miR.96.5p\"   \"hsa.miR.335.3p\" \n",
            " [5] \"hsa.miR.877.5p\"  \"hsa.miR.6131\"    \"hsa.miR.503.5p\"  \"hsa.miR.885.5p\" \n",
            " [9] \"hsa.miR.5096\"    \"hsa.miR.502.5p\"  \"hsa.miR.193a.5p\"\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "prev1 = read(\"/content/EST_all_top10_consensus_RIFgse108966_Biomarkers.csv\")"
      ],
      "metadata": {
        "id": "ZDxeax5HKtph"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "# Re-read prev1 to ensure it's available and correctly parsed\n",
        "# Assuming the miRNA names are in the first column of the CSV\n",
        "prev1_data <- read.csv(\"/content/EST_all_top10_consensus_RIFgse108966_Biomarkers.csv\")\n",
        "prev_mirnas <- as.character(prev1_data[[1]])\n",
        "\n",
        "# Find the common miRNAs between common_mirnas and prev_mirnas\n",
        "common_to_all <- intersect(common_mirnas, prev_mirnas)\n",
        "\n",
        "# Print the results\n",
        "print(paste(\"Number of common miRNAs in DESeq2, edgeR, and prev dataset:\", length(common_to_all)))\n",
        "if (length(common_to_all) > 0) {\n",
        "  print(\"Common miRNAs across all datasets:\")\n",
        "  print(common_to_all)\n",
        "} else {\n",
        "  print(\"No common miRNAs found across DESeq2, edgeR, and the prev dataset.\")\n",
        "}"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "8OS0yiroLN2p",
        "outputId": "284e05eb-b0ec-4352-bbeb-ec633f774bbe"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "[1] \"Number of common miRNAs in DESeq2, edgeR, and prev dataset: 0\"\n",
            "[1] \"No common miRNAs found across DESeq2, edgeR, and the prev dataset.\"\n"
          ]
        }
      ]
    }
  ],
  "metadata": {
    "colab": {
      "provenance": [],
      "authorship_tag": "ABX9TyOiX5vDrRWtC8oVz2TvoNR4"
    },
    "kernelspec": {
      "display_name": "R",
      "name": "ir"
    },
    "language_info": {
      "name": "R"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 0
}
