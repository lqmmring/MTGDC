# MTGDC
===========================================================================================
> MTGDC is an unsupervised clustering framework for single-cell RNA-seq data. MTGDC can learn global topological information between cells from multiple scales and use a simple and efficient tensor diffusion update algorithm to spread the high-order cell relationship graph to its neighbors until convergence to global stable state which preserves local and global cell topology structure.To achieve the purpose of mining potential similarity distributions among cells under a large amount of noise, we design a multi-scale affinity learning method to construct a fully connected graph between cells.

## Overview
![Overview](figures/overview.png)

## Table of content
- [Introduction](#Introduction)
- [Installation](#Installation)
    - [Requirement](#--required-installations)
    - [Data](#--Data)
    - [Files](#--Files)
- [Example Usage](#example-usage)
    - [Preprocessing](#--Preprocessing)
    - [AIDE](#--AIDE)
    - [RPH-kmeans](#--RPH-kmeans)
    - [Biological Analysis](#--biological-analysis)
    - [Example Results](#--examples-results)
    - [Scalability](#--scalability)
- [Citation](#citation-&-references)
- [Acknowledgment](#Acknowledgment)
- [Maintenance](#Maintenance)

## Introduction
The algorithm has the following mechanisms: Multi-scale Affinity Learning, Tensor Graph Diffusion Learning, and Mixture Operator.

- **Multi-scale Affinity Learning**: To achieve the purpose of mining potential similarity distributions among cells under a large amount of noise, we design a multi-scale affinity learning method to construct a fully connected graph between cells.
- **Tensor Graph Diffusion Learning**: For each affinity matrix, we propose an efficient tensor graph diffusion learning framework to learn high-order in context of cells with multi-scale affinity matrices.
- **Mixture Operator**: Finally, we mix the multi-scale tensor graph together to obtain the final complete high-order affinity matrix and apply it to spectral clustering. 

## Installation
### - Required Installations
The software is coded using MATLAB and is free for academic purposes.

- MATLAB: MATLAB R2019b
- R: Seurat, splatter, ggplot2

Part of the code are from following paper:
###### Bai S, Zhou Z, Wang J, et al. Ensemble diffusion for retrieval[C]//Proceedings of the IEEE International Conference on Computer Vision. 2017: 774-783.

### - Install

To use, please download the MTGDC folder and follow the instructions

### - Data

### - Real Data
we selected 12 public scRNA-seq datasets to verify the performance of clustering analysis

- [Pollen](https://www.ncbi.nlm.nih.gov/sra?term=SRP041736)
###### 1. Pollen, A.A., et al., Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. Nature biotechnology, 2014. 32(10): p. 1053.
- [Deng](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719)
###### 2. Deng, Q., et al., Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells. Science, 2014. 343(6167): p. 193-196.
- [Ginhoux](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60783)
###### 3. Schlitzer, A., et al., Identification of cDC1-and cDC2-committed DC progenitors reveals early lineage priming at the common DC progenitor stage in the bone marrow. Nature immunology, 2015. 16(7): p. 718.
- [Buettner](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2512)
###### 4. Buettner, F., et al., Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature biotechnology, 2015. 33(2): p. 155.
- [Ting](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51372)
###### 5. Ting, D.T., et al., Single-cell RNA sequencing identifies extracellular matrix gene expression by pancreatic circulating tumor cells. Cell reports, 2014. 8(6): p. 1905-1918.
- [Treutlin](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52583)
###### 6. Treutlein, B., et al., Reconstructing lineage hierarchies of the distal lung epithelium using single-cell RNA-seq. Nature, 2014. 509(7500): p. 371.
- [Kolod](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/)
###### 7. Kolodziejczyk, Aleksandra A., et al., Single Cell RNA-Sequencing of Pluripotent States Unlocks Modular Transcriptional Variation. Cell Stem Cell, 2015. 17(4): p. 471-485.
- [mECS](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74535)
###### 8. Angermueller, C., et al., Parallel single-cell sequencing links transcriptional and epigenetic heterogeneity. Nature methods, 2016. 13(3): p. 229-232.
- [Usoskin](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59739)
###### 9. Usoskin, D., et al., Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nature neuroscience, 2015. 18(1): p. 145-153.
- [Tasic](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585)
###### 10. Tasic, B., et al., Adult mouse cortical cell taxonomy revealed by single cell transcriptomics. Nature neuroscience, 2016. 19(2): p. 335-346.
- [Macosko](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63473)
###### 11. Macosko, E.Z., et al., Highly parallel genome-wide expression profiling of individual cells using nanoliter droplets. Cell, 2015. 161(5): p. 1202-1214.
- [Zeisel](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361)
###### 12. Zeisel, A., et al., Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science, 2015. 347(6226): p. 1138-1142.

The data is stored as the mat format in the folder [Data](https://github.com/lqmmring/MTGDC/tree/main/Data).

### - Simulated Data
we evaluated our method on simulated datasets. Synthetic datasets were simulated by the R package [*Splatter*](https://github.com/Oshlack/splatter-paper). Thed parameters of simulated data are provided in [splatter](https://github.com/lqmmring/MTGDC/blob/main/splatter.R). 

The data is stored as the mat format in the folder [Sim-Data](https://github.com/lqmmring/MTGDC/tree/main/Data/Sim-Data) 

### - Files

- [MTGDC](https://github.com/lqmmring/MTGDC/blob/main/MTGDC.m): main MTGDC algorithm consisting of the three steps.
- [NormalizeFea](https://github.com/lqmmring/MTGDC/blob/main/LIB/NormalizeFea.m): provide the normalized processing. 
- [adaptiveGaussian](https://github.com/lqmmring/MTGDC/blob/main/LIB/adaptiveGaussian.m): provide the multi-scale affinity learning.
- [IterativeDiffusionTPGKNN](https://github.com/lqmmring/MTGDC/blob/main/LIB/IterativeDiffusionTPGKNN.m):this code is an implementation of the diffusion process.
- [knnSparse](https://github.com/lqmmring/MTGDC/blob/main/LIB/knnSparse.m): sparse the affinity matrix by k-nearst neighbor.
- [mergeW](https://github.com/lqmmring/MTGDC/blob/main/LIB/mergeW.m): this code is an implementation of the mixture of diffusion affinity matrices.

### - Baselines
To verify the performance of our method (MTGDC), we compared it with some competitive baselines.

We selected several widely used scRNA-seq data clustering tools, including graph-based methods (Seurat and SNN-Cliq), ensemble-based methods (SC3, EMEP) and reprehensive learning-based methods (MPSSC, SIMILR, SMSC and SinNLRR). 

- Seurat-[code](https://github.com/satijalab/seurat)
###### 13. Grubman, A., et al., A single-cell atlas of entorhinal cortex from individuals with Alzheimer’s disease reveals cell-type-specific gene expression regulation. Nature neuroscience, 2019. 22(12): p. 2087-2097.
- SNN-Cliq-[code](http://bioinfo.uncc.edu/SNNCliq)
###### 14. Xu C, Su Z. Identification of cell types from single-cell transcriptomes using a novel clustering method[J]. Bioinformatics, 2015, 31(12): 1974-1980.
- SC3-[code](http://bioconductor.org/packages/SC3)
###### 15. Kiselev, V.Y., et al., SC3: consensus clustering of single-cell RNA-seq data. Nat Methods, 2017. 14(5): p. 483-486.
- EMEP-[code](https://github.com/lixt314/EMEP)
###### 16. Li, X., S. Zhang, and K.C. Wong, Single-cell RNA-seq Interpretations using Evolutionary Multiobjective Ensemble Pruning. Bioinformatics, 2018.
- MPSSC-[code](https://github.com/ishspsy/project/tree/master/MPSSC)
###### 17. Park, S. and H. Zhao, Spectral clustering based on learning similarity matrix. Bioinformatics, 2018. 34(12): p. 2069-2076.
- SIMILR-[code](https://github.com/BatzoglouLabSU/SIMLR)
###### 18. Wang, B., et al., SIMLR: A Tool for Large-Scale Genomic Analyses by Multi-Kernel Learning. Proteomics, 2018. 18(2).
- SMSC-[code](https://github.com/Cuteu/SMSC)
###### 19. Qi, R., et al., A spectral clustering with self-weighted multiple kernel learning method for single-cell RNA-seq data. Briefings in Bioinformatics, 2020.
- SinNLRR-[code](https://github.com/zrq0123/SinNLRR)
###### 20. Zheng, R., et al., SinNLRR: a robust subspace clustering method for cell type detection by nonnegative and low rank representation. Bioinformatics, 2019.

## Example Usage:
A demo is provided in [run](https://github.com/lqmmring/MTGDC/blob/main/run.m) file, showing details of data preprocessing, and clustering with MTGDC. Two datasets (`Mouse bladder` and `Mouse retina`) are also given.

### - Preprocessing
A pre-processed single-cell data is accepted, provided that it is log-transformed. 
The input is configured as n cells (rows) by m genes (columns).

### - AIDE
```python
# Load data:
# For small to medium size datasets (up to few hundred thousands of cells)
import pandas as pd
import numpy as np

# Make sure that the final input is n cells (rows) by m genes (cols)
sc_data = pd.read_csv("single_cell_dataset.csv", index_col=0)
sc_data = sc_data.values.astype('float32') # type = np.ndarray

# Configuring AIDE parameters:
from aide import AIDE, AIDEConfig
config = AIDEConfig()
# We may tune the following 4 parameters, but default values are usually sufficient.
config.pretrain_step_num = 1000 # Pretrain step
config.ae_drop_out_rate = 0.4 # Dropout rate
config.alpha = 12.0 # A parameter that determines the portion of AE loss vs MDS encoder loss
config.early_stop = True # Early stop (maximum step number = 20000, minimum step number = 4000)

# Running AIDE:
encoder = AIDE(name = "sc_test", save_folder = "sc_test_folder")
sc_embedding = encoder.fit_transform(sc_data, config=config)

# save embedding
np.savetxt("~/sc_embedding.txt", sc_embedding)
```

### - RPH-kmeans

```python
from rph_kmeans import RPHKMeans
# In the case that n_clusters is already known:
clt = RPHKMeans(n_init=10, n_clusters=10)
clt_labels = clt.fit_predict(sc_embedding)

# In the case that n_clusters is unknown: In order to automatically detect the number of clusters, 
# we implemented a weighted BIC value that determines the optimal k based on 'kneedle' point.

# Important Note: Please set the parameter max_point to a smaller number for small datasets (i.e. less than 2000 cells). 
# The max_point defaults to 2000 cells, and the RPH algorithm stops when the reduced number of cells is below max_point.
# In other words, RPH is not performed if the dataset is smaller than max_point.

max_point = 50 # Defaults to 2000

from rph_kmeans import select_k_with_bic
kmax = 30 # set the maximum number of k to explore
optimal_k, _, _ = select_k_with_bic(sc_embedding, kmax=kmax, point_reducer_kwargs={'max_point':max_point})

clt = RPHKmeans(n_init=10, n_clusters=optimal_k, max_point = max_point) # run RPH-kmeans with optimal_k to get the clustering results

clt_labels = clt.fit_predict(sc_embedding)

# Output results
np.savetxt("~/clt_labels.txt", clt_labels)
```

### - Biological Analysis

```r
library(scAIDE)

# load original data and clustering results:
sc_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # rows = genes, cols = cells
sc_clusters <- read.table("clt_labels.txt")$V1
# Evaluate wilcox rank sum test and log fold change values
eval_gene_markers <- store_markers(gene_expression_data, sc_clusters, threads = 8)
gene_names <- rownames(gene_expression_data)
# returns the list of markers for each cluster, with your specified threshold
cluster_markers_list <- curate_markers(eval_gene_markers, gene_names, wilcox_threshold=0.001, logfc_threshold=1.5)

# Cell type assignment probability according to the markers in the database
# panglao_marker_list: pre-processed list of markers for neural and immune cell types.
# returns a cluster (rows) by cell types (cols) matrix
celltype_prob <- calculate_celltype_prob(cluster_markers_list, panglao_marker_list, type = "jacc")
celltype_id <- rowMaxs(celltype_prob)

# Enrichment probability (based on hypergeometric distribution), this is to be compared with celltype_id to ensure that the number of marker genes detected is of statistical significance.
n_genes <- nrow(gene_expression_data)
# returns a cluster (rows) by cell types (cols) matrix, with p-value entries
enrichment_prob <- calculate_enrichment_prob(cluster_markers_list, panglao_marker_list, n_genes, type = "jacc")

######################################################################
# Visualizing marker genes:
# example marker list:
selected_marker_genes <- c("SOX2", "ALDOC", "CCND2", "OLIG1", "OLIG2")
gene_expression_subset <- gene_expression_data[match(tolower(selected_marker_genes), tolower(rownmaes(gene_expression_data))), ]
# Process the data for plots
processed_markers <- process_marker_expression(gene_expression_subset, sc_clusters)
# Specify the cell type order, and the gene order that you wish to plot
cell_levels <- unique(sc_clusters)
gene_levels <- selected_marker_genes
marker_plot <- plot_marker_expression(processed_markers, gene_levels=gene_levels, cell_levels=cell_levels)
```

### - Example results
The annotated labels for PBMC and Neural datasets are included in the folder 'Predicted labels'. The .RData files include the predicted annotated labels for these datasets. </br>
</br>
The following figures show the results for the PBMC 68k dataset and the 1.3 million neural dataset. 

<p align="center">
  <img src=figures/pbmc.png alt="pbmc" title="pbmc" align="center" height="300">
  <img src= figures/neural.png alt="neural" title="neural" align="center" height="300">
</p>

### - Scalability
The time taken to cluster 1.3 million cells (with roughly 20,000 genes) is less than 30 minutes, using 7GB of memory.

## File Illustration
- `sc_cluster`: Codes and results of clustering experiments using AIDE and RPH-kmeans.
- `baseline`: Codes and results of clustering experiments using baseline tools (e.g. DCA, MAGIC, scScope, scDeepCluster, ...).
- `demo`: A demo of data preprocessing, embedding with AIDE and clustering with RPH-kmeans.
- `scAIDE`: The R package of biological analysis.
- `figures`: Figures of README

## Citation & References

scAIDE: clustering of large-scale single-cell RNA-seq data reveals putative and rare cell types. NAR Genomics and Bioinformatics 2020.

References:
###### 1. Zheng, G. X. et al. (2017) Massively parallel digital transcriptional profiling of single cells. Nature Communications 8, 14049, doi:10.1038/ncomms14049
###### 2. Genomics, X. J. C. B. 1.3 million brain cells from E18 mice. (2017).
###### 3. Franzen, O., Gan, L. M. & Bjorkegren, J. L. M. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. Database (Oxford) 2019.

## Acknowledgment

The authors would like to appreciate the support and guidance from Dr. G.H. Wang and Dr. J. Li. 

## Maintenance

If there's any questions / problems about MTGDC, please feel free to contact Q.M. Liu - cslqm@hit.edu.cn. Thank you!
