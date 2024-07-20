# MIMIR
**MIMIR** is an approach for **M**odule **I**dentification via **M**ulti-source **I**ntegration for **R**econstructing Differentiation. Briefly, this approach groups genes into modules using both gene expression similarities from **single-cell RNA sequencing (scRNA-seq) data and functional similarities from functional annotation databases**. Each gene module contains genes that likely work together to drive a biological process during cell differentiation. MIMIR is part of our effort to extract biological insights from single-cell transcriptomic atlases and trajectories to understand how differentiation is regulated.

# Usage
Scripts in [functions](https://github.com/YiqunW/MIMIR/tree/main/functions) folder include all functions required for running MIMIR. To use these functions, download the R files in the folder and load them into your own script using the `source()` function. We recommend going through our examples in the [example_scripts](https://github.com/YiqunW/MIMIR/tree/main/example_scripts) folder to get started. To run the example scripts on your local machine, files in [example_data](https://github.com/YiqunW/MIMIR/tree/main/example_data) will be needed. You can download individual scripts and data files required, or download the entire repository (`Code` >> `Download ZIP`, or run `git clone https://github.com/YiqunW/MIMIR.git` in your terminal). 

# Steps
MIMIR involves 5 major steps:

### Step 1: Identify enriched genes during differentiation
Genes that are more highly expressed in the cell type of interest compared to background cell types can be identified using appropriate differential expression tests. In our [example](https://github.com/YiqunW/MIMIR/blob/main/example_scripts/step1_Identify_Enriched_Genes.md), we analysized timecourse scRNA-seq data for two cell types: the zebrafish notochord and hatching gland. **The output from this step required for the next steps are**:
1. A list of cell type enriched genes (e.g. [genes.noto.txt](https://github.com/YiqunW/MIMIR/blob/main/example_results/genes.noto.txt))
2. The gene expression matrix (normalized and in log space) with enriched genes as rows and pseudotime as columns (e.g. [notochord_GeneByPseudotime.csv](https://github.com/YiqunW/MIMIR/blob/main/example_data/notochord_GeneByPseudotime.csv))


### Step 2: Calculate expression similarities between enriched genes
First, a MIMIR object is created with the enriched genes by pseudotime expression matrix from the last step:
```r
noto.obj <- createMIMIR(name="Notochord", exp.data=as.matrix(noto.exp))
```
Then, expression similarities between genes are calculated:
```r
noto.obj@exp.dis <- exp_dis(x=noto.obj@exp.data, methods=c("cosine", "soft_cosine", "euclidean", "JS"))
noto.obj@exp.sim <- all_dist_to_sim(noto.obj@exp.dis, max.score = 0.95)
```
A range of expression similarity/ distance metrics can be used (e.g. Euclidean, (soft) cosine, Jensen-Shannon). See our [example](https://github.com/YiqunW/MIMIR/blob/main/example_scripts/step2_Calculate_expression_similarities.md). The resulting object (`noto.obj`) will be needed in later steps. 


### Step 3: Calculate functional similarities between enriched genes
Annotation data downloaded from the functional databases is needed in this step. In our [example](https://github.com/YiqunW/MIMIR/blob/main/example_scripts/step3_Calculate_functional_similarities.md), we used the following functiobal databases: [Gene Ontology](https://geneontology.org/), [Reactome](https://reactome.org/), [InterPro](https://www.ebi.ac.uk/interpro/), and [STRING database](https://string-db.org/cgi/download?sessionId=bykC2Can3gR6). The first three contain annotation terms that are hierarchically organized. We calculate the semantic similarities between genes based on each database:
```r
## reacAnno and reac.p2c are tables downloaded from reactome site
reactome = createAnno(db='Reactome', gene.anno = reacAnno, genes = noto.genes, bg.genes = background.genes,  hierarchy.df=reac.p2c)
reactome <- computeIC.anno(reactome) ##compute information content of annotations
reactome <- func_sim(reactome, genes=reactome@genes) ##compute pair-wise functional similarities between genes
```
The [STRING database](https://string-db.org/cgi/download?sessionId=bykC2Can3gR6) provides protein-protein interaction scores between genes, which can be used directly as similarity scores. Scores calculated from different databases are combined into one functional similarity score per gene pair:
```r
## collect similarities from different databases into a 3d array
anno.sim.noto <- stack_func_sim(list("String"=str.sim.matrix, "Reactome"=reactome@similarity, genes_use=noto.genes)
## combine similartities from different databases and add to the 3d array
anno.sim.noto = add_combined_scores(anno.sim=anno.sim.noto, 
                                    how.to=list("STRING+Reactome" = c("STRING", "Reactome")), add=T)
```
Add annotation similarity to the MIMIR object created in Step 2:
```r
noto.obj <- add_to_MIMIR(noto.obj,anno.sim.noto) ## add functional similarities
noto.obj <- add_anno_to_MIMIR(noto.obj, reactome) ## add gene-annotation tables
```


### Step 4: Cluster enriched genes using combined expression and functional similarities
In this step, expression and functional similarities are integrated. A few methods for similarity integration and clustering methods were tested in our [example](https://github.com/YiqunW/MIMIR/blob/main/example_scripts/step4_Cluster_genes_with_integrated_similarities.md). Expression and functional similarities can be integrated by multiplication (AND), summation (+), and adjusted combined probability (OR):
```r
sim <- integrate_all_exp_anno(noto.obj, method="AND") ## can change method to "+" or "OR"
noto.obj@integrated.sim <- sim
```
To cluster the genes with integrated similarities with a few clustering algorithms and parameters:
```r
noto.obj@clusters <- cluster_all(noto.obj@integrated.sim, method=c("louvain","infomap","leiden"), leiden_iter=50, leiden_res=c(2,4,6,8))
```
The object with clustering results will be used for the next step.


### Step 5: Annotate and curate the resulting gene modules
The previous steps generated a multitude of clustering results with different similarities and parameters tested. In this step, we will choose the best clustering result to serve as the basis of the gene modules (see [example](https://github.com/YiqunW/MIMIR/blob/main/example_scripts/step5_Check_clustering_results.md)). To do this, we inspect several metrics for each clustering result, such as the number and size of clusters, the mean within-cluster gene correlation, and more. When manually curated modules (with a subset of genes) or additional gold-standard gene annotations are available, they can be used to assess the clustering results too.

To calculate quality measures for clustering results:
```r
# module_subset is a table of manually curated gene modules that included a subset of noto.genes
# msigdb_c2 is a gene-annotation table obtained from the MSigDB database
noto.obj <- assess_result(noto.obj, manual_modules=module_subset, external_db_tbl=msigdb_c2)
```

To visualize the calculated quality measures:
```r
plotly_scatter(noto.obj@cluster_metrics, x="n_cluster", y="max_cluster",color="similarity_mode", hover_text = "similarity")
plotly_scatter(metric_tbl, x="pct_gene_in_cluster", y="AMI",color="similarity_mode") 
```

To pick a clustering result and inspect the expression and functional coherence in each cluster:
```r
metric_tbl=noto.obj@cluster_metrics
# pick the result with the highest agreement with the manual curation:
metric_use=metric_use[order(metric_use[["AMI"]], decreasing = T),][1,] 

# add the chosen clustering result to @gene.info slot
noto.obj<-add_cluster_result(noto.obj, metric_use)

# plot gene expression and the top 3 enriched annotations in each cluster
cluster_use=colnames(noto.obj@gene.info)[dim(noto.obj@gene.info)[2]]
cluster_plots = plotClusters(noto.obj, cluster_use=cluster_use, exp_use="smoothed.exp", save_pdf="../example_results/cluster_plots1.pdf") 
```
