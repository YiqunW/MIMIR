Step4: Cluster with integrated similarities
================
Yiqun Wang

In this example, we will cluster the notochord-enriched genes
(identified in Step1) based on their expression similarities (calculated
in Step2) as well as their functional similarities (calculated in
Step3).

``` r
source("../functions/combine_scores_and_cluster.R")
source("../functions/MIMIR_obj.R")
```

# 1 Load the MIMIR object that contains both expression and functional similarities

Object with expression similarities is constructed in step 2; funcitonal
similarities were calculated in step 3 and added to the same object.

``` r
## object with expression similarity scores and expression matrix:
noto.obj <- readRDS("../example_results/noto_mimir.rds")

## expression similarities available in the object
print(dimnames(noto.obj@exp.sim)[[3]]) 
```

    ## [1] "soft_cosine_sim" "cosine_sim"      "euclidean_sim"   "JS_sim"

``` r
## annotation similarities available in the object
print(dimnames(noto.obj@func.sim)[[3]])
```

    ##  [1] "STRING"                                 
    ##  [2] "STRING_ExpTxt"                          
    ##  [3] "GO_cc"                                  
    ##  [4] "GO_bp"                                  
    ##  [5] "GO_mf"                                  
    ##  [6] "Reactome"                               
    ##  [7] "Interpro"                               
    ##  [8] "STRING+GO+Reactome+Interpro"            
    ##  [9] "STRING_ExpTxt+GO_noCC+Reactome+Interpro"
    ## [10] "GO+Reactome+Interpro"                   
    ## [11] "GO"

# 2 Integrate expression and functional similarities

We can calculate the integrated score for each expression and functional
similarity pair, e.g. cos\_dist and STRING, or jsdiv and
string\_ExpTxt+GO+Reactome+Interpro. We will try different
expression-function combinations and different integration methods and
compare the clustering results derived from them.

``` r
## pairwise combinations of the expression and functional similarities listed in the following vectors will be used
exp_use = c("soft_cosine_sim", "euclidean_sim", "JS_sim")
anno_use = c("STRING+GO+Reactome+Interpro","STRING_ExpTxt+GO_noCC+Reactome+Interpro",
             "GO+Reactome+Interpro", "GO", "STRING")
```

## 2.1 1. Expression OR function: S = 1-(1-Se)(1-Sa)

Integrated score will be high if either expression or functional
similarity is high. For instance, a high expression similarity combined
with a low functional similarity will still give a high score.

``` r
sim <- integrate_all_exp_anno(noto.obj, exp_use=exp_use, anno_use=anno_use, method="OR",add=0, maxscale = 1)
## store in the mimir object
noto.obj@integrated.sim <- sim
## noto.obj@integrated.sim is a 3d array (gene x gene x similarity_method)
dimnames(noto.obj@integrated.sim)[[3]]
```

    ##  [1] "STRING+GO+Reactome+Interpro:OR:soft_cosine_sim"            
    ##  [2] "STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:soft_cosine_sim"
    ##  [3] "GO+Reactome+Interpro:OR:soft_cosine_sim"                   
    ##  [4] "GO:OR:soft_cosine_sim"                                     
    ##  [5] "STRING:OR:soft_cosine_sim"                                 
    ##  [6] "STRING+GO+Reactome+Interpro:OR:euclidean_sim"              
    ##  [7] "STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:euclidean_sim"  
    ##  [8] "GO+Reactome+Interpro:OR:euclidean_sim"                     
    ##  [9] "GO:OR:euclidean_sim"                                       
    ## [10] "STRING:OR:euclidean_sim"                                   
    ## [11] "STRING+GO+Reactome+Interpro:OR:JS_sim"                     
    ## [12] "STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:JS_sim"         
    ## [13] "GO+Reactome+Interpro:OR:JS_sim"                            
    ## [14] "GO:OR:JS_sim"                                              
    ## [15] "STRING:OR:JS_sim"

## 2.2 2. Expression AND function: S = Se\*Sa

Integrated score will be high if both expression and functional
similarities are high.

``` r
sim=integrate_all_exp_anno(noto.obj, exp_use=exp_use, anno_use=anno_use, method="AND",add=0.05, maxscale = 1) # add a small value to avoid 0 product when one of the similarity is 0 but the other is high
## add to the object
noto.obj@integrated.sim <- abind(noto.obj@integrated.sim, sim, along = 3)
```

## 2.3 3. Expression + function: S = Se + Sa

``` r
sim=integrate_all_exp_anno(noto.obj, exp_use=exp_use, anno_use=anno_use, method="+",add=0, maxscale = 1)
## add to the object
noto.obj@integrated.sim <- abind(noto.obj@integrated.sim, sim, along = 3)
```

# 3 Cluster

## 3.1 Cluster the genes using integrated similarities scores

In our `igraph_clus` function, we use functions from the
[igraph](https://r.igraph.org/) package to turn similarity matrices into
gene networks and cluster the genes accordingly, using either built-in
methods in igraph, or leiden algorithm from the
[leiden](https://cran.r-project.org/web/packages/leiden/index.html)
package.

``` r
## Try 3 different clustering methods: louvain, infomap, and leiden. Leiden needs additional resolution and partition method parameters
## try a few resolution parameters for leiden
## trim the similarity matrix to keep only strong connections (can help speed up the calculation and resolve over connections between genes)
res=c(2,3,4,5,6,7)
noto.obj@clusters <- cluster_all(noto.obj@integrated.sim, trim_adj=T, n_kp=120, 
                                 method=c("louvain","infomap","leiden"), leiden_iter=50, leiden_res=res)
```

    ## [1] "Custering using STRING+GO+Reactome+Interpro:OR:soft_cosine_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:soft_cosine_sim"
    ## [1] "Custering using GO+Reactome+Interpro:OR:soft_cosine_sim"
    ## [1] "Custering using GO:OR:soft_cosine_sim"
    ## [1] "Custering using STRING:OR:soft_cosine_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:OR:euclidean_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:euclidean_sim"
    ## [1] "Custering using GO+Reactome+Interpro:OR:euclidean_sim"
    ## [1] "Custering using GO:OR:euclidean_sim"
    ## [1] "Custering using STRING:OR:euclidean_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:OR:JS_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:OR:JS_sim"
    ## [1] "Custering using GO+Reactome+Interpro:OR:JS_sim"
    ## [1] "Custering using GO:OR:JS_sim"
    ## [1] "Custering using STRING:OR:JS_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:AND:soft_cosine_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:AND:soft_cosine_sim"
    ## [1] "Custering using GO+Reactome+Interpro:AND:soft_cosine_sim"
    ## [1] "Custering using GO:AND:soft_cosine_sim"
    ## [1] "Custering using STRING:AND:soft_cosine_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:AND:euclidean_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:AND:euclidean_sim"
    ## [1] "Custering using GO+Reactome+Interpro:AND:euclidean_sim"
    ## [1] "Custering using GO:AND:euclidean_sim"
    ## [1] "Custering using STRING:AND:euclidean_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:AND:JS_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:AND:JS_sim"
    ## [1] "Custering using GO+Reactome+Interpro:AND:JS_sim"
    ## [1] "Custering using GO:AND:JS_sim"
    ## [1] "Custering using STRING:AND:JS_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:+:soft_cosine_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:+:soft_cosine_sim"
    ## [1] "Custering using GO+Reactome+Interpro:+:soft_cosine_sim"
    ## [1] "Custering using GO:+:soft_cosine_sim"
    ## [1] "Custering using STRING:+:soft_cosine_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:+:euclidean_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:+:euclidean_sim"
    ## [1] "Custering using GO+Reactome+Interpro:+:euclidean_sim"
    ## [1] "Custering using GO:+:euclidean_sim"
    ## [1] "Custering using STRING:+:euclidean_sim"
    ## [1] "Custering using STRING+GO+Reactome+Interpro:+:JS_sim"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro:+:JS_sim"
    ## [1] "Custering using GO+Reactome+Interpro:+:JS_sim"
    ## [1] "Custering using GO:+:JS_sim"
    ## [1] "Custering using STRING:+:JS_sim"

## 3.2 Cluster with only expression or only functional similarities too for comparison

``` r
clus <- cluster_all(noto.obj@exp.sim[,,exp_use], trim_adj=T, n_kp=120, 
                                 method=c("louvain","infomap","leiden"), leiden_iter=50, leiden_res=res)
```

    ## [1] "Custering using soft_cosine_sim"
    ## [1] "Custering using euclidean_sim"
    ## [1] "Custering using JS_sim"

``` r
noto.obj@clusters = append(noto.obj@clusters, clus)

clus <- cluster_all(noto.obj@func.sim[,,anno_use], trim_adj=T, n_kp=120, 
                                 method=c("louvain","infomap","leiden"), leiden_iter=50, leiden_res=res)
```

    ## [1] "Custering using STRING+GO+Reactome+Interpro"
    ## [1] "Custering using STRING_ExpTxt+GO_noCC+Reactome+Interpro"
    ## [1] "Custering using GO+Reactome+Interpro"
    ## [1] "Custering using GO"
    ## [1] "Custering using STRING"

``` r
noto.obj@clusters = append(noto.obj@clusters, clus)
```

### 3.2.1 Save the object with clustering results

``` r
saveRDS(noto.obj, "../example_results/noto_mimir.rds")
```

### 3.2.2 Check how cluster results look

``` r
noto.obj=readRDS("../example_results/noto_mimir.rds")
head(noto.obj@clusters$`GO+Reactome+Interpro:AND:euclidean_sim`) 
```

    ##         louvain infomap leiden_2_RBConfigurationVertex
    ## ABHD15A       5       3                             11
    ## ACBD3        10       2                              2
    ## ACKR4B        3       7                              3
    ## ACSBG2        3      10                              3
    ## ACTB2         1      15                              4
    ## ACVR1BA       5       3                              7
    ##         leiden_3_RBConfigurationVertex leiden_4_RBConfigurationVertex
    ## ABHD15A                             15                             11
    ## ACBD3                               19                             25
    ## ACKR4B                              14                             14
    ## ACSBG2                               1                              1
    ## ACTB2                                9                             13
    ## ACVR1BA                              4                              8
    ##         leiden_5_RBConfigurationVertex leiden_6_RBConfigurationVertex
    ## ABHD15A                             23                             24
    ## ACBD3                               32                             36
    ## ACKR4B                              15                             12
    ## ACSBG2                               8                              6
    ## ACTB2                               19                             20
    ## ACVR1BA                             11                             14
    ##         leiden_7_RBConfigurationVertex
    ## ABHD15A                             27
    ## ACBD3                               40
    ## ACKR4B                              15
    ## ACSBG2                              11
    ## ACTB2                               14
    ## ACVR1BA                             19

*Each column contains the cluster id of genes for a specific clustering
method.*
