Step4: Cluster with combined similarities
================
Yiqun Wang

In this example, we will cluster the notochord-enriched genes
(identified in Step1) based on their expression similarities (calculated
in Step2) as well as their functional similarities (calculated in
Step3).

``` r
source("../functions/combine_scores_and_cluster.R")
```

    ## Loading required package: abind

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## Loading required package: leiden

    ## conda environment r-reticulate installed

    ## python modules igraph and leidenalg installed

## Load previously calculated expression and functional similarities

``` r
## Expression similarity scores:
exp_sim = readRDS("../example_results/exp_similarities.rds")
exp_sim=exp_sim$noto
print(dimnames(exp_sim)[[3]]) ## distances were already converted to similarities dispite the names
```

    ## [1] "cos_dist" "soft_cos" "euc"      "canb"     "jsdiv"    "jsdis"

``` r
## Functional similarity scores:
anno_sim = readRDS("../example_results/notochord_functional_similarities.RDS")
print(dimnames(anno_sim)[[3]])
```

    ##  [1] "string_noDB+GO+Reactome+Interpro"     
    ##  [2] "string_noDB+GO_noCC+Reactome+Interpro"
    ##  [3] "string_ExpTxt+GO+Reactome+Interpro"   
    ##  [4] "GO+Reactome+Interpro"                 
    ##  [5] "GO"                                   
    ##  [6] "STRING"                               
    ##  [7] "STRING_noDB"                          
    ##  [8] "STRING_ExpTxt"                        
    ##  [9] "GO_cc"                                
    ## [10] "GO_bp"                                
    ## [11] "GO_mf"                                
    ## [12] "Reactome"                             
    ## [13] "Interpro"

## Integrate expression and functional similarities

We can calculate the integrated score for each expression and functional
similarity pair, e.g. cos\_dist and STRING, or jsdiv and
string\_ExpTxt+GO+Reactome+Interpro. We will try different
expression-function combinations and different integration methods and
compare the clustering results derived from them.

``` r
## will calculate combinations of the following similarities
exp_use = c("soft_cos", "euc", "canb", "jsdis")
anno_use = c("string_noDB+GO+Reactome+Interpro", "string_noDB+GO_noCC+Reactome+Interpro", "string_ExpTxt+GO+Reactome+Interpro",
             "GO+Reactome+Interpro", "GO", "STRING", "Reactome")
```

### 1. Expression OR function: S = 1-(1-Se)(1-Sa)

Integrated score will be high if either expression or functional
similarity is high. For instance, a high expression similarity combined
with a low functional similarity will still give a high score.

``` r
simOR=integrate_all_exp_anno(exp_sim, anno_sim, exp_use=exp_use, anno_use=anno_use, method="OR",add=0, maxscale = 1)
## simOR is a 3d array
dim(simOR)
```

    ## [1] 806 806  28

``` r
dimnames(simOR)[[3]]
```

    ##  [1] "string_noDB+GO+Reactome+Interpro:OR:soft_cos"     
    ##  [2] "string_noDB+GO_noCC+Reactome+Interpro:OR:soft_cos"
    ##  [3] "string_ExpTxt+GO+Reactome+Interpro:OR:soft_cos"   
    ##  [4] "GO+Reactome+Interpro:OR:soft_cos"                 
    ##  [5] "GO:OR:soft_cos"                                   
    ##  [6] "STRING:OR:soft_cos"                               
    ##  [7] "Reactome:OR:soft_cos"                             
    ##  [8] "string_noDB+GO+Reactome+Interpro:OR:euc"          
    ##  [9] "string_noDB+GO_noCC+Reactome+Interpro:OR:euc"     
    ## [10] "string_ExpTxt+GO+Reactome+Interpro:OR:euc"        
    ## [11] "GO+Reactome+Interpro:OR:euc"                      
    ## [12] "GO:OR:euc"                                        
    ## [13] "STRING:OR:euc"                                    
    ## [14] "Reactome:OR:euc"                                  
    ## [15] "string_noDB+GO+Reactome+Interpro:OR:canb"         
    ## [16] "string_noDB+GO_noCC+Reactome+Interpro:OR:canb"    
    ## [17] "string_ExpTxt+GO+Reactome+Interpro:OR:canb"       
    ## [18] "GO+Reactome+Interpro:OR:canb"                     
    ## [19] "GO:OR:canb"                                       
    ## [20] "STRING:OR:canb"                                   
    ## [21] "Reactome:OR:canb"                                 
    ## [22] "string_noDB+GO+Reactome+Interpro:OR:jsdis"        
    ## [23] "string_noDB+GO_noCC+Reactome+Interpro:OR:jsdis"   
    ## [24] "string_ExpTxt+GO+Reactome+Interpro:OR:jsdis"      
    ## [25] "GO+Reactome+Interpro:OR:jsdis"                    
    ## [26] "GO:OR:jsdis"                                      
    ## [27] "STRING:OR:jsdis"                                  
    ## [28] "Reactome:OR:jsdis"

``` r
saveRDS(simOR,"../example_results/notochord_OR_integrated_sim.rds", compress = "xz")
```

### 2. Expression AND function: S = Se\*Sa

Integrated score will be high if both expression and functional
similarities are high.

``` r
simAND=integrate_all_exp_anno(exp_sim, anno_sim, exp_use=exp_use, anno_use=anno_use, method="AND",add=0.05, maxscale = 1)
## simOR is a 3d array
dimnames(simAND)[[3]]
```

    ##  [1] "string_noDB+GO+Reactome+Interpro:AND:soft_cos"     
    ##  [2] "string_noDB+GO_noCC+Reactome+Interpro:AND:soft_cos"
    ##  [3] "string_ExpTxt+GO+Reactome+Interpro:AND:soft_cos"   
    ##  [4] "GO+Reactome+Interpro:AND:soft_cos"                 
    ##  [5] "GO:AND:soft_cos"                                   
    ##  [6] "STRING:AND:soft_cos"                               
    ##  [7] "Reactome:AND:soft_cos"                             
    ##  [8] "string_noDB+GO+Reactome+Interpro:AND:euc"          
    ##  [9] "string_noDB+GO_noCC+Reactome+Interpro:AND:euc"     
    ## [10] "string_ExpTxt+GO+Reactome+Interpro:AND:euc"        
    ## [11] "GO+Reactome+Interpro:AND:euc"                      
    ## [12] "GO:AND:euc"                                        
    ## [13] "STRING:AND:euc"                                    
    ## [14] "Reactome:AND:euc"                                  
    ## [15] "string_noDB+GO+Reactome+Interpro:AND:canb"         
    ## [16] "string_noDB+GO_noCC+Reactome+Interpro:AND:canb"    
    ## [17] "string_ExpTxt+GO+Reactome+Interpro:AND:canb"       
    ## [18] "GO+Reactome+Interpro:AND:canb"                     
    ## [19] "GO:AND:canb"                                       
    ## [20] "STRING:AND:canb"                                   
    ## [21] "Reactome:AND:canb"                                 
    ## [22] "string_noDB+GO+Reactome+Interpro:AND:jsdis"        
    ## [23] "string_noDB+GO_noCC+Reactome+Interpro:AND:jsdis"   
    ## [24] "string_ExpTxt+GO+Reactome+Interpro:AND:jsdis"      
    ## [25] "GO+Reactome+Interpro:AND:jsdis"                    
    ## [26] "GO:AND:jsdis"                                      
    ## [27] "STRING:AND:jsdis"                                  
    ## [28] "Reactome:AND:jsdis"

``` r
saveRDS(simAND,"../example_results/notochord_AND_integrated_sim.rds", compress = "xz")
```

### 3. Expression + function: S = Se + Sa

``` r
simPlus=integrate_all_exp_anno(exp_sim, anno_sim, exp_use=exp_use, anno_use=anno_use, method="+",add=0, maxscale = 1)
## simOR is a 3d array
dimnames(simPlus)[[3]]
```

    ##  [1] "string_noDB+GO+Reactome+Interpro:+:soft_cos"     
    ##  [2] "string_noDB+GO_noCC+Reactome+Interpro:+:soft_cos"
    ##  [3] "string_ExpTxt+GO+Reactome+Interpro:+:soft_cos"   
    ##  [4] "GO+Reactome+Interpro:+:soft_cos"                 
    ##  [5] "GO:+:soft_cos"                                   
    ##  [6] "STRING:+:soft_cos"                               
    ##  [7] "Reactome:+:soft_cos"                             
    ##  [8] "string_noDB+GO+Reactome+Interpro:+:euc"          
    ##  [9] "string_noDB+GO_noCC+Reactome+Interpro:+:euc"     
    ## [10] "string_ExpTxt+GO+Reactome+Interpro:+:euc"        
    ## [11] "GO+Reactome+Interpro:+:euc"                      
    ## [12] "GO:+:euc"                                        
    ## [13] "STRING:+:euc"                                    
    ## [14] "Reactome:+:euc"                                  
    ## [15] "string_noDB+GO+Reactome+Interpro:+:canb"         
    ## [16] "string_noDB+GO_noCC+Reactome+Interpro:+:canb"    
    ## [17] "string_ExpTxt+GO+Reactome+Interpro:+:canb"       
    ## [18] "GO+Reactome+Interpro:+:canb"                     
    ## [19] "GO:+:canb"                                       
    ## [20] "STRING:+:canb"                                   
    ## [21] "Reactome:+:canb"                                 
    ## [22] "string_noDB+GO+Reactome+Interpro:+:jsdis"        
    ## [23] "string_noDB+GO_noCC+Reactome+Interpro:+:jsdis"   
    ## [24] "string_ExpTxt+GO+Reactome+Interpro:+:jsdis"      
    ## [25] "GO+Reactome+Interpro:+:jsdis"                    
    ## [26] "GO:+:jsdis"                                      
    ## [27] "STRING:+:jsdis"                                  
    ## [28] "Reactome:+:jsdis"

## Cluster the genes using integrated similarities scores

In our `igraph_clus` function, we use functions from the
[igraph](https://r.igraph.org/) package to turn similarity matrices into
gene networks and cluster the genes accordingly, using either built-in
methods in igraph, or leiden algorithm from the
[leiden](https://cran.r-project.org/web/packages/leiden/index.html)
package.

``` r
## Initiate a list object to store cluster results from different similarity scores
clus_noto=list()

## Try 3 different clustering methods: louvain, infomap, and leiden. Leiden needs additional resolution and partition method parameters
## try a few resolution parameters for leiden
res=c(2,2.5,3,3.5,4,5,7)

## Cluster
## Clustering with different scores and clustering methods can be parallelized to speed up the calculation.
### 1. with OR integrated scores. 
for(s in dimnames(simOR)[[3]]){
  print(s)
  sim_use=trim_adj(simOR[,,s], n_kp=120) # trim the adjacency matrix such that each gene keeps only their top 80 connections. This will speed up clustering and help discretize the overly connected gene network
  clus_noto[[s]]=igraph_cls(sim_use, method=c("louvain","infomap","leiden"),leiden_res=res, 
                                leiden_iter = 50, leiden_par=c("RBConfigurationVertexPartition"), seed=1)
}

### 2. with AND integrated scores. 
for(s in dimnames(simAND)[[3]]){
  print(s)
  sim_use=trim_adj(simAND[,,s], n_kp=120)
  clus_noto[[s]]=igraph_cls(sim_use, method=c("louvain","infomap","leiden"),leiden_res=res,
                                leiden_iter = 50, leiden_par=c("RBConfigurationVertexPartition"), seed=1)
}

### 3. with + integrated scores
for(s in dimnames(simPlus)[[3]]){
  print(s)
  sim_use=trim_adj(simPlus[,,s], n_kp = 120)
  clus_noto[[s]]=igraph_cls(sim_use, method=c("louvain","infomap","leiden"),leiden_res=res,
                                leiden_iter = 50, leiden_par=c("RBConfigurationVertexPartition"), seed=1)
}
```

### Cluster with only expression or only functional similarities too for comparison

``` r
for(s in dimnames(exp_sim)[[3]]){
  print(s)
  sim_use=trim_adj(exp_sim[,,s], n_kp = 120)
  clus_noto[[s]]=igraph_cls(sim_use, method=c("louvain","infomap","leiden"),leiden_res=res,
                                leiden_iter = 50, leiden_par=c("RBConfigurationVertexPartition"), seed=1)
}

for(s in dimnames(anno_sim)[[3]]){
  print(s)
  sim_use=trim_adj(anno_sim[,,s], n_kp = 120)
  clus_noto[[s]]=igraph_cls(sim_use, method=c("louvain","infomap","leiden"),leiden_res=res,
                                leiden_iter = 50, leiden_par=c("RBConfigurationVertexPartition"), seed=1)
}
```

### Save all clustering results

``` r
saveRDS(clus_noto, "../example_results/notochord_all_clusters.rds")
```

### Check how cluster results look

``` r
clus_noto=readRDS("../example_results/notochord_all_clusters.rds")
head(clus_noto$`GO:OR:soft_cos`) ## each matrix contains cluster ids
```

    ##         louvain infomap leiden_2_RBConfigurationVertex
    ## ABHD15A       1       3                              1
    ## ACBD3         5       1                              4
    ## ACKR4B        1       3                              1
    ## ACSBG2        2       6                              3
    ## ACTB2         5       1                              4
    ## ACVR1BA       2       7                              3
    ##         leiden_2.5_RBConfigurationVertex leiden_3_RBConfigurationVertex
    ## ABHD15A                                1                              1
    ## ACBD3                                  3                             15
    ## ACKR4B                                 1                              1
    ## ACSBG2                                 8                             13
    ## ACTB2                                  3                              3
    ## ACVR1BA                                9                              8
    ##         leiden_3.5_RBConfigurationVertex leiden_4_RBConfigurationVertex
    ## ABHD15A                                1                              1
    ## ACBD3                                 17                             17
    ## ACKR4B                                 1                              1
    ## ACSBG2                                11                             10
    ## ACTB2                                  3                              5
    ## ACVR1BA                                8                              7
    ##         leiden_5_RBConfigurationVertex leiden_7_RBConfigurationVertex
    ## ABHD15A                              1                              2
    ## ACBD3                               18                             29
    ## ACKR4B                               1                              1
    ## ACSBG2                              11                             11
    ## ACTB2                               10                             85
    ## ACVR1BA                              5                              7
