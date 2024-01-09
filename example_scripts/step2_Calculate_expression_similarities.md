Step2: Calculate Expression Similarities
================
Yiqun Wang

In this example, we will calculate the expression similarities between
enriched genes identified from step 1 in the zebrafish notochord and
hatching gland. Several similarity metrics will be calculated.

``` r
source("../functions/expression_similarity_calculations.R")
```

## Read in expression matrix for notochord and hatching gland

The expression matrices were extracted from the URD reconstructed
transcriptional trajectories for hatching gland and notochord from 3 to
24 hpf. The trajectory data (URD object DraftTree.rds) can be found in
the example\_data folder. Cell names in the expression matrices are
replaced with their pseudotimes.

``` r
## read in read normalized log scaled expression for notochord (noto) and hatching gland (pcp, prechordal plate)
noto.exp = read.csv("../example_data/notochord_GeneByPseudotime.csv",as.is = T, stringsAsFactors = F,check.names = F,row.names = 1)
pcp.exp = read.csv("../example_data/hatching_gland_GeneByPseudotime.csv",as.is = T, stringsAsFactors = F,check.names = F,row.names = 1)

print("Will calculate expression similarities between:")
```

    ## [1] "Will calculate expression similarities between:"

``` r
print(paste0("1. ", dim(noto.exp)[1], " notochord enriched genes in ", dim(noto.exp)[2], " notochord cells"))
```

    ## [1] "1. 806 notochord enriched genes in 1221 notochord cells"

``` r
print(paste0("2. ",dim(pcp.exp)[1], " hatching gland enriched genes in ", dim(pcp.exp)[2], " hatching gland cells"))
```

    ## [1] "2. 802 hatching gland enriched genes in 1343 hatching gland cells"

``` r
## store expression matrices in a list
exp.data=list()
exp.data[["exp_noto"]]=as.matrix(noto.exp)
exp.data[["exp_pcp"]]=as.matrix(pcp.exp)

## calculate a min-max scaled version of the expression matrices
exp.data$scl_exp_noto=scale_trim(exp.data$exp_noto,min.sub = 0.15,max.scl = 0.01)
exp.data$scl_exp_pcp=scale_trim(exp.data$exp_pcp,min.sub = 0.15,max.scl = 0.01)
```

## Calculate expression similarities (distances)

``` r
## initiate 3D arrays to store results for each cell type
library(abind)
pcp.genes=rownames(pcp.exp)
noto.genes=rownames(noto.exp)
exp.dist.pcp=array(0,dim=c(length(pcp.genes),length(pcp.genes),1),
                        dimnames = list(pcp.genes, pcp.genes, "cos_dist"))
exp.dist.noto=array(0,dim=c(length(noto.genes),length(noto.genes),1),
                        dimnames = list(noto.genes, noto.genes, "cos_dist"))

## Calculate pair-wise expression distances between genes and add to the 3D arrays
## Cosine
cos_dist=cosineDist(exp.data$exp_pcp)
exp.dist.pcp[,,"cos_dist"]=cos_dist[pcp.genes,pcp.genes]
cos_dist=cosineDist(exp.data$exp_noto)
exp.dist.noto[,,"cos_dist"]=cos_dist[noto.genes,noto.genes]

## Soft-cosine
exp.matrix=exp.data$exp_noto
s.blur=simM(pt = as.numeric(colnames(exp.matrix)),sd = 0.01,cutoff = 4)
soft_cos=soft.cos.dist(exp.matrix,s.blur)
soft_cos=soft_cos[noto.genes,noto.genes]
exp.dist.noto=abind(exp.dist.noto,soft_cos,make.names=T)

exp.matrix=exp.data$exp_pcp
s.blur=simM(pt = as.numeric(colnames(exp.matrix)),sd = 0.01,cutoff = 4)
soft_cos=soft.cos.dist(exp.matrix,s.blur)
soft_cos=soft_cos[pcp.genes,pcp.genes]
exp.dist.pcp=abind(exp.dist.pcp,soft_cos,make.names=T)

## Euclidean (min-max scaled)
euc=as.matrix(dist(exp.data$scl_exp_noto,method = 'euclidean'))
euc=euc[noto.genes,noto.genes]
exp.dist.noto=abind(exp.dist.noto,euc,make.names=T)

euc=as.matrix(dist(exp.data$scl_exp_pcp,method = 'euclidean'))
euc=euc[pcp.genes,pcp.genes]
exp.dist.pcp=abind(exp.dist.pcp,euc,make.names=T)

## canberra distance: |u-v|/(|u|+|v|) (min-max scaled)
canb=as.matrix(dist(exp.data$scl_exp_noto,method = 'canberra'))
canb=canb[noto.genes,noto.genes]
exp.dist.noto=abind(exp.dist.noto,canb,make.names=T)

canb=as.matrix(dist(exp.data$scl_exp_pcp,method = 'canberra'))
canb=canb[pcp.genes,pcp.genes]
exp.dist.pcp=abind(exp.dist.pcp,canb,make.names=T)

## Jensen-Shannon divergence
library(philentropy)
exp.matrix=exp.data$exp_noto
jsdiv=JSD(exp.matrix,est.prob = 'empirical')
```

    ## Metric: 'jensen-shannon' using unit: 'log2'; comparing: 806 vectors.

``` r
rownames(jsdiv)=rownames(exp.matrix)
colnames(jsdiv)=rownames(exp.matrix)
jsdiv=jsdiv[noto.genes,noto.genes]
exp.dist.noto=abind(exp.dist.noto,jsdiv,make.names=T)

exp.matrix=exp.data$exp_pcp
jsdiv=JSD(exp.matrix,est.prob = 'empirical')
```

    ## Metric: 'jensen-shannon' using unit: 'log2'; comparing: 802 vectors.

``` r
rownames(jsdiv)=rownames(exp.matrix)
colnames(jsdiv)=rownames(exp.matrix)
jsdiv=jsdiv[pcp.genes,pcp.genes]
exp.dist.pcp=abind(exp.dist.pcp,jsdiv,make.names=T)

## Jensen-Shannon distance
jsdis=sqrt(exp.dist.noto[,,"jsdiv"])
exp.dist.noto=abind(exp.dist.noto,jsdis,make.names=T)

jsdis=sqrt(exp.dist.pcp[,,"jsdiv"])
exp.dist.pcp=abind(exp.dist.pcp,jsdis,make.names=T)
```

## Convert expression distances to similarities

``` r
max.score=0.95 ## similarity scores will be between 0 to 0.95
## Notochord
exp.dist.use=exp.dist.noto
exp.sim.noto=exp.dist.use*0
for (exp.dist in dimnames(exp.dist.use)[[3]]){
  expM=exp.dist.use[,,exp.dist]
  expM=expM*(expM>0)
  expM=1-(expM/max(expM))
  expM=expM*max.score
  exp.sim.noto[,,exp.dist]=expM
}

## Hatching gland
exp.dist.use=exp.dist.pcp
exp.sim.pcp=exp.dist.use*0
for (exp.dist in dimnames(exp.dist.use)[[3]]){
  expM=exp.dist.use[,,exp.dist]
  expM=expM*(expM>0)
  expM=1-(expM/max(expM))
  expM=expM*max.score
  exp.sim.pcp[,,exp.dist]=expM
}
```

### Save expression similarities

``` r
saveRDS(list('noto'=exp.sim.noto, 'pcp'=exp.sim.pcp),"../example_results/exp_similarities.rds")
```