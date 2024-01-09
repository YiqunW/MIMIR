Step3: Calculate functional similarities
================
Yiqun Wang

In this example, we will calculate the functional similarities between
enriched genes identified from step 1 in the zebrafish notochord and
hatching gland. Several functional annotation databases will be used,
including Gene Ontology, Reactome, Interpro, and STRING.

``` r
source("../functions/functional_similarity_calculations.R")
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: GOSemSim

    ## 

    ## GOSemSim v2.16.1  For help: https://guangchuangyu.github.io/GOSemSim
    ## 
    ## If you use GOSemSim in published research, please cite:
    ## [36m-[39m Guangchuang Yu. Gene Ontology Semantic Similarity Analysis Using GOSemSim. In: Kidder B. (eds) Stem Cell Transcriptional Networks. Methods in Molecular Biology, 2020, 2117:207-215. Humana, New York, NY. doi:10.1007/978-1-0716-0301-7_11
    ## [36m-[39m Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu, Shengqi Wang. GOSemSim: an R package for measuring semantic similarity among GO terms and gene products Bioinformatics 2010, 26(7):976-978. doi:10.1093/bioinformatics/btq064

    ## Loading required package: GO.db

## Read in enriched genes and background genes

``` r
noto.genes=scan("../example_results/genes.noto.txt",what="character") #genes enriched in the notochord (calculated in step1)
pcp.genes=scan("../example_results/genes.pcp.txt",what="character") #genes enriched in the hatching gland (calculated in step1)
bg.genes=scan("../example_data/genes.bg.2x.txt", what="character") #union of genes enriched in all other cell types
```

## Gather annotation data tables needed

Gene Ontololy, Reactome, Interpro, and HGNC annotate genes with
functional terms or identifiers. The functional terms in these databases
are organized hierarchically. We will calculate the similarities between
the functional terms that are associated with different genes and use
that as the functional similarities between genes. To calculate the
similarities between functional terms in a hierarchical system, we will
use methods for calculating semantic similarity as described in [Yu et
al. 2010](doi:10.1093/bioinformatics/btq064) and their R package
[GOSemSim](https://bioconductor.org/packages/release/bioc/html/GOSemSim.html).
Essentially, 2 types of information are needed: 1. functional annotation
terms associated with each gene; 2. the hierarchical relation between
functional terms. We subsetted the gene to annotation term tables to
contain only genes expressed in the embryo (union of genes enriched in
each cell type, unique(c(noto.genes, pcp.genes, bg.genes))\`).

### Reactome

``` r
# 1. A reactome identifier to zfin symbol table
### This table is constructed from the file downloaded from Reactome website at https://reactome.org/download/current/Ensembl2Reactome.txt
### The downloaded table contains ensemble protein ids, their reactome pathway identifiers and descriptions.
### We mapped the ensemble ids to zfin symbols that match our scRNA-seq data
### The following table contains genes' zfin symbols (column 'gene') and their reactome identifiers (column 'id')
reacAnno <- read.table("../example_data/Reactome2ZfinSymbol.txt", header = 1, stringsAsFactors = F)

### convert reactome ids to that of humans (to facilitate id hierachical relationship mapping, since annotations for human ids are more complete)
reacAnno.orig=reacAnno
reacAnno$id=sub("DRE","HSA",reacAnno$id)

# 2. The parent-child relation between Reactome identifiers
### This table is downloaded from the Reactome website at https://reactome.org/download/current/ReactomePathwaysRelation.txt
reac.p2c <- read.table("../example_data/ReactomePathwaysRelation.txt")
colnames(reac.p2c) <- c("parent","child")
## keep only human (HSA) and zebrafish (DRE) annotations
reac.p2c <- reac.p2c[c(grep("HSA",reac.p2c$parent),grep("DRE",reac.p2c$parent)),]#4216 left
```

### Gene Ontology (GO)

``` r
# 1. A GO term ID to zfin symbol table. An additional column is included to specify whether a term is "biological process", "cellular component", or "molecular function".
### This table can be constructed with various methods. The zebrafish Gene Ontology Data can be downloaded from https://zfin.org/downloads, which contains the relevant information. Alternatively, the biomaRt R package can be used to get GO ids for all genes. 
goAnno <- read.csv("../example_data/go_table_axial_sep9_2021.csv", stringsAsFactors = F)

# 2. The hierarchical relation for GO terms are available from the package AnnotationDbi and will be automatically called in downstream functions
```

### InterPro

``` r
# 1. A table that matches gene symbols with interpro annotations
### This table can be constructed with various methods. The zebrafish ZFIN Marker associations to InterPro protein data can be downloaded from https://zfin.org/downloads. Alternatively, the biomaRt R package can be used to get interpro ids for all genes. 
interpAnno <- read.csv("../example_data/Interpro2ZfinSymbol.csv", stringsAsFactors = F)

# 2. Get parent-child table for interpro ids (relationship file downloaded from interpro site: https://www.ebi.ac.uk/interpro/download/)
### Read and format the downloaded file into a dataframe with 2 columns ("parent" and "child")
fileName <- "../example_data/interpro_ParentChildTreeFile.txt"
conn <- file(fileName,open="r")
interpRelation <-readLines(conn)
close(conn)

## organize parent-child info into a list (list("parent1"=children1, "parent2"=children2, ...))
interp.offspring=list()
level_indicator="--"
term_start="IPR"
while(length(interpRelation)>0){
  p_nodes_ind=which(startsWith(interpRelation,term_start))
  p_nodes_ind=c(p_nodes_ind,length(interpRelation)+1)
  for(i in 1:(length(p_nodes_ind)-1)){
    ind=p_nodes_ind[i]
    ind.next=p_nodes_ind[i+1]
    if((ind.next-ind)>1){
      p_node=unlist(strsplit(interpRelation[ind],"::"))[1]
      p_node=gsub("-","",p_node)
      child_lines=interpRelation[ind:(ind.next-1)]
      child_nodes=unlist(strsplit(child_lines,"-"))
      child_nodes=unlist(strsplit(child_nodes,"::"))
      child_nodes=child_nodes[which(startsWith(child_nodes,"IPR"))]
      interp.offspring[[p_node]]=setdiff(child_nodes,p_node)
    }
  }
  interpRelation=interpRelation[-p_nodes_ind[1:(length(p_nodes_ind)-1)]]
  term_start=paste0(level_indicator,term_start)
}

## additionally format as a table for input to the similarity calculation function
interp.p2c=stack(interp.offspring)
interp.p2c=interp.p2c[,c(2,1)]
colnames(interp.p2c)=c("parent","child")
interp.p2c=as.data.frame(as.matrix(interp.p2c))
```

## Calculate functional similarities between genes

This involves 2 steps: 1. Calculate the information content (IC) of each
annotation term 2. Use the IC and term hierarchical relation to
calculate semantic similarities between genes’ annotated terms Functions
for computing IC and Semantic similarities are adapted from the
[GOSemSim package](https://github.com/YuLab-SMU/GOSemSim) package, which
were designed for working with only GO terms.

### 1. Calculate infomation content (IC) of each annotation term (use all genes in the genome or all genes expressed in the embryo)

This step requires both the parent-child relationship between terms and
the overall frequency of each term. The frequencies of terms will be
calculated from the gene-annotation tables, which should contain not
only the enriched genes, but also background genes (e.g. all genes in
the genome, or all genes expressed in embryo). In our example, we used
the union of genes enriched in each cell type.

``` r
## Reactome
reac.offspring = offspring.list(reac.p2c) ## format parent-child relation dataframe into a list (list("parent1"=children1, "parent2"=children2, ...))
reac.offspring = reac.offspring[which(!is.na(reac.offspring))]
reac.ic <- computeIC.anno(reacAnno, cat=NULL, Offsprings=reac.offspring)
### use IC calculated for human
reac.ic.h=reac.ic[grep("HSA",names(reac.ic))]

## GO
### IC for terms in Molecular Funciton, Cellular Component, and Biological Function will be calculated separately
go.ic=list()
for(cat in c("BP","CC","MF")){
  go.ic[[cat]] <- computeIC.anno(goAnno, cat, Offsprings=NULL) # set Offsprings=NULL to use the GO relations from library(AnnotationDbi)
}


# Interpro
interp.ic <- computeIC.anno(interpAnno, cat=NULL, Offsprings=interp.offspring)
```

#### Example calculations of functional similarities

We can now calculate functional similarities between pairs of genes
based on each functional annotation database. The function
`geneSim.anno()` returns both the similarity score (between 0 and 1), as
well as the annotation terms associated with each gene.

``` r
## intepro 
geneSim.anno("SOX32","XBP1",interpAnno, interp.ic, measure="Jiang", cat=NULL,anno.p2c = interp.p2c,combine = "BMA")
```

    ## $geneSim
    ## [1] 0
    ## 
    ## $GO1
    ## [1] "IPR009071" "IPR036910"
    ## 
    ## $GO2
    ## [1] "IPR004827"

``` r
geneSim.anno("ATF6","XBP1",interpAnno, interp.ic, measure="Jiang", cat=NULL,anno.p2c = interp.p2c,combine = "BMA")
```

    ## $geneSim
    ## [1] 0.6667
    ## 
    ## $GO1
    ## [1] "IPR004827" "IPR029801"
    ## 
    ## $GO2
    ## [1] "IPR004827"

``` r
## go
# cat="MF"
# geneSim.anno("SOX32","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
# geneSim.anno("ATF6","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
# cat="CC"
# geneSim.anno("SOX32","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
# geneSim.anno("ATF6","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
cat="BP"
geneSim.anno("SOX32","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
```

    ## $geneSim
    ## [1] 0.3808
    ## 
    ## $GO1
    ##  [1] "GO:0007507" "GO:0001570" "GO:0001525" "GO:0006355" "GO:0030154"
    ##  [6] "GO:0009653" "GO:0045944" "GO:0001706" "GO:0042663" "GO:0001946"
    ## [11] "GO:0061035" "GO:0003262" "GO:0007492" "GO:0045893" "GO:0042074"
    ## [16] "GO:0035701" "GO:0007493" "GO:0010468" "GO:0001889" "GO:0007498"
    ## [21] "GO:0043534" "GO:0003007" "GO:0007368" "GO:0035050" "GO:0001885"
    ## [26] "GO:0060216"
    ## 
    ## $GO2
    ## [1] "GO:0006355" "GO:0006990" "GO:0035188" "GO:0006366" "GO:0001889"
    ## [6] "GO:0048785"

``` r
geneSim.anno("ATF6","XBP1",goAnno[which(goAnno$category==cat),], go.ic[[cat]], measure="Jiang", cat=cat,combine = "BMA")
```

    ## $geneSim
    ## [1] 0.59
    ## 
    ## $GO1
    ## [1] "GO:0006355" "GO:0045944" "GO:0030968" "GO:0006357" "GO:0001666"
    ## 
    ## $GO2
    ## [1] "GO:0006355" "GO:0006990" "GO:0035188" "GO:0006366" "GO:0001889"
    ## [6] "GO:0048785"

``` r
## reactome
geneSim.anno("SOX32","XBP1",reacAnno, reac.ic.h, measure="Jiang", cat=NULL,anno.p2c = reac.p2c, combine = "BMA")
```

    ## $geneSim
    ## [1] 0
    ## 
    ## $GO1
    ## [1] "R-HSA-3769402"
    ## 
    ## $GO2
    ## [1] "R-HSA-381038" "R-HSA-381070" "R-HSA-381183"

``` r
geneSim.anno("ATF6","XBP1",reacAnno, reac.ic.h, measure="Jiang", cat=NULL,anno.p2c = reac.p2c, combine = "BMA")
```

    ## $geneSim
    ## [1] 0.4888
    ## 
    ## $GO1
    ## [1] "R-HSA-381033"
    ## 
    ## $GO2
    ## [1] "R-HSA-381038" "R-HSA-381070" "R-HSA-381183"

### 2. Calculate pairwise gene similarities based on the annotations

Similarities are calculated for both genes enriched in cell types of
interest and background genes (genes enriched in other cell types).
Similarity between background genes will be used to calculate the priors
when combining annotation similarities from different databases. If only
one database will be used, background genes can be omitted to speed up
calculation.

``` r
## Reactome 
### create an empty similarity matrix
reac.axial.genes=unique(reacAnno$gene)
reac.axial.sim=matrix(0,nrow = length(reac.axial.genes),ncol=length(reac.axial.genes),dimnames = list(reac.axial.genes,reac.axial.genes))
### This step is very slow if running the for loop as is. We recommend subsetting the job and run multiple threads in parallel.
genes=reac.axial.genes
for(i in 1:(length(genes)-1)){
  if(i%%50==0){
     print(paste0(i,"genes calculated"))
  }
  genei=genes[i]
  for(j in (i+1):length(genes)){
    genej=genes[j]
    res=geneSim.anno(genei,genej, reacAnno, reac.ic.h, measure="Jiang", cat=NULL,anno.p2c = reac.p2c, combine = "BMA")
    reac.axial.sim[genei,genej]=res$geneSim
    reac.axial.sim[genej,genei]=res$geneSim
  }
}


## GO
### Create empty similarity matrices
go.axial.genes=unique(goAnno$gene)
go.axial.sim=list()
go.axial.sim[["BP"]]=matrix(0,nrow = length(go.axial.genes),ncol=length(go.axial.genes),dimnames = list(go.axial.genes,go.axial.genes))
go.axial.sim[["MF"]]=matrix(0,nrow = length(go.axial.genes),ncol=length(go.axial.genes),dimnames = list(go.axial.genes,go.axial.genes))
go.axial.sim[["CC"]]=matrix(0,nrow = length(go.axial.genes),ncol=length(go.axial.genes),dimnames = list(go.axial.genes,go.axial.genes))
### Calculate pairwise gene similarities 
### This step is very slow if running the for loop as is. We recommend subsetting the job and run multiple threads in parallel.
genes=go.axial.genes
for(cat in names(go.axial.sim)){
  goAnno.use=goAnno[which(goAnno$category==cat),]
  print(cat)
  for(i in 1:(length(genes)-1)){
    if(i%%50==0){
      print(paste0(i,"genes calculated"))
    }
    genei=genes[i]
    for(j in (i+1):length(genes)){
      genej=genes[j]
      res=geneSim.anno(genei,genej, goAnno.use, go.ic[[cat]], measure="Jiang", cat=cat, combine = "BMA") # hierarchical relation between GO terms will be automatically calculated in function (for non-GO databases, hierarchical relation needs to be supplied to this function)
      go.axial.sim[[cat]][genei,genej]=res$geneSim
      go.axial.sim[[cat]][genej,genei]=res$geneSim
    }
  }
}


## Interpro
### create an empty simlarity matrix
interp.axial.genes=unique(interpAnno$gene)
interp.axial.sim=matrix(0,nrow = length(interp.axial.genes),ncol=length(interp.axial.genes),dimnames = list(interp.axial.genes,interp.axial.genes))
### This step is very slow if running the for loop as is. We recommend subsetting the job and run multiple threads in parallel.
genes=interp.axial.genes
for(i in 1:(length(genes)-1)){
  if(i%%50==0){
     print(paste0(i,"genes calculated"))
  }
  genei=genes[i]
  for(j in (i+1):length(genes)){
    genej=genes[j]
    res=geneSim.anno(genei,genej, interpAnno, interp.ic, measure="Jiang", cat=NULL,anno.p2c = interp.p2c, combine = "BMA")
    interp.axial.sim[genei,genej]=res$geneSim
    interp.axial.sim[genej,genei]=res$geneSim
  }
}
```

## Save functional similarity results and annotation tables used

``` r
## diagonal values are all 0 right now
all.go=list("similarity"=go.axial.sim, "IC"=go.ic, "annotation"=goAnno)
all.reac=list("similarity"=reac.axial.sim, "IC"=reac.ic.h, "annotation"=reacAnno, "parent_child"=reac.p2c)
all.interp=list("similarity"=interp.axial.sim, "IC"=interp.ic, "annotation"=interpAnno, "parent_child"=interp.p2c)
saveRDS(list("GO"=all.go,"Reactome"=all.reac,"Interprot"=all.interp),"../example_results/SemanticSimilarities_Sep2021.rds")
```

## Extract functional similarities from STRING database

### read in the STRING protein-protein interaction (PPI) info file for zebrafish

STRING PPI tables can be downloaded from
<https://string-db.org/cgi/download.pl> (we used
7955.protein.links.full.v11.5.txt.gz). The downloaded table contains
Ensembl protein ids. We converted them into gene symbols that match what
we use in our scRNA-seq data and subsetted the table to contain only the
notochord/hatching gland enriched genes and the background genes.

``` r
str_zf_symb=read.csv("../example_data/string_table_GeneSymbol_sep_2021.csv",stringsAsFactors = F)
## print header
colnames(str_zf_symb)
```

    ##  [1] "gene1"                    "gene2"                   
    ##  [3] "neighborhood"             "neighborhood_transferred"
    ##  [5] "fusion"                   "cooccurence"             
    ##  [7] "homology"                 "coexpression"            
    ##  [9] "coexpression_transferred" "experiments"             
    ## [11] "experiments_transferred"  "database"                
    ## [13] "database_transferred"     "textmining"              
    ## [15] "textmining_transferred"   "combined_score"

The last column (“combined\_score”) is the default combined score from
STRING database.

### Compute STRING PPI scores with selected evidence channels

``` r
## 1. All channels except "neighborhood" and "database".
## The "database" scores in STRING are derived from a number of functional annotation databases including GO and Reactome. We found that the STRING "database" scores are more sparse compared to the annotation similarities we calculated above. For instance, XBP1 and ATF6 has a database score of 0 in STRING, but semantic similarity calculation returns 0.67, 0.59, 0.48 using interpro, GO, and reactome, respectively. We could omit the database channel from STRING and supplement functional similarity scores with the semantic similarities calculated earlier. 
set_noDB=str_comb(str_zf_symb, protein_cols=c("gene1","gene2"), 
                  evi=c("fusion","cooccurence","coexpression", "coexpression_transferred", "experiments", "experiments_transferred",
                        "textmining","textmining_transferred", "homology"))


## 2. Only "experiments" and "textmining"
set_exp_txt=str_comb(str_zf_symb, protein_cols=c("gene1","gene2"), 
                     evi=c("experiments", "experiments_transferred","textmining","textmining_transferred", "homology"))

## gather all combined scores into one dataframe
str_combined=str_zf_symb[,c("gene1","gene2")]
str_combined["combined_score"]=str_zf_symb["combined_score"]/1000
str_combined["combined_noDB"]=set_noDB["combined_score"]
str_combined["combined_exp_txt"]= set_exp_txt["combined_score"]

write.csv(str_combined, "../example_results/STRING_scores.csv",quote = F,row.names = F)
head(str_combined)
```

    ##    gene1  gene2 combined_score combined_noDB combined_exp_txt
    ## 1 CCDC80 COL4A5          0.173     0.1734077            0.154
    ## 2 CCDC80 FSTL1A          0.177     0.1770000            0.041
    ## 3 CCDC80  ECSIT          0.161     0.1620000            0.162
    ## 4 CCDC80 COL5A1          0.472     0.4728117            0.364
    ## 5 CCDC80  ACTN1          0.156     0.1570000            0.157
    ## 6 CCDC80   PAX9          0.434     0.4340000            0.434

### Format pair-wise STRING score dataframe into symmetric similarity matrices

``` r
## only keep the notochord or hatching gland enriched genes
ind1=which(str_combined$gene1%in%union(pcp.genes,noto.genes))
ind2=which(str_combined$gene2%in%union(pcp.genes,noto.genes))
str_axial_ind=intersect(ind1,ind2)
## construct similarity matrix for each combined score (combined_score, combined_noDB, combined_exp_txt)
str_score.matrices=list()
for(chan in setdiff(colnames(str_combined),c("gene1","gene2"))){
  str_score.matrices[[chan]]=str_ppi_2tbl(str_combined[str_axial_ind,],thres = 0,values = chan, g1.col="gene1", g2.col="gene2", method="max")
}
```

### Save STRING similarity matrices

``` r
saveRDS(str_score.matrices, "../example_results/STRING_score_matrices.rds")
```
