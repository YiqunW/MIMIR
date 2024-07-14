Step3: Calculate functional similarities
================
Yiqun Wang

In this example, we will calculate the functional similarities between
enriched genes identified from step 1 in the zebrafish notochord and
hatching gland. Several functional annotation databases will be used,
including Gene Ontology, Reactome, Interpro, and STRING.

Relevant information from the annotation databases will need to be
downloaded as the input for this step.

Functional similarities between enriched genes will be calculated based
on information from each database individually, then these similarities
will be integrated into a combined score. The individual and combined
functional similarities will be organized as a 3d array (gene x gene x
database) and saved into the MIMIR object created in the last step,
which contains expression matrix and expression similarity array.

``` r
source("../functions/functional_similarity_calculations.R")
source("../functions/MIMIR_obj.R")
```

# 1 Read in tissue enriched and background genes

``` r
noto.genes=scan("../example_results/genes.noto.txt",what="character") #genes enriched in the notochord (calculated in step1). 
pcp.genes=scan("../example_results/genes.pcp.txt",what="character") #genes enriched in the hatching gland (calculated in step1). 
bg.genes=scan("../example_data/genes.bg.2x.txt", what="character") #union of genes enriched in all other cell types
```

# 2 Load annotation data tables needed

Gene Ontololy, Reactome, Interpro, and HGNC annotate genes with
functional terms or identifiers. The functional terms in these databases
are organized hierarchically. Essentially, 2 types of information are
needed in this notebook: 1. functional annotation terms associated with
each gene; 2. the hierarchical relation between functional terms.

The data tables loaded below were subsetted to contain annotations
associated with only the genes expressed in the embryo (union of genes
enriched in each cell type, i.e., the union of `noto.genes`,
`pcp.genes`, and `bg.genes`.

## 2.1 Reactome

### 2.1.1 load the reactome identifier to zfin symbol table

This table is constructed from the file downloaded from Reactome website
at <https://reactome.org/download/current/Ensembl2Reactome.txt>

The downloaded table contains ensemble protein ids, their reactome
pathway identifiers and descriptions. We mapped the ensemble ids to zfin
symbols that match our scRNA-seq data. `reacAnno` below contains genes’
zfin symbols (column `gene`) and their reactome identifiers (column
`id`)

``` r
reacAnno <- read.table("../example_data/Reactome2ZfinSymbol.txt", header = 1, stringsAsFactors = F,sep = "\t")
### convert reactome ids to that of humans (to facilitate id hierachical relationship mapping, since annotations for human ids are more complete)
reacAnno.orig=reacAnno
reacAnno$id=sub("DRE","HSA",reacAnno$id)
head(reacAnno)
```

    ##              id   gene                                              description
    ## 1 R-HSA-6811434  STX18            COPI-dependent Golgi-to-ER retrograde traffic
    ## 2  R-HSA-450513  ZFP36 Tristetraprolin (TTP, ZFP36) binds and destabilizes mRNA
    ## 3 R-HSA-1483248 PI4K2B                     Synthesis of PIPs at the ER membrane
    ## 4 R-HSA-1660499 PI4K2B                 Synthesis of PIPs at the plasma membrane
    ## 5 R-HSA-1660514 PI4K2B                  Synthesis of PIPs at the Golgi membrane
    ## 6 R-HSA-1660516 PI4K2B         Synthesis of PIPs at the early endosome membrane

### 2.1.2 load the parent-child relation between Reactome identifiers

This table is downloaded from the Reactome website at
<https://reactome.org/download/current/ReactomePathwaysRelation.txt>

`reac.p2c` below contains two columns: a `parent` column contains parent
reactome annotaion ids, and a `child` column contains child id for the
parent. If a parent id has multiple children terms, it will occupy
multiple rows.

``` r
reac.p2c <- read.table("../example_data/ReactomePathwaysRelation.txt")
colnames(reac.p2c) <- c("parent","child")
## keep human (HSA) annotations
reac.p2c <- reac.p2c[grep("HSA",reac.p2c$parent),]#2555 left
head(reac.p2c)
```

    ##            parent         child
    ## 9978 R-HSA-109581  R-HSA-109606
    ## 9979 R-HSA-109581  R-HSA-169911
    ## 9980 R-HSA-109581 R-HSA-5357769
    ## 9981 R-HSA-109581   R-HSA-75153
    ## 9982 R-HSA-109582  R-HSA-140877
    ## 9983 R-HSA-109582  R-HSA-202733

### 2.1.3 Create an Anno object to store reactome annotation information

``` r
reactome = createAnno(db='Reactome', gene.anno = reacAnno, gene.col='gene', id.col='id', desp.col = 'description',
                      genes = union(noto.genes, pcp.genes), bg.genes = bg.genes, 
                      hierarchy.df=reac.p2c, parent.col=1, child.col=2)
```

    ## 592 input genes are not found in the gene.anno table

    ## 3818 input background genes are not found in the gene.anno table

## 2.2 Gene Ontology (GO)

``` r
# 1. A GO term ID to zfin symbol table. An additional column is included to specify whether a term is "biological process", "cellular component", or "molecular function".
### This table can be constructed with various methods. The zebrafish Gene Ontology Data can be downloaded from https://zfin.org/downloads, which contains the relevant information. Alternatively, the biomaRt R package can be used to get GO ids for all genes. 
goAnno <- read.csv("../example_data/go_table_axial_sep9_2021.csv", stringsAsFactors = F)

# 2. The hierarchical relation for GO terms are available from the package AnnotationDbi and will be automatically called in downstream functions
# inspect GO table
head(goAnno)
```

    ##           id       gene category
    ## 1 GO:0031105      SEPT3       CC
    ## 2 GO:0016021       PCNX       CC
    ## 3 GO:0019215        NES       MF
    ## 4 GO:0006888 ZGC:103697       BP
    ## 5 GO:0046872      NRP2A       MF
    ## 6 GO:0046872 ZGC:113886       MF
    ##                                                 description
    ## 1                                            septin complex
    ## 2                            integral component of membrane
    ## 3                             intermediate filament binding
    ## 4 endoplasmic reticulum to Golgi vesicle-mediated transport
    ## 5                                         metal ion binding
    ## 6                                         metal ion binding

``` r
# create Anno object
go = createAnno(db='GO', gene.anno = goAnno, gene.col='gene', id.col='id', cat.col='category', desp.col = 'description',
                      genes = union(noto.genes, pcp.genes), bg.genes = bg.genes)
```

    ## 34 input genes are not found in the gene.anno table

    ## 333 input background genes are not found in the gene.anno table

## 2.3 InterPro

``` r
# 1. A table that matches gene symbols with interpro annotations
### This table can be constructed with various methods. The zebrafish ZFIN Marker associations to InterPro protein data can be downloaded from https://zfin.org/downloads. Alternatively, the biomaRt R package can be used to get interpro ids for all genes. 
interpAnno <- read.csv("../example_data/Interpro2ZfinSymbol.csv", stringsAsFactors = F)

# 2. Get parent-child table for interpro ids (relationship file downloaded from interpro site: https://www.ebi.ac.uk/interpro/download/)
### Read and format the downloaded file into a dataframe with 2 columns ("parent" and "child")
interp.p2c=format_interpro_p2c("../example_data/interpro_ParentChildTreeFile.txt")

# inspect interpro tables
head(interpAnno)
```

    ##      gene        id                                           description
    ## 1 FAM162A IPR009432                   Protein of unknown function DUF1075
    ## 2 TMEM177 IPR026620                             Transmembrane protein 177
    ## 3    OSTC IPR021149 Oligosaccharyl transferase complex, subunit OST3/OST6
    ## 4    OSTC IPR042416        Oligosaccharyltransferase complex subunit OSTC
    ## 5  VPS37C IPR009851                         Modifier of rudimentary, Modr
    ## 6  VPS37C IPR037202                                 ESCRT assembly domain

``` r
head(interp.p2c)
```

    ##      parent     child
    ## 1 IPR000008 IPR014020
    ## 2 IPR000008 IPR033884
    ## 3 IPR000008 IPR037300
    ## 4 IPR000008 IPR037301
    ## 5 IPR000008 IPR037302
    ## 6 IPR000008 IPR037303

``` r
# creat anno object
interpro = createAnno(db='Interpro', gene.anno = interpAnno, gene.col='gene', id.col='id', desp.col = 'description',
                      genes = union(noto.genes, pcp.genes), bg.genes = bg.genes, 
                      hierarchy.df=interp.p2c, parent.col=1, child.col=2)
```

    ## 40 input genes are not found in the gene.anno table

    ## 290 input background genes are not found in the gene.anno table

# 3 Calculate functional similarities between genes

We will calculate the similarities between the functional terms that are
associated with different genes and use that as the functional
similarities between genes. To calculate the similarities between
functional terms in a hierarchical system, we will use methods for
calculating semantic similarity as described in [Yu et
al. 2010](doi:10.1093/bioinformatics/btq064) and their R package
[GOSemSim](https://bioconductor.org/packages/release/bioc/html/GOSemSim.html).
This involves 2 steps: 1. Calculate the information content (IC) of each
annotation term 2. Use the IC and term hierarchical relation to
calculate semantic similarities between genes’ annotated terms Functions
for computing IC and Semantic similarities are adapted from the
[GOSemSim package](https://github.com/YuLab-SMU/GOSemSim) package, which
were designed for working with only GO terms.

## 3.1 Calculate infomation content (IC) of each annotation term

This step requires both the parent-child relationship between terms and
the overall frequency of each term. The frequencies of terms will be
calculated from the `gene.anno` dataframes in the annotation objects,
which should contain not only the enriched genes, but also background
genes (e.g. all genes in the genome, or all genes expressed in embryo).
In our example, we used the union of genes enriched in each cell type.

``` r
## Reactome
reactome <- computeIC.anno(reactome)

## GO
### IC for terms in Molecular Funciton, Cellular Component, and Biological Function will be calculated separately
for(cat in c("BP","CC","MF")){
  go <- computeIC.anno(go, cat = cat) # the function will use GO relations from library(AnnotationDbi)
}

# Interpro
interpro <- computeIC.anno(interpro)
```

## 3.2 Calculate pairwise gene similarities based on individual annotation databases

Similarities are calculated for both genes enriched in cell types of
interest and background genes (genes enriched in other cell types).
Similarity between background genes will be used to calculate the priors
when combining annotation similarities from different databases. Since
this step is slow, a subset of background genes can be used to speed up
calculation. If only one database will be used, background genes can be
omitted.

### 3.2.1 Examples of pair-wise functional similarities

We can now calculate functional similarities between pairs of genes
based on each functional annotation database. The function
`geneSim.anno()` returns both the similarity score (between 0 and 1), as
well as the annotation terms associated with each gene.

``` r
## intepro 
geneSim.anno("SOX32","XBP1",interpro, measure="Jiang", cat=NULL, combine = "BMA")
```

    ## $geneSim
    ## [1] 0
    ## 
    ## $Gene1Annotations
    ## [1] "IPR009071" "IPR036910"
    ## 
    ## $Gene2Annotations
    ## [1] "IPR004827"

``` r
geneSim.anno("ATF6","XBP1",interpro, measure="Jiang", cat=NULL, combine = "BMA")
```

    ## $geneSim
    ## [1] 0.6667
    ## 
    ## $Gene1Annotations
    ## [1] "IPR004827" "IPR029801"
    ## 
    ## $Gene2Annotations
    ## [1] "IPR004827"

``` r
## go
cat="BP"
geneSim.anno("SOX32","XBP1",go, measure="Jiang", cat=cat,combine = "BMA")
```

    ## $geneSim
    ## [1] 0.3808
    ## 
    ## $Gene1Annotations
    ##  [1] "GO:0007507" "GO:0001570" "GO:0001525" "GO:0006355" "GO:0030154"
    ##  [6] "GO:0009653" "GO:0045944" "GO:0001706" "GO:0042663" "GO:0001946"
    ## [11] "GO:0061035" "GO:0003262" "GO:0007492" "GO:0045893" "GO:0042074"
    ## [16] "GO:0035701" "GO:0007493" "GO:0010468" "GO:0001889" "GO:0007498"
    ## [21] "GO:0043534" "GO:0003007" "GO:0007368" "GO:0035050" "GO:0001885"
    ## [26] "GO:0060216"
    ## 
    ## $Gene2Annotations
    ## [1] "GO:0006355" "GO:0006990" "GO:0035188" "GO:0006366" "GO:0001889"
    ## [6] "GO:0048785"

``` r
geneSim.anno("ATF6","XBP1",go, measure="Jiang", cat=cat,combine = "BMA")
```

    ## $geneSim
    ## [1] 0.59
    ## 
    ## $Gene1Annotations
    ## [1] "GO:0006355" "GO:0045944" "GO:0030968" "GO:0006357" "GO:0001666"
    ## 
    ## $Gene2Annotations
    ## [1] "GO:0006355" "GO:0006990" "GO:0035188" "GO:0006366" "GO:0001889"
    ## [6] "GO:0048785"

``` r
## reactome
geneSim.anno("SOX32","XBP1",reactome, measure="Jiang", cat=NULL, combine = "BMA")
```

    ## $geneSim
    ## [1] 0
    ## 
    ## $Gene1Annotations
    ## [1] "R-HSA-3769402"
    ## 
    ## $Gene2Annotations
    ## [1] "R-HSA-381038" "R-HSA-381070" "R-HSA-381183"

``` r
geneSim.anno("ATF6","XBP1",reactome, measure="Jiang", cat=NULL, combine = "BMA")
```

    ## $geneSim
    ## [1] 0.4888
    ## 
    ## $Gene1Annotations
    ## [1] "R-HSA-381033"
    ## 
    ## $Gene2Annotations
    ## [1] "R-HSA-381038" "R-HSA-381070" "R-HSA-381183"

### 3.2.2 Calculate similarities for all pairs of genes and store the results as matrices in the Anno objects

``` r
## Reactome 
reactome=func_sim(reactome, genes=reactome@genes)

## GO (calculate for each category separately)
for(cat in c("BP","MF","CC")){
  go=func_sim(go, cat = cat, genes=go@genes)
}

## Interpro
interpro=func_sim(interpro)
```

## 3.3 Save the objects with functional similarities and annotation tables used

``` r
saveRDS(reactome, "../example_results/reactome.rds")
saveRDS(go, "../example_results/go.rds")
saveRDS(interpro, "../example_results/interpro.rds")
```

## 3.4 Extract functional similarities from STRING database

## 3.5 read in the STRING protein-protein interaction (PPI) info file for zebrafish

STRING PPI tables can be downloaded from
<https://string-db.org/cgi/download.pl> (we used
7955.protein.links.full.v11.5.txt.gz). The downloaded table contains
Ensembl protein ids. We converted them into gene symbols that match what
we use in our scRNA-seq data and subsetted the table to contain only the
notochord/hatching gland enriched genes and the background genes.

``` r
str_zf_symb=read.csv("../example_data/string_table_GeneSymbol_sep_2021.csv",stringsAsFactors = F)
## print header
head(str_zf_symb)
```

    ##    gene1  gene2 neighborhood neighborhood_transferred fusion cooccurence
    ## 1 CCDC80 COL4A5            0                        0      0           0
    ## 2 CCDC80 FSTL1A            0                        0      0           0
    ## 3 CCDC80  ECSIT            0                        0      0           0
    ## 4 CCDC80 COL5A1            0                        0      0           0
    ## 5 CCDC80  ACTN1            0                        0      0           0
    ## 6 CCDC80   PAX9            0                        0      0           0
    ##   homology coexpression coexpression_transferred experiments
    ## 1        0            0                       63           0
    ## 2        0            0                      177           0
    ## 3        0            0                        0           0
    ## 4        0           53                      195           0
    ## 5        0            0                        0           0
    ## 6        0            0                        0           0
    ##   experiments_transferred database database_transferred textmining
    ## 1                       0        0                    0          0
    ## 2                       0        0                    0          0
    ## 3                       0        0                    0          0
    ## 4                       0        0                    0          0
    ## 5                       0        0                    0          0
    ## 6                       0        0                    0        434
    ##   textmining_transferred combined_score
    ## 1                    154            173
    ## 2                      0            177
    ## 3                    162            161
    ## 4                    364            472
    ## 5                    157            156
    ## 6                      0            434

*The column `$combined_score` is the default combined score from the
STRING database.*

## 3.6 Compute STRING PPI scores with selected evidence channels

The “database” scores in STRING are derived from a number of functional
annotation databases including GO and Reactome. We found that the STRING
“database” scores are more sparse compared to the annotation
similarities we calculated above. We will omit the database channel from
STRING and supplement functional similarity scores with the semantic
similarities calculated earlier. Expression channels will be omitted too
since we’ll use expression similarity calculated from our
transcriptional trajectories. The selected channels will be integrated
into a single score according to STRING’s evidence integration method
explained here:
<https://string-db.org/help//faq/#how-are-the-scores-computed>

``` r
evidence_use=c("experiments", "experiments_transferred","textmining","textmining_transferred", "homology")
genes_use=union(noto.genes, pcp.genes) 
str_ExpTxt = str_sim_M(string_tbl=str_zf_symb, evi=evidence_use, genes_use=genes_use, protein_cols=c("gene1","gene2"))

## Construct an additional string similarity matrix with the default combined scores
str_default = str_sim_M(string_tbl=str_zf_symb, evi="default", genes_use=genes_use, protein_cols=c("gene1","gene2"))
```

*`str_ExpTxt` and `str_default` are square gene x gene metrices with
similarity values (i.e. interaction scores) between 0 and 1*

## 3.7 Save STRING similarity matrices

``` r
str_scores=list("STRING"=str_default, "STRING_ExpTxt"=str_ExpTxt)
saveRDS(str_scores, "../example_results/STRING_score_matrices.rds")
```

# 4 Combinations functional similarity scores calculated from different databases

## 4.1 Read in the functional similarity scores calculated from individual databases earlier in this notebook

``` r
str_scores = readRDS("../example_results/STRING_score_matrices.rds")
go = readRDS("../example_results/go.rds")
reactome = readRDS("../example_results/reactome.rds")
interpro = readRDS("../example_results/interpro.rds")
```

## 4.2 Organize functional similarities into a 3D array to facilitate downstream processing

To make scores calculated from different databases more comparable, we
will adjust the scores from GO, Reactome, and Interpro with a “prior”,
which is the expected similarity for any random pair of genes. We
included the union of genes highly expressed in each cell type in the
embryo during our similarity calculations (`genes` and `bg.genes`). This
provides a global gene set for estimating similarities between random
gene pairs. We will calculate the mean similarity across all genes and
subtract this value from the similarity matrix. The STRING dataset is
already prior adjusted so don’t need to go through this procedure.

### 4.2.1 Construct the 3d functional similarity array for notochord enriched genes

``` r
anno.sim.noto=stack_func_sim(str_scores, adj_prior=F, genes_use=noto.genes, verbose = F) # string scores were already prior adjusted

anno.sim.noto=stack_func_sim(list("GO_cc"=go@cat.similarity$CC, "GO_bp"=go@cat.similarity$BP, "GO_mf"=go@cat.similarity$MF, 
                                  "Reactome"=reactome@similarity, "Interpro"=interpro@similarity), 
                             genes_use=noto.genes, adj_prior=T, add_to=anno.sim.noto, verbose = T)
```

    ## [1] "Adding GO_cc"
    ## [1] "Mean (prior) value to subtract = 0.443171335266195"
    ## [1] "Adding GO_bp"
    ## [1] "Mean (prior) value to subtract = 0.159372105837839"
    ## [1] "Adding GO_mf"
    ## [1] "Mean (prior) value to subtract = 0.257580234529552"
    ## [1] "Adding Reactome"
    ## [1] "Mean (prior) value to subtract = 0.0182937903553978"
    ## [1] "Adding Interpro"
    ## [1] "Mean (prior) value to subtract = 0.0126200897377966"
    ## [1] "Result array dimensions: 806 x 806 x 5"
    ## [1] "Sources used: GO_cc, GO_bp, GO_mf, Reactome, Interpro"
    ## [1] "Value range: 0 0.987379910262203"

### 4.2.2 Construct functional similarity array for hatching gland enriched genes

Set `genes_use = pcp.genes` instead of `noto.genes`.

``` r
anno.sim.pcp=stack_func_sim(str_scores, adj_prior=F, genes_use=pcp.genes, verbose = F) # string scores were already prior adjusted

anno.sim.pcp=stack_func_sim(list("GO_cc"=go@cat.similarity$CC, "GO_bp"=go@cat.similarity$BP, "GO_mf"=go@cat.similarity$MF, 
                                  "Reactome"=reactome@similarity, "Interpro"=interpro@similarity), 
                             genes_use=pcp.genes, adj_prior=T, add_to=anno.sim.pcp, verbose = F)
```

## 4.3 Combine similarities derived from different databases

The selected similarity measures will be integrated into a single score
by the evidence integration method used by STRING database:
<https://string-db.org/help//faq/#how-are-the-scores-computed>

``` r
## Sepecify which sources to combine in the format of a list, names are combined score name, entries are vectors of the source names to combine (subsets of the 3rd dimention names in the anno.sim arrays)
comb = list("STRING+GO+Reactome+Interpro" = c("STRING", "GO_cc", "GO_bp", "GO_mf", "Reactome", "Interpro"),
            "STRING_ExpTxt+GO_noCC+Reactome+Interpro" = c("STRING_ExpTxt", "GO_bp", "GO_mf", "Reactome", "Interpro"),
            "GO+Reactome+Interpro" = c("GO_cc", "GO_bp", "GO_mf", "Reactome", "Interpro"),
            "GO" = c("GO_cc", "GO_bp", "GO_mf"))

anno.sim.noto = add_combined_scores(anno.sim=anno.sim.noto, how.to=comb, add=T, verbose = T)
```

    ## [1] "Result array dimensions: 806 x 806 x 11"
    ## [1] "Sources used: STRING, STRING_ExpTxt, GO_cc, GO_bp, GO_mf, Reactome, Interpro, STRING+GO+Reactome+Interpro, STRING_ExpTxt+GO_noCC+Reactome+Interpro, GO+Reactome+Interpro, GO"
    ## [1] "Value range: 0 0.999999776549193"

``` r
anno.sim.pcp = add_combined_scores(anno.sim=anno.sim.pcp, how.to=comb, add=T, verbose = F)
```

## 4.4 Save the combined functional similarities into the MIMIR objects created in the previous step

Each object already contained the expression similarities
(`obj@exp.sim`) and the *enriched genes x pseudotime* expression matrix
(`obj@exp.data`).

``` r
noto.obj <- readRDS("../example_results/noto_mimir.rds")
noto.obj <- add_to_MIMIR(noto.obj,anno.sim.noto)
noto.obj <- add_anno_to_MIMIR(noto.obj, reactome)
noto.obj <- add_anno_to_MIMIR(noto.obj, go)
noto.obj <- add_anno_to_MIMIR(noto.obj, interpro)
## save updated MIMIR object
saveRDS(noto.obj, "../example_results/noto_mimir.rds")

## repeat for the hatching gland object
pcp.obj <- readRDS("../example_results/hg_mimir.rds")
pcp.obj <- add_to_MIMIR(pcp.obj,anno.sim.pcp)
pcp.obj <- add_anno_to_MIMIR(pcp.obj, reactome)
pcp.obj <- add_anno_to_MIMIR(pcp.obj, go)
pcp.obj <- add_anno_to_MIMIR(pcp.obj, interpro)
## save updated MIMIR object
saveRDS(pcp.obj, "../example_results/hg_mimir.rds")
```
