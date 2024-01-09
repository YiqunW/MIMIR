Step1: Define Enriched Genes in Cell Types of Interest
================
Yiqun Wang

In this example, we will identify genes that are enriched in the
zebrafish notochord and hatching gland during specification and
differentiation. Enriched genes are defined as the ones that show
stronger expression in the notochord and hatching gland compared to
other cell types. Two single-cell RNA-seq datasets will be used: 1.
whole embryo data from 3-12hpf, and 2. anterior embryo data from
14-24hpf.

## Define enriched genes from the whole embryo dataset

This dataset is from [Farrell et
al. 2018](https://pubmed.ncbi.nlm.nih.gov/29700225/). Single-cell
transcriptomes have been assembled into developmental trajectories using
the [URD package](https://github.com/farrellja/URD).

### Load data

``` r
suppressPackageStartupMessages(library(URD))
object <- readRDS("../../MIMIR large file/object_6_tree.rds")
plotTree(object,label.segments = T)
```

![](step1_Identify_Enriched_Genes_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Segment the developmental trajectories into several pseudotime windows

``` r
## Define functions to segment trajectories into equal pseudotime windows and to combine small adjacent segments
seg_tree <- function(object, num.windows=10, seg.rm=c("81","82","83")){
  pt.seg.tbl=object@tree$segment.pseudotime.limits
  pt.seg.tbl=pt.seg.tbl[setdiff(rownames(pt.seg.tbl),seg.rm),]
  ref.windows=seq(min(pt.seg.tbl[,"start"],na.rm=T),max(pt.seg.tbl[,"end"],na.rm = T),length.out = num.windows+1)
  node.tbl=c()
  for(seg in rownames(pt.seg.tbl)){
    seg.cells <- cellsInCluster(object, "segment", seg)
    seg.cells.pt=object@pseudotime[seg.cells,"pseudotime"]
    #seg.lim=as.numeric(pt.seg.tbl[seg,])
    seg.lim=c(min(seg.cells.pt,na.rm=T),max(seg.cells.pt,na.rm=T))
    ref.pt.ind=intersect(which(ref.windows>=seg.lim[1]),which(ref.windows<=seg.lim[2]))
    if(length(ref.pt.ind)==0){
      seg.names.tbl=c(rownames(node.tbl),seg)
      node.tbl=rbind(node.tbl,c(seg.lim,length(seg.cells)))
      rownames(node.tbl)=seg.names.tbl
    }else{
      add.pt=0
      if(min(ref.pt.ind)>1){
        if(sum(ref.windows[(min(ref.pt.ind)-1):min(ref.pt.ind)]) > seg.lim[1]*2){
          add.pt=add.pt+1
        }
      }
      if(max(ref.pt.ind)<length(ref.windows)){
        if(sum(ref.windows[max(ref.pt.ind):(max(ref.pt.ind)+1)])< seg.lim[2]*2){
          add.pt=add.pt+1
        }
      }
      num.pt=length(ref.pt.ind)+add.pt
      if(num.pt<=2){
        seg.names.tbl=c(rownames(node.tbl),seg)
        node.tbl=rbind(node.tbl,c(seg.lim,length(seg.cells)))
        rownames(node.tbl)=seg.names.tbl
      }else{
        seg.pt=seq(seg.lim[1],seg.lim[2],length.out = num.pt)
        for(i in 1:(length(seg.pt)-1)){
          st=seg.pt[i]
          ed=seg.pt[i+1]
          cells.in=intersect(which(seg.cells.pt>=st),which(seg.cells.pt<=ed))
          seg.names.tbl=c(rownames(node.tbl),paste0(seg,"_",i))
          node.tbl=rbind(node.tbl,c(st,ed,length(cells.in)))
          rownames(node.tbl)=seg.names.tbl
        }
      }
    }
  }
  colnames(node.tbl)=c("pt.start","pt.end","num.cells")
  return(node.tbl)
}
combine_ends <- function(seg_tbl,min.cells=40){
  combined_tbl=seg_tbl
  segs=unlist(lapply(rownames(seg_tbl), function(x) unlist(strsplit(x,"_"))[1]))
  for(seg in unique(segs)){
    sub_tbl=seg_tbl[which(segs==seg),,drop=F]
    if(dim(sub_tbl)[1]>1){
      num.cells=as.numeric(sub_tbl[,"num.cells"])
      below.ind=num.cells<min.cells
      sum.to.up=0
      if(below.ind[1]){
        sum.from.top=cumsum(num.cells)
        sum.to=min(which(sum.from.top>=min.cells))
        combined_tbl[rownames(sub_tbl)[1],]=c(sub_tbl[1,"pt.start"],sub_tbl[sum.to,"pt.end"],sum(sub_tbl[1:sum.to,"num.cells"]))
        inds=unlist(lapply(rownames(sub_tbl)[2:sum.to], function(x) unlist(strsplit(x,"_"))[2]))
        seg.name=paste0(rownames(sub_tbl)[1],"_",paste(inds,collapse = "_"))
        rownames(combined_tbl)[which(rownames(combined_tbl)==rownames(sub_tbl)[1])]=seg.name
        combined_tbl=combined_tbl[-which(rownames(combined_tbl)%in%rownames(sub_tbl)[1:sum.to]),]
        sum.to.up=sum.to
      }
      if(rev(below.ind)[1]){
        num.rows=dim(sub_tbl)[1]
        sum.from.bot=rev(cumsum(rev(num.cells)))
        sum.to=max(which(sum.from.bot>=min.cells))
        if(sum.to>sum.to.up){
          combined_tbl[rownames(sub_tbl)[sum.to],]=c(sub_tbl[sum.to,"pt.start"],sub_tbl[num.rows,"pt.end"],sum(sub_tbl[sum.to:num.rows,"num.cells"]))
          inds=unlist(lapply(rownames(sub_tbl)[(sum.to+1):num.rows], function(x) unlist(strsplit(x,"_"))[2]))
          seg.name=paste0(rownames(sub_tbl)[sum.to],"_",paste(inds,collapse = "_"))
          rownames(combined_tbl)[which(rownames(combined_tbl)==rownames(sub_tbl)[sum.to])]=seg.name
          combined_tbl=combined_tbl[-which(rownames(combined_tbl)%in%rownames(sub_tbl)[sum.to:num.rows]),]
        }
      }
    }
  }
  return(combined_tbl)
}
combine_segs <- function(seg_tbl,segs2combine){
  combined_tbl=seg_tbl
  for(seg.list in segs2combine){
    inds=unlist(lapply(seg.list[2:length(seg.list)], function(x) unlist(strsplit(x,"_"))[2:length(unlist(strsplit(x,"_")))]))
    seg.name=paste0(seg.list[1],"_",paste(inds,collapse = "_"))
    sub_tbl=seg_tbl[seg.list,]
    combined_tbl[seg.list[1],]=c(min(sub_tbl[,"pt.start"]),max(sub_tbl[,"pt.end"]),sum(sub_tbl[,"num.cells"]))
    rownames(combined_tbl)[which(rownames(combined_tbl)==seg.list[1])]=seg.name
    combined_tbl=combined_tbl[-which(rownames(combined_tbl)%in%seg.list),]
  }
  return(combined_tbl)
}
```

``` r
## Segment the tree
seg_tbl=seg_tree(object,num.windows=15,seg.rm=c("81","82","83")) # segments 81, 82, and 83 represent early, unspecified cells. They are removed from this analysis 

## fuse nodes that contain fewer than min.cells cells to adjacent segments
seg_tbl.com1=combine_ends(seg_tbl,min.cells = 60)

## manually combine small segments not fused in last step
segs2combine=list(c("17_3","17_4_5_6"),c("38_2","38_3","38_4","38_5","38_6","38_7"),
                  c("40_5","40_6","40_7","40_8"),c("53_1_2","53_3"))

seg_tbl.combined=combine_segs(seg_tbl.com1,segs2combine)

## Add new segment information to object@group.ids
new.id=rownames(object@group.ids)
names(new.id)=rownames(object@group.ids)
new.id[1:length(new.id)]="not.assigned"
for(seg in rownames(seg_tbl.combined)){
  rand_int=sample(10:99,1)
  group.name=paste0(rand_int,"*",seg)
  ori.seg=unlist(strsplit(seg,"_"))[1]
  seg.cells <- cellsInCluster(object, "segment", ori.seg)
  seg.cells.pt=object@pseudotime[seg.cells,"pseudotime"]
  group.cells=intersect(seg.cells[which(seg.cells.pt>seg_tbl.combined[seg,"pt.start"])],seg.cells[which(seg.cells.pt<=seg_tbl.combined[seg,"pt.end"])])
  new.id[group.cells]=group.name
}
object@group.ids["new.nodes"]=new.id[rownames(object@group.ids)]

object@tree$segment.names <- object@tree$segment.names.short

# Define axial mesoderm cells and the background (non-axial) cells 
axial.cells <- cellsInCluster(object, "segment", c("79", "29", "32"))
nonaxial.cells <- setdiff(unlist(object@tree$cells.in.segment), axial.cells)
noto.cells <- cellsInCluster(object, "segment", c("79", "32"))
pcp.cells <- cellsInCluster(object, "segment", c("79", "29"))

# Grab new.nodes along a lineage and also pseudotime match across other populations of cell
new.nodes.use <- function(object, cells.in.lineage, cells.pseudotime.match, pseudotime="pseudotime") {
  nnu <- setdiff(unique(object@group.ids[cells.in.lineage, "new.nodes"]), "not.assigned")
  nnu.nr <- unlist(lapply(strsplit(nnu, "\\*"), function(x) x[2]))
  nnu <- nnu[order(seg_tbl.combined[nnu.nr,"pt.start"])]
  nnu.nr <- nnu.nr[order(seg_tbl.combined[nnu.nr,"pt.start"])]
  cells.in.nodes <- lapply(nnu, function(n) rownames(object@group.ids)[which(object@group.ids$new.nodes == n)])
  cells.pt.matched <- lapply(nnu.nr, function(n) {
    pt.match <- object@pseudotime[cells.pseudotime.match,pseudotime]
    return(cells.pseudotime.match[which(pt.match >= seg_tbl.combined[n,"pt.start"] & pt.match <= seg_tbl.combined[n,"pt.end"])])
  })
  return(list(
    group.ids=nnu,
    seg.ids=nnu.nr,
    seg.info=seg_tbl.combined[nnu.nr,],
    cells.in.nodes=cells.in.nodes,
    cells.pt.matched=cells.pt.matched
  ))
}

noto.nodes <- new.nodes.use(object, noto.cells, nonaxial.cells)
pcp.nodes <- new.nodes.use(object, pcp.cells, nonaxial.cells)

## Visualize the segmentation on the tree
object@meta$noto.seg=0
for(i in noto.nodes$seg.info[,'pt.start']){
  cells=rownames(object@pseudotime)[which(object@pseudotime$pseudotime>i)]
  object@meta[cells,"noto.seg"]=object@meta[cells,"noto.seg"]+1
}
plotTree(object,"noto.seg",cell.alpha = 0.9)
```

![](step1_Identify_Enriched_Genes_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
object@meta$pcp.seg=0
for(i in pcp.nodes$seg.info[,'pt.start']){
  cells=rownames(object@pseudotime)[which(object@pseudotime$pseudotime>i)]
  object@meta[cells,"pcp.seg"]=object@meta[cells,"pcp.seg"]+1
}
plotTree(object,"pcp.seg",cell.alpha = 0.9)
```

![](step1_Identify_Enriched_Genes_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

### Run precision-recall (AUCPR) to find markers of notochord and hatching gland

``` r
## AUCPR To Find Markers -------------
# Compare each axial mesoderm window to its non-axial pseudotime-matched cells
noto.markers <- lapply(1:length(noto.nodes$group.ids), function(i) {
  message(paste0(Sys.time(), ": ", i))
  markers <- markersAUCPR(object, cells.1=noto.nodes$cells.in.nodes[[i]], cells.2=noto.nodes$cells.pt.matched[[i]], effect.size = 0.25, frac.must.express=0.1)
  thresh <- aucprThreshold(cells.1=noto.nodes$cells.in.nodes[[i]], cells.2=noto.nodes$cells.pt.matched[[i]], factor = 1, max.auc = Inf)
  names(markers)[4:7] <- c("posFrac_noto", "posFrac_rest", "nTrans_noto", "nTrans_rest")
  return(markers)
})

pcp.markers <- lapply(1:length(pcp.nodes$group.ids), function(i) {
  message(paste0(Sys.time(), ": ", i))
  markers <- markersAUCPR(object, cells.1=pcp.nodes$cells.in.nodes[[i]], cells.2=pcp.nodes$cells.pt.matched[[i]], effect.size = 0.25, frac.must.express=0.1)
  thresh <- aucprThreshold(cells.1=pcp.nodes$cells.in.nodes[[i]], cells.2=pcp.nodes$cells.pt.matched[[i]], factor = 1, max.auc = Inf)
  names(markers)[4:7] <- c("posFrac_pcp", "posFrac_rest", "nTrans_pcp", "nTrans_rest")
  return(markers)
})

saveRDS(list("noto"=noto.markers,"pcp"=pcp.markers), file='../example_results/markers3-12hpf.rds')
```

``` r
markers=readRDS('../example_results/markers3-12hpf.rds')
noto.markers=markers$noto
pcp.markers=markers$pcp
```

### Subset results to identify enriched genes with the following criteria:

1.  AUCPR is 2x random expectation
2.  log fold change is &gt;0.3
3.  Marker either passes in the last window or at least 2 windows.

``` r
noto.markers.good <- lapply(noto.markers, function(x) {
  x[x$AUCPR.ratio >= 2 & x$nTrans_noto-x$nTrans_rest > 0.3,]
})
noto.markers.2.fc3.2x <- names(which(table(unlist(lapply(noto.markers.good, rownames))) > 1))
noto.markers.2.fc3.last <- rownames(noto.markers.good[[length(noto.markers.good)]])
noto.markers.2.fc3.2xorlast <- unique(c(noto.markers.2.fc3.2x, noto.markers.2.fc3.last))
noto.genes <- noto.markers.2.fc3.2xorlast

pcp.markers.good <- lapply(pcp.markers, function(x) {
  x[x$AUCPR.ratio >= 2 & (x$nTrans_pcp-x$nTrans_rest) > 0.3,]
})
pcp.markers.2.fc3.2x <- names(which(table(unlist(lapply(pcp.markers.good, rownames))) > 1))
pcp.markers.2.fc3.last <- rownames(pcp.markers.good[[length(pcp.markers.good)]])
pcp.markers.2.fc3.2xorlast <- unique(c(pcp.markers.2.fc3.2x, pcp.markers.2.fc3.last))
pcp.genes <- pcp.markers.2.fc3.2xorlast
```

## Define enriched genes from the anterior embryo dataset

This dataset is from [Raj et
al. 2020](https://pubmed.ncbi.nlm.nih.gov/33068532/). Axial mesoderm
cells (notochord and hatching gland) have been extracted from this data
and combined with the previous dataset to reconstruct extended axial
mesoderm trajectories. \#\#\# Load data

``` r
bushra.full <- readRDS("../../MIMIR large file/BushraDataFullURD.rds") # sc-transcriptomes from 14-24hpf anterior embryos
axial.tree <- readRDS("../example_data/DraftTree.rds") # full trajectories for notochord and hatching gland using axial mesoderm data from both datasets
```

### Find markers for 14-24hpf notochord and hatching gland using AUCPR

``` r
# Make sure both objects have same stage assignments
axial.tree@group.ids$hpf <- axial.tree@meta$hpf
bushra.full@group.ids$hpf <- gsub("h", "", bushra.full@meta$stage)

stages.consider <- list(
  c("12", "14"),
  c("16", "18"),
  c("20", "24")
)

noto.markers <- lapply(stages.consider, function(stage) {
  message(paste0(Sys.time(), ": ", stage))
  branch <- "Notochord"
  branch.short <- "noto"
  cells.in <- intersect(intersect(cellsAlongLineage(axial.tree, branch), cellsInCluster(axial.tree, "hpf", stage)), colnames(bushra.full@logupx.data))
  cells.out <- setdiff(cellsInCluster(bushra.full, "hpf", stage), cellsAlongLineage(axial.tree, branch))
  markers <- markersAUCPR(bushra.full, cells.1=cells.in, cells.2=cells.out, effect.size = 0.25, frac.must.express=0.1)
  thresh <- aucprThreshold(cells.1=cells.in, cells.2=cells.out, factor = 1, max.auc = Inf)
  names(markers)[4:7] <- c(paste0("posFrac_", branch.short), "posFrac_rest", paste0("nTrans_", branch.short), "nTrans_rest")
  return(markers)
})

pcp.markers <- lapply(stages.consider, function(stage) {
  message(paste0(Sys.time(), ": ", stage))
  branch <- "Prechordal Plate"
  branch.short <- "pcp"
  cells.in <- intersect(intersect(cellsAlongLineage(axial.tree, branch), cellsInCluster(axial.tree, "hpf", stage)), colnames(bushra.full@logupx.data))
  cells.out <- setdiff(cellsInCluster(bushra.full, "hpf", stage), cellsAlongLineage(axial.tree, branch))
  markers <- markersAUCPR(bushra.full, cells.1=cells.in, cells.2=cells.out, effect.size = 0.25, frac.must.express=0.1)
  thresh <- aucprThreshold(cells.1=cells.in, cells.2=cells.out, factor = 1, max.auc = Inf)
  names(markers)[4:7] <- c(paste0("posFrac_", branch.short), "posFrac_rest", paste0("nTrans_", branch.short), "nTrans_rest")
  return(markers)
})

saveRDS(list("noto"=noto.markers,"pcp"=pcp.markers), file='../example_results/markers14-24hpf.rds')
```

``` r
markers=readRDS('../example_results/markers14-24hpf.rds')
noto.markers=markers$noto
pcp.markers=markers$pcp
```

### Subset results to identify enriched genes

*Criteria for subsetting the results are determined by altering
thresholds and visually checking the expression of newly added or lost
genes plotted on the developmental trajectories.*

``` r
noto.markers.good <- lapply(noto.markers, function (nm) {
  rownames(nm)[which(nm$AUCPR.ratio >= 7.5 & (nm$nTrans_noto-nm$nTrans_rest) >= 0.3)]
})
pcp.markers.good <- lapply(pcp.markers, function (pm) {
  rownames(pm)[which(pm$AUCPR.ratio >= 7.5 & (pm$nTrans_pcp-pm$nTrans_rest) >= 0.3)]
})

noto.markers.good.2 <- names(which(table(unlist(noto.markers.good)) > 1))

pcp.markers.good.1 <- names(which(table(unlist(pcp.markers.good)) == 1))
pcp.markers.good.2 <- names(which(table(unlist(pcp.markers.good)) > 1))
pcp.markers.good.1.new <- setdiff(pcp.markers.good.1, pcp.genes)
pcp.markers.good.2.new <- setdiff(pcp.markers.good.2, pcp.genes)

pcp.markers.good.1.new.highAUCPR <- intersect(
  pcp.markers.good.1.new,
  unlist(lapply(1:3, function(i) {
    pm <- pcp.markers[[i]][pcp.markers.good[[i]],]
    rownames(pm)[which(pm$AUCPR >= 0.1)]
  }))
)
```

## Combine enriched genes identified from both datasets and save the result

``` r
pcp.genes <- sort(unique(c(pcp.genes, pcp.markers.good.2, pcp.markers.good.1.new.highAUCPR)))
noto.genes <- sort(unique(c(noto.genes, noto.markers.good.2)))
print(paste0("# notochord enriched genes: ",length(noto.genes)))
```

    ## [1] "# notochord enriched genes: 806"

``` r
print(paste0("# hatching gland enriched genes: ",length(pcp.genes)))
```

    ## [1] "# hatching gland enriched genes: 802"

``` r
write(pcp.genes, file="../example_results/genes.pcp.txt")
write(noto.genes, file="../example_results/genes.noto.txt")
```