require(URD)
## Define functions to segment trajectories into roughly equal pseudotime windows (branching points will be preserved as segmentation points)
seg_tree <- function(object, num.windows=10, seg.rm=NULL){
  ## cut the pseudotime in num.window chunks
  pt.seg.tbl=object@tree$segment.pseudotime.limits
  pt.seg.tbl=pt.seg.tbl[setdiff(rownames(pt.seg.tbl),seg.rm),]
  ref.windows=seq(min(pt.seg.tbl[,"start"],na.rm=T),max(pt.seg.tbl[,"end"],na.rm = T),length.out = num.windows+1)
  ## ref.windows is a vector containing num.window+1 values that evenly divide the pseudotime line into num.windows chunks
  node.tbl=c()
  ## Each tree segment between 2 branching points will be processed. 
  ## If a segment is small and doesn't contain a pseudotime division point (as in ref.windows), the segment will be preserved in the new segmentation.
  ## If a segment is large and contain at least 1 pseudotime division points, the segment will be divided into 2 or more qual chunks 
  for(seg in rownames(pt.seg.tbl)){
    seg.cells <- cellsInCluster(object, "segment", seg)
    seg.cells.pt=object@pseudotime[seg.cells,"pseudotime"]
    #seg.lim=as.numeric(pt.seg.tbl[seg,])
    seg.lim=c(min(seg.cells.pt,na.rm=T),max(seg.cells.pt,na.rm=T)) #min and max pseudotime in the tree segment
    ref.pt.ind=intersect(which(ref.windows>=seg.lim[1]),which(ref.windows<=seg.lim[2])) #indices of the cutting points in the segment
    if(length(ref.pt.ind)==0){
      ## if the segment doesn't contain a pseudotime cutting point, then include it as a new segment in the segmentation table
      seg.names.tbl=c(rownames(node.tbl),seg)
      node.tbl=rbind(node.tbl,c(seg.lim,length(seg.cells)))
      rownames(node.tbl)=seg.names.tbl
    }else{
      ## If the segment contains at least 1 pseudotime cutting point, cut the segment into equal pseudotime chunks 
      ## If n pseudotime cutting points are in the segment, the default is to divide the segment into n-1 equal chunks
      ## Add 1 chunk if the first cutting point is a lot later than the start of the segment
      ## Similarly add 1 chunk if the last cutting point is a lot earlier than the end of the segment
      add.pt=0
      if(min(ref.pt.ind)>1){
        # if the chunk in the segment before the first cutting point is larger than half of a standard pseudotime window, need 1 more new segment
        if(mean(ref.windows[(min(ref.pt.ind)-1):min(ref.pt.ind)]) > seg.lim[1]){
          add.pt=add.pt+1
        }
      }
      if(max(ref.pt.ind)<length(ref.windows)){
        # if the chunk after the last cutting point is larger than half of a standard pseudotime window, need 1 more new segment
        if(mean(ref.windows[max(ref.pt.ind):(max(ref.pt.ind)+1)])< seg.lim[2]){
          add.pt=add.pt+1
        }
      }
      num.pt=length(ref.pt.ind)+add.pt # divide the segment into num.pt-1 of new chunks
      if(num.pt<=2){
        ## keep the entire segment as a whole
        seg.names.tbl=c(rownames(node.tbl),seg)
        node.tbl=rbind(node.tbl,c(seg.lim,length(seg.cells)))
        rownames(node.tbl)=seg.names.tbl
      }else{
        ## cut the segment into num.pt-1 equal pseudotime chunks
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
  ## return a data table with the start pseudotime, end pseudotime, and number of cells for each new segment
  ## rownames preserves old segment names, but with _n added to the end for segments that are splitted.
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

## add the new segmentation as a new column in group.ids of the URD object
## seg_tbl is the output of the seg_tree() function
add_seg_id <- function(object, seg_tbl, new_group="seg_nodes"){
  new.id=rownames(object@group.ids)
  names(new.id)=rownames(object@group.ids)
  new.id[1:length(new.id)]="not.assigned"
  for(seg in rownames(seg_tbl)){
    rand_int=sample(10:99,1)
    group.name=paste0(rand_int,"*",seg)
    ori.seg=unlist(strsplit(seg,"_"))[1]
    seg.cells <- cellsInCluster(object, "segment", ori.seg)
    seg.cells.pt=object@pseudotime[seg.cells,"pseudotime"]
    group.cells=intersect(seg.cells[which(seg.cells.pt>seg_tbl[seg,"pt.start"])],seg.cells[which(seg.cells.pt<=seg_tbl[seg,"pt.end"])])
    new.id[group.cells]=group.name
  }
  object@group.ids[new_group]=new.id[rownames(object@group.ids)]
  return(object)
}

# Grab new.nodes along a lineage and also pseudotime match across other populations of cell
new.nodes.match <- function(object, cells.in.lineage, bg.cells, pseudotime="pseudotime", seg_use="seg_nodes", update_group="matched_seg", plt=T) {
  lineage.segs=object@group.ids[cells.in.lineage,seg_use,drop=F]
  segs=setdiff(unique(lineage.segs[,seg_use]), "not.assigned")
  object@group.ids[,update_group]="not.assigned"
  start_pt=c()
  for(s in segs){
    cells.in.seg=rownames(lineage.segs)[which(lineage.segs[,seg_use]==s)]
    pt_range=range(object@pseudotime[cells.in.seg,'pseudotime'], na.rm=T)
    ind.in.range=which(object@pseudotime[,'pseudotime']>=pt_range[1])
    ind.in.range=intersect(ind.in.range,which(object@pseudotime[,'pseudotime']<=pt_range[2]))
    object@group.ids[ind.in.range, update_group]=s
    start_pt=c(start_pt, pt_range[1])
  }
  start_pt=order(start_pt)
  start_pt_order=start_pt*0
  start_pt_order[start_pt]=c(1:length(start_pt))
  names(start_pt_order)=segs
  start_pt_order['not.assigned']="not.assigned"
  new_names=start_pt_order[object@group.ids[ ,update_group]]
  object@group.ids[, update_group]=new_names
  if(plt){
    library(RColorBrewer)
    print(plotTree(object, update_group, cell.alpha = 0.9, discrete.colors = brewer.pal(n=length(unique(object@group.ids[,update_group])),name = "Set3")))
  }
  return(object)
}

markersAUCPR_segs <- function(object, group_id, cells.in.lineage, bg.cells, effect.size=0.25, 
                              frac.must.express=0.1, thres_factor=1, max.auc=Inf){
  segs=setdiff(unique(object@group.ids[,group_id]),"not.assigned")
  nsegs=length(segs)
  lineage.markers <- lapply(1:length(segs), function(i) {
    message(paste0(Sys.time(), ": calculating markers for ", i, " out of ", nsegs, " time windows."))
    seg=segs[i]
    cells.in.seg=cellsInCluster(object, group_id, seg)
    cells.in=intersect(cells.in.lineage, cells.in.seg)
    cells.out=intersect(bg.cells, cells.in.seg)
    markers <- markersAUCPR(object, cells.1=cells.in, cells.2=cells.out, effect.size = effect.size, frac.must.express=frac.must.express)
    thresh <- aucprThreshold(cells.1=cells.in, cells.2=cells.out, factor = thres_factor, max.auc = max.auc)
    names(markers)[4:7] <- c("posFrac", "posFrac_rest", "nTrans", "nTrans_rest")
    return(markers)
  })
  names(lineage.markers)=segs
  return(lineage.markers)
}

