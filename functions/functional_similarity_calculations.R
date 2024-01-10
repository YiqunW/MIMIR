require(AnnotationDbi)
require(GOSemSim)
require(GO.db)
require(abind)

computeIC.anno <- function(goAnno, cat=NULL,Offsprings=NULL) {
  ## goAnno, dataframe with columns gene, id, evidence, category
  if(!is.null(cat)){
    goAnno=goAnno[goAnno$category == cat,]
  }
  goAnno_id=goAnno$id
  gocount=table(goAnno_id)
  gonames=names(gocount)
  gocount=as.vector(gocount)
  names(gocount)=gonames
  
  if(!is.null(cat)){
    Offsprings <- switch(cat,
                         MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                         BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                         CC = AnnotationDbi::as.list(GOCCOFFSPRING))
  }else{
    if(is.null(Offsprings)){
      print("Need to supply Offspring as a list")
      return(NA)
    }
  }
  cnt <- gocount[gonames] + sapply(gonames, function(i) sum(gocount[Offsprings[[i]]], na.rm=TRUE))
  names(cnt) <- gonames
  ## cnt is the number of descendents of each go_id (including the occurance of the go_id)
  ## gocount is the occurance of each go_id
  ## the probabilities of occurrence of GO terms in a specific corpus.
  p <- cnt/sum(gocount)
  ## IC of GO terms was quantified as the negative log likelihood.
  IC <- -log(p)
  return(IC)
}

offspring.list <- function(p2c){
  c.list=list()
  ps=unique(p2c[,1])
  cs=unique(p2c[,2])
  c.list[setdiff(cs,ps)]=NA
  for(p in unique(ps)){
    c.list[[p]]=get_children(p2c,p)
    if(length(c.list[[p]])==0){
      c.list[[p]]=NA
    }
  }
  return(c.list)
}

get_children <- function(p2c, p){
  all_c=c()
  M_cs=p2c[which(p2c[,1]==p),2]
  if(length(M_cs)>0){
    all_c=unique(c(all_c,M_cs))
    M_cs=intersect(M_cs,p2c[,1])
    if(length(M_cs)>0){
      for(M_c in M_cs){
        all_c=unique(c(all_c,get_children(p2c,M_c)))
      }
    }
  }
  return(all_c)
}

ancestor.list <- function(p2c, children){
  p.list=list()
  for(c in unique(children)){
    if(c%in%p2c[,2]){
      ps=get_parents(p2c,c)
      if(length(ps)>0){
        p.list[[c]]=ps
      }
    }
  }
  return(p.list)
}

get_parents <- function(p2c, c){
  all_p=c()
  M_ps=p2c[which(p2c[,2]==c),1]
  if(length(M_ps)>0){
    all_p=unique(c(all_p,M_ps))
    M_ps=intersect(M_ps,p2c[,2])
    if(length(M_ps)>0){
      for(M_p in M_ps){
        all_p=unique(c(all_p,get_parents(p2c,M_p)))
      }
    }
  }
  return(all_p)
}

geneSim.anno <- function(gene1, gene2, goAnno, IC, measure=c("Jiang","Wang"), combine="BMA",cat=NULL,anno.p2c=NULL,anc.all=NULL){
  ## IC should be named vector
  go1 <- goAnno[which(goAnno$gene==gene1),]
  go2 <- goAnno[which(goAnno$gene==gene2),]
  if(!is.null(cat)){
    go1 <- go1[which(go1$category==cat),]
    go2 <- go2[which(go2$category==cat),]
    goAnno=goAnno[which(goAnno$category==cat),]
  }
  go1=go1$id
  go2=go2$id
  go1=unique(go1[!is.na(go1)])
  go2=unique(go2[!is.na(go2)]) #unique shouldn't be neccessary
  if(length(go1) == 0 || length(go2) == 0){
    return (NA)
  }
  if(is.null(cat)){
    if(is.null(anno.p2c)){
      if(is.null(anc.all)){
        print("please supply parent to child table for non-GO annotations.")
        return(NA)
      }else{
        anc=anc.all[union(go1,go2)]
      }
    }else{
      if(is.null(anc.all)){
        anc=ancestor.list(anno.p2c,union(go1,go2))
      }else{
        anc=anc.all[union(go1,go2)]
      }
    }
  }
  res <- annoSim(go1, go2, goAnno=goAnno, IC=IC, measure=measure, combine=combine, cat=cat, anc=anc)
  return (list(geneSim=res, GO1=go1, GO2=go2))
}

annoSim <- function(GO1, GO2, goAnno, IC, measure=c("Jiang","Wang"), combine="BMA", cat=NULL, anc=NULL, precision=4){
  scores <- termSim.anno(GO1, GO2, goAnno, IC, method=measure, cat=cat, anc=anc)
  #print(scores)
  res <- combineScores.anno(scores, combine,digits=precision)
  return(round(res,digits = precision))
}

termSim.anno <- function(go1, go2, goAnno, IC, method=c("Wang","Resnik","Rel","Jiang","Lin"), cat=NULL, anc=NULL) {
  if(!method%in%c("Wang","Resnik","Rel","Jiang","Lin")){
    print("method has to be one of Wang, Resnik, Rel, Jiang, Lin")
    return(NA)
  }
  if(all(is.na(go1)) || all(is.na(go2))){
    return (NA)
  }
  go1 <- unique(go1[!is.na(go1)])
  go2 <- unique(go2[!is.na(go2)])
  
  if ( method %in% c("Resnik", "Jiang", "Lin", "Rel") ) {
    return(infoContentMethod.anno(go1, go2, goAnno, IC, method=method, cat=cat, anc=anc))
  } else if ( method == "Wang" ) {
    #return(wangMethod(go1, go2, cat=cat))
    print("function for graph based method not ready.")
    return(NA)
  }
}

getAncestors <- function(ont) {
  Ancestors <- switch(ont,
                      MF = "GOMFANCESTOR",
                      BP = "GOBPANCESTOR",
                      CC = "GOCCANCESTOR",
                      DO = "DO.db::DOANCESTOR"
  )
  if (ont == "DO") {
    db <- "DO.db"
    ## require(db, character.only=TRUE)
    requireNamespace(db)
  }
  return (eval(parse(text=Ancestors)))
}

infoContentMethod.anno <- function(ID1, ID2, goAnno, method="Jiang", IC=NULL, cat=NULL, anc=NULL) {
  ## IC is biased
  ## because the IC of a term is dependent on its children but not on its parents.
  if (length(IC) == 0) {
    stop("IC data not found, calculate with computeIC()")
  }
  if(!is.null(cat)){
    anc <- AnnotationDbi::as.list(getAncestors(cat))
    allid <- union(ID1, ID2)
    anc <- anc[intersect(allid,names(anc))]
    anc <- anc[!vapply(anc, function(x) is.null(x)||all(is.na(x)), logical(1))]
  }else if(is.null(anc)){
    print("Please supply ancestor lists.")
    return(NA)
  }else{
    allid <- union(ID1, ID2)
    for(id in allid){
      if(!"all"%in%anc[[id]]){
        anc[[id]]=c(anc[[id]],"all")
      }
      anc=anc[names(anc[which(!is.na(names(anc)))])]
    }
  }
  ## anc is a list of all parent nodes for the ids
  ## each entry is named with an id (input go term)
  return (GOSemSim:::infoContentMethod_cpp(ID1, ID2, anc, IC, method, ont="all")) #ont is used to set the root IC to 0. If not one of "BP,CC,MF. it wouldn't change the input IC.
}


# infoContentMethod_cpp <- function(id1_, id2_, anc_, ic_, method_, ont_) {
#   .Call('_GOSemSim_infoContentMethod_cpp', PACKAGE = 'GOSemSim', id1_, id2_, anc_, ic_, method_, ont_)
# }
# 

combineScores.anno <- function(SimScores, combine, digits=4) {
  ####!! modified from the original function to properly calculate cases with 1 row/column 
  if (length(combine) == 0) {  #if not define combine
    return(round(SimScores, digits=digits))
  }
  
  ## if combine was defined...
  if(!sum(!is.na(SimScores))) return (NA)
  
  if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) { ##if 1 row/column, 
    if (combine == "avg") {
      return(round(mean(SimScores, na.rm=TRUE), digits=digits))
    } else if(combine=="max" || combine=="rcmax"){
      return (round(max(SimScores, na.rm=TRUE), digits=digits))
    }else{ ## added by YW for combine=BMA
      return(round((sum(SimScores, na.rm=TRUE)+max(SimScores, na.rm=TRUE))/(1+length(SimScores)), digits=digits))
    }
  }
  
  ##remove all na rows and cols
  row.na.idx <- apply(SimScores, 1, function(i) all(is.na(i)))
  if (any(row.na.idx)) {
    SimScores <- SimScores[-which(row.na.idx), ]
  }
  
  if (! is.null(dim(SimScores)) ) {
    col.na.idx <- apply(SimScores, 2, function(i) all(is.na(i)))
    if (any(col.na.idx)) {
      SimScores <- SimScores[ , -which(col.na.idx)]
    }
  }
  if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
    if (combine == "avg") {
      return(round(mean(SimScores, na.rm=TRUE), digits=digits))
    } else if(combine=="max" || combine=="rcmax"){
      return (round(max(SimScores, na.rm=TRUE), digits=digits))
    } else{ ## added by YW
      return(round((sum(SimScores, na.rm=TRUE)+max(SimScores, na.rm=TRUE))/(1+length(SimScores)), digits=digits))
    }
  }
  
  ## if the remain is a matrix
  if (combine == "avg") {
    result   <- mean(SimScores, na.rm=TRUE)
  } else if (combine == "max") {
    result   <- max(SimScores, na.rm=TRUE)
  } else if (combine == "rcmax") {
    rowScore <- mean(apply(SimScores, 1, max, na.rm=TRUE))
    colScore <- mean(apply(SimScores, 2, max, na.rm=TRUE))
    result   <- max(rowScore, colScore)
  } else if (combine == "rcmax.avg" || combine == "BMA") {
    result   <- sum( apply(SimScores, 1, max, na.rm=TRUE),
                     apply(SimScores, 2, max, na.rm=TRUE)
    ) / sum(dim(SimScores))
  }
  
  return (round(result, digits=digits))
}

str_comb <- function(str_tbl_full, protein_cols=c("protein1","protein2"),
                     evi=c("neighborhood","fusion","cooccurence","coexpression","experimental","database","textmining")){
  ## input is string table not scaled by 1/1000
  comb_tbl=str_tbl_full[,c(protein_cols,"combined_score")]
  comb_tbl[,"combined_score"]=0
  if("homology"%in%evi){
    homology=as.numeric(str_tbl_full[,"homology"])/1000
    evi=setdiff(evi,"homology")
  }else{
    homology=as.numeric(str_tbl_full[,evi[1]]*0)
    if(length(intersect(evi, c("textmining","textmining_transferred","cooccurence")))>0){
      str_tbl_full[,intersect(evi, c("textmining","textmining_transferred","cooccurence"))]=0
    }
  }
  str_tbl_full=(str_tbl_full[,evi]-41) ## 41 is the prior used by STRING
  str_tbl_full=str_tbl_full*(str_tbl_full>0) ## scores smaller than prior are set to 0 (instead of a negative number)
  str_tbl_full=str_tbl_full/(1000-41) ## compute away prior completely and convert to probabilities (>0, <1)
  
  ## homology correction only required for textmining and cooccurence
  homo0_ind=which(homology==0)
  # if homology is 0, or all channels require homology adjustment is 0, no adjustment needs to be made
  homo0_ind=union(homo0_ind,
                  which(apply(str_tbl_full[,intersect(evi, c("textmining","textmining_transferred","cooccurence")),drop=F],1,sum)==0))
  str_tbl_minus = 1 - str_tbl_full
  comb_tbl[homo0_ind,"combined_score"]=1-apply(str_tbl_minus[homo0_ind,,drop=F],1,prod)*(1-0.041)
  # for the rest, homology adjustment is needed before computing combined scores
  for(i in setdiff(1:dim(str_tbl_full)[1],homo0_ind)){
    vec=str_tbl_minus[i,]
    piS_minus=1
    homo=homology[i]
    if(length(setdiff(evi,c("textmining","textmining_transferred","cooccurence")))>0){
      piS_minus=piS_minus*prod(vec[setdiff(evi, c("textmining","textmining_transferred","cooccurence"))])
    }
    if("textmining"%in%evi && "textmining_transferred"%in%evi){
      piS_minus=piS_minus*(homo+(1-homo)*prod(vec[c("textmining","textmining_transferred")]))
    }else{
      piS_minus=piS_minus*prod(vec[intersect(evi,c("textmining","textmining_transferred","cooccurence"))]*(1-homo)+homo)
    }
    comb_tbl[i,"combined_score"]=1 - piS_minus*(1-0.041)
  }
  return(comb_tbl)
}

## define a function that converts gene-gene interaction table from long to symmetric form
str_ppi_2tbl <- function(str_ppi, thres, values="combined_score",g1.col="protein1", g2.col="protein2",method="max"){
  all_g=union(str_ppi[[g1.col]],str_ppi[[g2.col]])
  str_ppi.tbl.max=matrix(0,nrow=length(all_g), ncol = length(all_g), dimnames = list(sort(all_g),sort(all_g)))
  unique_pairs=unique(str_ppi[,c(g1.col,g2.col)])
  for(i in 1:dim(unique_pairs)[1]){
    p1=as.character(unique_pairs[i,g1.col])
    p2=as.character(unique_pairs[i,g2.col])
    pair_ind=intersect(which(str_ppi[[g1.col]]==p1),which(str_ppi[[g2.col]]==p2))
    scores=str_ppi[pair_ind,values]
    if(length(scores)>0){
      scores=scores[which(scores>=thres)]
      if(length(scores)>0){
        if(method=="max"){
          str_ppi.tbl.max[p1,p2]=max(scores)
        }else if(method=="mean"){
          str_ppi.tbl.max[p1,p2]=mean(scores)
        }else if(method=="median"){
          str_ppi.tbl.max[p1,p2]=median(scores)
        }
        str_ppi.tbl.max[p2,p1]=str_ppi.tbl.max[p1,p2]
      }
    }
  }
  all.zeros=which(apply(str_ppi.tbl.max,1,sum)==0)
  if(length(all.zeros)>0){
    str_ppi.tbl.max=str_ppi.tbl.max[-all.zeros,-all.zeros]
  }
  return(str_ppi.tbl.max)
}

subtract_mean <- function(anno.matrix, prt_mean=F, prior=NULL){
  if(is.null(prior)){
    bg.ave=mean(anno.matrix[upper.tri(anno.matrix)])
  }else{
    bg.ave=prior
  }
  anno.matrix=anno.matrix-bg.ave
  anno.matrix=anno.matrix*(anno.matrix>0)
  if(prt_mean){print(paste0("Mean (prior) value to subtract = ", bg.ave))}
  return(anno.matrix)
}

stack_func_sim <- function(simM_list, genes_use=NULL, adj_prior=F, add_to=NULL){
  if(is.null(genes_use)&is.null(add_to)){
    stop("Need to provide either genes_use as a vector of genes or add_to as a similarity matrix")
  }
  if(is.null(genes_use)){
    genes_use=dimnames(add_to)[[1]]
  }else{
    if(!is.null(add_to)){
      if(setequal(genes_use, dimnames(add_to)[[1]])){
        genes_use=dimnames(add_to)[[1]]
      }else{
        if(length(intersect(genes_use, dimnames(add_to)[[1]]))>0){
          genes_use=dimnames(add_to)[[1]]
          print(paste0("Only ", length(intersect(genes_use, dimnames(add_to)[[1]])))," genes overalp between genes_use and genes in the add_to matrix. Ignoring gene_use and using genes in the add_to matrix.")
        }else{
          stop("There's no overlap between genes_use and genes in the add_to matrix. Consider setting one of them to NULL")
        }
      }
    }
  }
  anno_names=names(simM_list)
  anno3d=array(0,dim=c(length(genes_use),length(genes_use),length(anno_names)),
               dimnames = list(genes_use, genes_use, anno_names))
  for(anno in anno_names){
    simM=simM_list[[anno]]
    if(is.null(adj_prior) || adj_prior){
      if(is.numeric(adj_prior)){
        simM=subtract_mean(simM, prior=adj_prior)
      }else{
        simM=subtract_mean(simM, prior=NULL, prt_mean = T)
      }
    }
    comm.genes=intersect(genes_use,rownames(simM))
    anno3d[comm.genes, comm.genes, anno]=simM[comm.genes, comm.genes]
  }
  if(is.null(add_to)){
    return(anno3d)
  }else{
    return(abind(add_to, anno3d))
  }
}

combine_anno_scores <- function(all_links, channel.use, prior=NULL){
  if(length(channel.use)==1){
    if(is.null(prior)){
      return(all_links[,,channel.use])
    }else{
      adj.M=all_links[,,channel.use]-prior
      return(adj.M*(adj.M>0))
    }
  }else{
    if(is.null(prior)){
      return(apply(all_links[,,channel.use],c(1,2),function(x) 1-prod(1-x)))
    }else{
      adj.M=all_links[,,channel.use]-prior
      adj.M=adj.M*(adj.M>0)
      return(apply(adj.M,c(1,2),function(x) 1-prod(1-x)))
    }
  }
}

combine_anno_scores <- function(all_links, channel.use, prior=NULL){
  ## check if input similarities scores are all between 0 and 1
  minmax=range(all_links[,,channel.use])
  if(minmax[1]<0 | minmax[2]>1){stop("Abort. Similarity score range out of bound [0,1].")}
  if(length(channel.use)==1){
    if(is.null(prior)){
      return(all_links[,,channel.use])
    }else{
      adj.M=all_links[,,channel.use]-prior
      return(adj.M*(adj.M>0))
    }
  }else{
    if(is.null(prior)){
      return(apply(all_links[,,channel.use],c(1,2),function(x) 1-prod(1-x)))
    }else{
      adj.M=all_links[,,channel.use]-prior
      adj.M=adj.M*(adj.M>0)
      return(apply(adj.M,c(1,2),function(x) 1-prod(1-x)))
    }
  }
}