require(AnnotationDbi)
require(GOSemSim)
require(GO.db)
require(abind)

#' Convert the ParentChildTreeFile downloaded from interprot (https://www.ebi.ac.uk/interpro/download/) into a dataframe with two columns.
#'
#' @param fileName String. Path to the downloaded txt file.
#'               
#' @return A dataframe with two columns. The 'parent' column contains parent annotation ID, while the 'child' column conatins the corresponding child ID.
#'
#' @export
format_interpro_p2c <- function(fileName="../example_data/interpro_ParentChildTreeFile.txt"){
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
  
  ## format as a table for input to the similarity calculation function
  interp.p2c=stack(interp.offspring)
  interp.p2c=interp.p2c[,c(2,1)]
  colnames(interp.p2c)=c("parent","child")
  interp.p2c=as.data.frame(as.matrix(interp.p2c)) 
  return(interp.p2c)
}

#' Compute the information content of the annotation terms from a database. Modified from the computeIC function in the GOSemSim package.
#'
#' @param goAnno dataframe with at least two columns: "gene" - gene symbols, "id" - annotation terms associated with the gene.
#'               If a gene is associated with multiple annotations, it will have multiple rows in the dataframe, each row with 1 annotation id.
#'               For GO annotations, an additional column named "category" might be expected if the dataframe contains GO terms from multiple categories 
#'               (among MF, BP, and CC).
#' @param cat optional. "MF", "BP", or "CC". Required only if input is GO annotations with a "category" column. Specifies the category in which ICs of annotations should be calculated. 
#' @param Offsprings A dataframe or list specifying the hierarchical relationship between annotation terms included in the goAnno dataframe.
#'                   If Offsprings is a dataframe, the first column must contain the parent ids, and the second the child ids.
#'                   If Offsprings is a list, each entry should be named with the parent id and contain the children ids as a vector.
#'                 
#' @return A named numeric vector with information contents of each annotation ids.
#'
#' @export
computeIC.anno <- function(object, cat=NULL) {
  #goAnno, cat=NULL, Offsprings=NULL
  goAnno=object@gene.anno
  Offsprings=object@hierarchy.df
  if(dim(Offsprings)[1]==0){
    Offsprings=object@hierarchy.l
    if(length(Offsprings)==0){
      Offsprings=NULL
    }
  }
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
      print("Need to supply Offspring as a list or dataframe")
      return(NA)
    }else if(class(Offsprings)!="list"){
      Offsprings = offspring.list(Offsprings) ## format parent-child relation dataframe into a list (list("parent1"=children1, "parent2"=children2, ...))
      Offsprings = Offsprings[which(!is.na(Offsprings))]
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
  #return(IC)
  if(is.null(cat)){
    object@IC=IC
  }else{
    object@cat.IC[[cat]]=IC
  }
  return(object)
}

#' Takes a parent-child dataframe and convert into a list.
#'
#' @param p2c dataframe with two columns. The first column contains the parent ids, the second contains the child ids.
#'                 
#' @return A list where each entry is named with the parent id and contains the children ids as a vector.
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

#' Compute the functional similarity between two genes based on their annotations from one functional annotation database. Modified from functions in the GOSemSim package.
#'
#' @param gene1, a gene symbol that can be found in the "gene" column in the goAnno dataframe.
#' @param gene2, a gene symbol that can be found in the "gene" column in the goAnno dataframe.
#' @param goAnno dataframe with at least two columns: "gene" - gene symbols, "id" - annotation terms associated with the gene.
#'               If a gene is associated with multiple annotations, it will have multiple rows in the dataframe, each row with 1 annotation id.
#'               For GO annotations, an additional column named "category" might be expected if the dataframe contains GO terms from multiple categories 
#'               (among MF, BP, and CC).
#' @param cat optional. "MF", "BP", or "CC". Required only if input is GO annotations with a "category" column. Specifies the category in which similarities should be calculated. 
#' @param IC named vector. names are functional annotation ids that matches those in the "id" column of the goAnno dataframe.
#'           Values are each id's information content (calculated with computeIC.anno()).
#' @param anno.p2c dataframe with two columns. The first column contains the parent annotation ids, the second contains the child ids.
#'                 If set to NULL, either anc.all need to be provided as a list (can be calculated with the ancester.list() function),
#'                 or GO annotation must be used with a specific cat. This provides the hierarchical relationship used in calculating semantic similarities. 
#' @param anc.all A list of annotation ids. each entry is named with an id and contains a vector of all its parent ids. 
#'                Not required (set to NULL) when anno.p2c is provided or if GO annotation is used and cat is one of "MF", "BP", and "CC".                  
#' 
#' @return A list of three entries: "geneSim" - similarity score between gene1 and gene2; 
#'         "Gene1Annotations" - annotation ids associated with gene1; 
#'         "Gene2Annotations" - annotation ids associated with gene2; 
#'
#' @export
geneSim.anno <- function(gene1, gene2, object, measure=c("Jiang","Wang"), combine="BMA",cat=NULL){
  #goAnno, IC, ,anno.p2c=NULL,anc.all=NULL
  ## IC should be named vector
  goAnno=object@gene.anno
  anno.p2c=object@hierarchy.df
  anc.all=NULL
  if(dim(anno.p2c)[1]==0){
    anno.p2c=NULL
    anc.all=object@hierarchy.l
    if(length(anc.all)==0){
      anc.all=NULL
    }
  }
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
    return (list(geneSim=NA, Gene1Annotations=go1, Gene2Annotations=go2))
  }
  if(is.null(cat)){
    IC=object@IC
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
  }else{
    IC=object@cat.IC[[cat]]
  }
  res <- annoSim(go1, go2, goAnno=goAnno, IC=IC, measure=measure, combine=combine, cat=cat, anc=anc)
  return (list(geneSim=res, Gene1Annotations=go1, Gene2Annotations=go2))
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


library(parallel)
#' Compute the pairwise gene functional similarities 
#'
#' @param gene_anno dataframe with at least two columns: "gene" - gene symbols, "id" - annotation terms associated with the gene.
#'                  If a gene is associated with multiple annotations, it will have multiple rows in the dataframe, each row with 1 annotation id.
#'                  For GO annotations, an additional column named "category" might be expected if the dataframe contains GO terms from multiple categories 
#'                  (among MF, BP, and CC).
#' @param cat optional. "MF", "BP", or "CC". Required only if input is GO annotations with a "category" column. Specifies the category in which similarities should be calculated. 
#' @param IC named vector. names are functional annotation ids that matches those in the "id" column of the goAnno dataframe.
#'           Values are each id's information content (calculated with computeIC.anno()).
#' @param p2c dataframe with two columns. The first column contains the parent annotation ids, the second contains the child ids.
#'            If set to NULL, GO annotation must be used with cat equals one of "MF", "BP", and "CC". 
#'            This provides the hierarchical relationship used in calculating semantic similarities. 
#' @param genes A vector of genes for which pairwise functional similarities will be calculated. 
#'              Should be a subset of the genes in the "gene" column of the gene_anno dataframe.             
#' 
#' @return A square matrix (gene x gene) with gene names as rownames and colnames, and functional similarities as values.
#'
#' @export
func_sim <- function(object, cat=NULL, genes=NULL){
  #gene_anno, p2c, IC
  if(is.null(genes)){
    genes=union(object@genes, object@bg.genes)
  }else{
    genes=intersect(genes,object@gene.anno$gene)
  }
  gene_anno=object@gene.anno
  
  my_function <- function(gene1, gene2, object.use=object, cat.use=cat){
    res <- geneSim.anno(gene1,gene2, object=object.use, measure="Jiang", cat=cat.use, combine = "BMA")
    return(res$geneSim)
  }
  n_genes <- length(genes)
  # Define a helper function for parallel computation
  compute_row <- function(i) {
    sapply(i:n_genes, function(j) my_function(genes[i], genes[j]))
  }
  # Use mclapply for parallel row computation (on Unix-like systems)
  result_list <- mclapply(1:n_genes, compute_row, mc.cores = detectCores()-1)
  
  anno.sim <- matrix(0, nrow = n_genes, ncol = n_genes)
  # Fill the lower triangle with the vectors in result list
  anno.sim[lower.tri(anno.sim, diag = TRUE)] <- unlist(result_list)
  anno.sim <- anno.sim + t(anno.sim)
  rownames(anno.sim)=genes
  colnames(anno.sim)=genes
  
  # Clean up: assign 1 to diagonal, make sure values are between 0 and 1
  diag(anno.sim) <- 1
  anno.sim[is.na(anno.sim)] <- 0
  anno.sim[anno.sim<0] <- 0
  anno.sim[anno.sim>1] <- 1
  if(is.null(cat)){
    object@similarity=anno.sim
    object@prior = mean(anno.sim[upper.tri(anno.sim)])
  }else{
    object@cat.similarity[[cat]]=anno.sim
    object@cat.prior[[cat]] = mean(anno.sim[upper.tri(anno.sim)])
  }
  return(object)
}

#' Calculate similarity scores between two sets of annotation ids. Adapted from functions from the GOSemSim package.
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

#' Re-compute gene-gene interaction scores using custom evidence channels from the STRING database  
#'
#' @param str_tbl_full dataframe downloaded from the STRING database website (https://string-db.org/cgi/download.pl).
#' @param protein_cols a vector of length 2 specifying the names of the two columns containing gene (or protein) ids.
#' @param evi a vector specifying the evidence channels to use for computing combined interaction scores. 
#'            Should be a subset of the colnames in the str_tbl_full (excluding "combined_score").
#' 
#' @return A dataframe with 3 columns: gene1, gene2 (colnames same as specified in protein_cols), and combined interaction score ("combined_score").
#'
#' @export
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

## Convert the STRING database gene-gene interaction table to a symmetric similarity matrix
library(dplyr)
#library(igraph)
str_ppi_2tbl <- function(str_ppi, thres, values="combined_score",g1.col="protein1", g2.col="protein2",method="max"){
  ind_kp=which(str_ppi[,values]>thres)
  str_ppi=str_ppi[ind_kp,]
  all_g=union(str_ppi[[g1.col]],str_ppi[[g2.col]])
  str_ppi=str_ppi[,c(g1.col, g2.col, values)]
  colnames(str_ppi)=c("gene1","gene2","combined_score")
  if(method=="max"){
    unique_values <- str_ppi %>% group_by(gene1, gene2) %>% summarize(Score_use = max(combined_score),.groups = 'drop')
  }else if(method=="mean"){
    unique_values <- str_ppi %>% group_by(gene1, gene2) %>% summarize(Score_use = mean(combined_score),.groups = 'drop')
  }else if(method=="median"){
    unique_values <- str_ppi %>% group_by(gene1, gene2) %>% summarize(Score_use = median(combined_score),.groups = 'drop')
  }

  # g <- graph_from_data_frame(d = unique_values, directed = FALSE)
  # str_ppi.M <- as.matrix(as_adjacency_matrix(g, attr = "Score_use"))
  
  str_ppi.M=matrix(0,nrow=length(all_g), ncol = length(all_g), dimnames = list(sort(all_g),sort(all_g)))
  for(i in 1:dim(unique_values)[1]){
    p1=as.character(unique_values[i,g1.col])
    p2=as.character(unique_values[i,g2.col])
    str_ppi.M[p1,p2]=unique_values[["Score_use"]][i]
    str_ppi.M[p2,p1]=unique_values[["Score_use"]][i]
  }
  all.zeros=which(apply(str_ppi.M,1,sum)==0)
  if(length(all.zeros)>0){
    str_ppi.M=str_ppi.M[-all.zeros,-all.zeros]
  }
  return(str_ppi.M)
}

str_sim_M <- function(string_tbl, evi, genes_use, protein_cols=c("gene1","gene2"), score_col="combined_score"){
  if(!is.null(genes_use)){
    ind1=which(string_tbl[,protein_cols[1]]%in%genes_use)
    ind2=which(string_tbl[,protein_cols[2]]%in%genes_use)
    ind=intersect(ind1,ind2)
    if(length(ind)==0){
      stop("String table doesn't contain input genes in genes_use.")
    }
    string_tbl=string_tbl[ind,]
  }
  if(length(evi)==1 && evi=="default"){
    string_tbl[score_col]=string_tbl[score_col]/1000
    return(str_ppi_2tbl(string_tbl,thres = 0,values = score_col, g1.col=protein_cols[1], g2.col=protein_cols[2], method="max"))
  }else{
    if(length(evi)==1){
      str_combined=string_tbl[evi]/1000
    }else{
      str_combined=str_comb(string_tbl, protein_cols=protein_cols, evi=evi)
      str_combined=str_combined["combined_score"]
    }
    string_tbl[score_col]=str_combined
    return(str_ppi_2tbl(string_tbl,thres = 0,values = score_col, g1.col=protein_cols[1], g2.col=protein_cols[2], method="max"))
  }
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

stack_func_sim <- function(simM_list, genes_use="all", adj_prior=F, add_to=NULL, verbose=T){
  if(is.null(genes_use)){
    if(!is.null(add_to)){
      genes_use=dimnames(add_to)[[1]]
    }else{
      message("Using union of genes across input matrices...")
      genes_use="all"
    }
  }
  if(length(genes_use)==1 && genes_use=="all"){
    genes_use=c()
    for(M in simM_list){
      genes_use=unique(c(genes_use,rownames(M)))
    }
    if(!is.null(add_to)){
      genes_use=unique(c(genes_use, dimnames(add_to)[[1]]))
    }
  }else{
    if(!is.null(add_to)){
      if(length(intersect(genes_use, dimnames(add_to)[[1]]))==0){
        stop("There's no overlap between genes_use and genes in the add_to matrix. Abort.")
      }
    }
  }
  genes_use=sort(genes_use)
  anno_names=names(simM_list)
  anno3d=array(0,dim=c(length(genes_use),length(genes_use),length(anno_names)),
               dimnames = list(genes_use, genes_use, anno_names))
  for(anno in anno_names){
    if(verbose){
      print(paste0('Adding ', anno))
    }
    simM=simM_list[[anno]]
    if(is.null(adj_prior) || adj_prior){
      if(is.numeric(adj_prior)){
        simM=subtract_mean(simM, prior=adj_prior)
      }else{
        simM=subtract_mean(simM, prior=NULL, prt_mean = verbose)
      }
    }
    comm.genes=intersect(genes_use,rownames(simM))
    anno3d[comm.genes, comm.genes, anno]=simM[comm.genes, comm.genes]
  }
  if(is.null(add_to)){
    if(verbose){
      print(paste0("Result array dimensions: ", paste(dim(anno3d),collapse = " x ")))
      print(paste0("Sources used: ",paste(dimnames(anno3d)[[3]], collapse = ", ")))
      print(paste(c("Value range:", range(anno3d)), collapse = " "))
    }
    return(anno3d)
  }else{
    add_to_use=array(0,dim=c(length(genes_use),length(genes_use),dim(add_to)[3]),
                     dimnames = list(genes_use, genes_use, dimnames(add_to)[[3]]))
    comm.genes=intersect(dimnames(add_to)[[1]], genes_use)
    add_to_use[comm.genes,comm.genes, dimnames(add_to)[[3]]] = add_to[comm.genes,comm.genes,]
    if(verbose){
      print(paste0("Result array dimensions: ", paste(dim(anno3d),collapse = " x ")))
      print(paste0("Sources used: ",paste(dimnames(anno3d)[[3]], collapse = ", ")))
      print(paste(c("Value range:", range(anno3d)), collapse = " "))
    }
    return(abind(add_to_use, anno3d))
  }
}

combine_anno_scores <- function(anno.sim, channel.use, prior=NULL){
  ## check if input similarities scores are all between 0 and 1
  minmax=range(anno.sim[,,channel.use])
  if(minmax[1]<0 | minmax[2]>1){stop("Abort. Similarity score range out of bound [0,1].")}
  if(length(channel.use)==1){
    if(is.null(prior)){
      return(anno.sim[,,channel.use])
    }else{
      adj.M=anno.sim[,,channel.use]-prior
      return(adj.M*(adj.M>0))
    }
  }else{
    if(is.null(prior)){
      return(apply(anno.sim[,,channel.use],c(1,2),function(x) 1-prod(1-x)))
    }else{
      adj.M=anno.sim[,,channel.use]-prior
      adj.M=adj.M*(adj.M>0)
      return(apply(adj.M,c(1,2),function(x) 1-prod(1-x)))
    }
  }
}

# Combine scores from different annotation sources
# 
# This function takes functional similarities derived from different annotation databases and combines them.
# 
# Arguments:
# - anno.sim: A 3d array of size n_genes x n_genes x n_databases. Storing pairwise functional similarities derived from each database.
# - how.to: A named list of vectors specifying which scores to combine (names matches the 3rd dimnames of the input score array). Entry names will be used as the 3rd dimension names of the array storing the combined scores.
# - add: Boolean. If TRUE, bind the newly computed combined score array with the input anno.sim array and return this extended score array. If FALSE, return the newly computed combined score array only.
# 
# Returns:
# - A 3d array with combined (and individual, if add=T) similarity scores.
# 
# Examples:
# add_combined_scores(anno.sim, how.to = list("combination1"=c("GO_BP","Interpro"), "combination2"=c("GO_BP","GO_MF")))
add_combined_scores <- function(anno.sim, how.to, add=T, verbose=T){
  anno_combined=list()
  for(comb_name in names(how.to)){
    anno_combined[[comb_name]]=combine_anno_scores(anno.sim, how.to[[comb_name]])
  }
  combined_anno=do.call(abind,c(anno_combined,along=3))
  if(add){
    combined_anno=abind(anno.sim, combined_anno)
  }
  if(verbose){
    print(paste0("Result array dimensions: ", paste(dim(combined_anno),collapse = " x ")))
    print(paste0("Sources used: ",paste(dimnames(combined_anno)[[3]], collapse = ", ")))
    print(paste(c("Value range:", range(combined_anno)), collapse = " "))
  }
  return(combined_anno)
}

## define a function to add the functional similarity array to mimir object
add_to_MIMIR <- function(mimir.obj, anno.sim){
  gene_order=rownames(mimir.obj@exp.sim) #make sure the new array is in compatible with the expression similarity array already in the object
  mimir.obj@func.sim <- anno.sim[gene_order, gene_order,]
  return(mimir.obj)
}

## define a function to add the gene-annotation-description table in the anno objects to mimir object
add_anno_to_MIMIR <- function(mimir.obj, anno.obj, filter.genes=T){
  if(filter.genes){
    genes.keep.ind=which(anno.obj@gene.anno$gene%in%mimir.obj@genes)
    mimir.obj@annotations[[anno.obj@db]]=anno.obj@gene.anno[genes.keep.ind,]
  }else{
    mimir.obj@annotations[[anno.obj@db]]=anno.obj@gene.anno
  }
  return(mimir.obj)
}
