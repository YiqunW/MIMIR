require(abind)
require(igraph)
require(leiden)
integrate_exp_anno <- function(expM, annoM, method=c("AND","OR","+"), genes=c("inner","outer"), minmax=NULL){
  ## Test if both matrices contain the same genes
  if(!setequal(rownames(expM), rownames(annoM))){
    print("Genes in the two matrices don't match")
    if(genes=="inner"){
      g_use=intersect(rownames(expM), rownames(annoM))
      expM=expM[g_use, g_use]
      annoM=annoM[g_use, g_use]
    }else if(genes=="outer"){
      g_use=union(rownames(expM), rownames(annoM))
      M=matrix(0,nrow = length(g_use), ncol = length(g_use), dimnames = list(g_use, g_use))
      M[rownames(expM),colnames(expM)]=expM
      expM=M
      M=matrix(0,nrow = length(g_use), ncol = length(g_use), dimnames = list(g_use, g_use))
      M[rownames(annoM),colnames(annoM)]=annoM
    }else{
      stop("Sepcify which genes to include in the integrated matrix with genes='inner' or 'outer'")
    }
  }else{
    annoM=annoM[rownames(expM),colnames(expM)]
  }
  ## Can scale the similarity matrices before integration
  if(!is.null(minmax)){
    if(minmax[1]<0 | minmax[2]>1 | minmax[1]>=minmax[2]){
      stop("minmax = c(min, max), min < max, both need to be within the bounds of [0,1].")
    }else{
      expM=(expM/max(expM))*(minmax[2]-minmax[1])+minmax[1]
      annoM=(annoM/max(annoM))*(minmax[2]-minmax[1])+minmax[1]
    }
  }
  if(method=="AND"){
    simM=expM*annoM
    return(simM)
  }else if(method=="OR"){
    simM=1-(1-expM)*(1-annoM)
  }else if(method=="+"){
    simM=expM+annoM
  }else{stop("method needs to be 'AND','OR', or '+'")}
  ## force diagonal to 1
  diag(simM) <- 1
  return(simM)
}

integrate_all_exp_anno <- function(exp_sim, anno_sim, exp_use=NULL, anno_use=NULL, method="AND", minmax=c(0.05, 1), combine_genes="outer", prt_process=F){
  if(is.null(exp_use)){
    exp_use=dimnames(exp_sim[[3]])
  }
  if(is.null(anno_use)){
    anno_use=dimnames(anno_sim[[3]])
  }
  genes.use=dimnames(exp_sim)[[1]]
  sim_integ=NULL
  for(e in exp_use){
    expM=exp_sim[,,e]
    for(a in anno_use){
      score_name=paste0(a,":",method,":",e)
      if(prt_process){print(paste0("Calculating ",score_name))}
      annoM=anno_sim[,,a]
      if(is.null(sim_integ)){
        sim_integ=array(0,dim=c(length(genes.use),length(genes.use),1),
                        dimnames = list(genes.use, genes.use, score_name))
        simM=integrate_exp_anno(expM, annoM, method=method, minmax=minmax, genes=combine_genes)
        sim_integ[,,score_name]=simM[genes.use,genes.use]
      }else{
        simM=integrate_exp_anno(expM, annoM, method=method, minmax=minmax, genes=combine_genes)
        sim_integ=abind(sim_integ, simM[genes.use,genes.use], new.names = list(genes.use,genes.use,c(dimnames(sim_integ)[[3]], score_name)))
      }
    }
  }
  return(sim_integ)
}

igraph_cls <- function(adj_tbl, method=NULL,leiden_res=8, leiden_iter=-1, leiden_par="RBConfigurationVertexPartition",...){
  diag(adj_tbl) <- 0
  if(is.null(method)){
    cls_methods=c("walktrap","louvain","infomap","label_prop","leiden","fast_greedy") #"spinglass", too slow "optimal", has error
  }else{
    cls_methods=intersect(method,c("components","walktrap","louvain","infomap","label_prop","leiden","spinglass","eigen","optimal", "edge_btw","fast_greedy"))
  }
  cls_names=cls_methods
  if("leiden"%in%cls_methods){
    cls_names=setdiff(cls_names,"leiden")
    for(res in leiden_res){
      for(parti in leiden_par){
        cls_names=c(cls_names,paste("leiden",res,strsplit(parti,"Partition")[[1]],sep = "_"))
      }
    }
  }
  cls_tbl=matrix(NA, nrow = dim(adj_tbl)[1], ncol=length(cls_names), dimnames = list(rownames(adj_tbl), cls_names))
  g <- graph.adjacency(adj_tbl, mode="undirected", weighted=TRUE)
  if("walktrap"%in%cls_methods){
    walktrap=cluster_walktrap(g, weights = E(g)$weight, steps = 3, merges = T, modularity = TRUE, membership = TRUE)
    cls_tbl[walktrap$names,"walktrap"]=membership(walktrap)
  }
  if("components"%in%cls_methods){
    connected=clusters(g)
    cls_tbl[names(connected$membership),"components"]=membership(connected) 
  }
  if("louvain"%in%cls_methods){
    louvain=cluster_louvain(g, weights = NULL)
    cls_tbl[louvain$names,"louvain"]=louvain$membership
  }
  if("leiden"%in%cls_methods){
    for(res in leiden_res){
      for(parti in leiden_par){
        lei=leiden(adj_tbl, resolution_parameter=res, n_iterations = leiden_iter, partition_type = parti,...) #"RBERVertexPartition", c("RBConfigurationVertexPartition", "ModularityVertexPartition", "RBERVertexPartition", "CPMVertexPartition" "SurpriseVertexPartition")
        cls_tbl[rownames(adj_tbl),paste("leiden",res,strsplit(parti,"Partition")[[1]],sep = "_")]=lei
      }
    }
  }
  if("eigen" %in% cls_methods){
    l_eigen=cluster_leading_eigen(g)
    cls_tbl[l_eigen$names,"eigen"]=membership(l_eigen)
  }
  if("infomap"%in%cls_methods){
    infomap=cluster_infomap(g, e.weights = NULL, v.weights = NULL, nb.trials = 100,modularity = TRUE)
    cls_tbl[infomap$names,"infomap"]=membership(infomap)
  }
  if("label_prop"%in%cls_methods){
    prop=cluster_label_prop(g, weights = NULL, initial = NULL, fixed = NULL)
    cls_tbl[prop$names,"label_prop"]=membership(prop)
  }
  if("spinglass"%in%cls_methods){
    print("spinglass")
    if(!"components"%in%cls_methods){
      connected=clusters(g)
    }
    if(length(unique(membership(connected)))==1){
      spin=round(dim(adj_tbl)[1]/10)
      print("1 component")
      spinglass= cluster_spinglass(g, spins = spin, parupdate = FALSE, start.temp = 1, stop.temp = 0.01,
                                   cool.fact = 0.99, update.rule ="simple", gamma = 1, implementation = "orig")
      cls_tbl[spinglass$names,"spinglass"]=membership(spinglass)
    }else{
      print("multiple components")
      k=0
      for(i in unique(membership(connected))){
        genei=names(connected$membership)[which(membership(connected)==i)]
        if(length(genei)<5){
          cls_tbl[genei,"spinglass"]=k+1
          k=k+1
        }else{
          gi=graph.adjacency(adj_tbl[genei,genei], mode="undirected", weighted=TRUE)
          #print(length(genei))
          spin=max(round(length(genei)/10),2)
          spinglass = cluster_spinglass(gi, spins = spin, parupdate = FALSE, start.temp = 1, stop.temp = 0.01,
                                        cool.fact = 0.99, update.rule ="simple", gamma = 1, implementation = "orig")
          cls_tbl[spinglass$names,"spinglass"]=membership(spinglass)+k
          k=max(cls_tbl[spinglass$names,"spinglass"])
        }
      }
    }
  }
  if("optimal"%in%cls_methods){
    optimal=cluster_optimal(g)
    cls_tbl[optimal$names,"optimal"]=membership(optimal)
  }
  if("edge_btw"%in%cls_methods){
    weights=max(E(g)$weight)*1.1-E(g)$weight
    edge_btw=cluster_edge_betweenness(g, weights = weights, directed = F, edge.betweenness = F, merges = F,
                                      bridges = F, modularity = F, membership = TRUE)
    cls_tbl[edge_btw$names,"edge_btw"]=membership(edge_btw)
  }
  if("fast_greedy"%in%cls_methods){
    fast_greedy=cluster_fast_greedy(g, merges = F, modularity = F)
    cls_tbl[fast_greedy$names,"fast_greedy"]=membership(fast_greedy)
  }
  return(cls_tbl)
}