library(mclust)
library(aricode)
#' Compare Clustering Results
#'
#' This function calculates AMI and other metrics from the clustering results.
#'
#' @param clus_s A matrix where each column stores cluster membership from a clustering method.
#' @param module_tbl A dataframe with columns 'gene' and 'level', representing the true module assignments.
#'
#' @return A matrix with AMI values for each clustering method.
#'
#' @export
comp_par <- function(clus_s, module_tbl){
  g_use=intersect(rownames(clus_s),module_tbl$gene)
  module_tbl=module_tbl[module_tbl[,"gene"]%in%g_use,]
  ## initiate a matrix to store results
  results_m=matrix(0,nrow=1,ncol=dim(clus_s)[2])
  colnames(results_m)=colnames(clus_s) #each column in clus_s stores cluster membership from a clustering method
  rownames(results_m)=c("AMI")
  for(c in colnames(results_m)){
    ## calculate AMI
    results_m["AMI",c]=AMI(module_tbl[,'level'], clus_s[module_tbl[,'gene'],c])
  }
  return(results_m)
}

#' Fisher's Exact Test for Enrichment
#'
#' This function performs a Fisher's exact test to identify enriched annotations.
#'
#' @param genes A vector of genes to test for enrichment.
#' @param gXm A binary matrix indicating gene-annotation relationships.
#' @param p.adj.method Method for adjusting p-values. Default is "fdr".
#' @param p.adj_only Boolean. If TRUE, only adjusted p-values are returned.
#'
#' @return A dataframe with enrichment results.
#'
#' @export
fisher.exact.enrichment <- function(genes, gXm, p.adj.method="fdr", p.adj_only=FALSE) {
  # find the annotations in the group
  gXm=gXm[which(apply(gXm,1,sum)>0),]
  gXm=gXm==1
  genes=intersect(genes,rownames(gXm))
  if(length(genes)>0){
    annos.test=colnames(gXm)[which(apply(gXm[genes,,drop=F],2,sum)>0)]
    # pre-allocate dataframe to store test results
    enrich=data.frame(row.names = annos.test, stringsAsFactors = F)
    enrich$p=NA
    if(!p.adj_only){
      enrich$genes.in.selection=""
      enrich$genes.in.tissue=""
    }
    all.genes=rownames(gXm)
    out.genes=setdiff(all.genes,genes)
    # for each annotation, calculate enrichment
    for(anno in annos.test){
      all.with.anno = all.genes[gXm[,anno]]
      in.with.anno = intersect(all.with.anno, genes)
      out.with.anno = intersect(all.with.anno, out.genes)
      in.without.anno = setdiff(genes, in.with.anno)
      out.without.anno = setdiff(out.genes, out.with.anno)
      a=length(in.with.anno)
      b=length(out.with.anno)
      c=length(in.without.anno)
      d=length(out.without.anno)
      ft <- fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2),alternative="greater")
      enrich[anno,"p"]=ft$p.value
      if(!p.adj_only){
        enrich[anno,c("genes.in.selection","genes.in.tissue")]=c(paste0(in.with.anno,collapse = ", "),paste0(out.with.anno,collapse = ", "))
      }
    }
    # calculate fdr adjusted p values
    enrich$p.adj=p.adjust(enrich$p, p.adj.method)
    enrich=enrich[order(enrich$p.adj,decreasing = F),]
    return(enrich)
  }else{
    return(NULL)
  }
}

#' Convert Dataframe to Binary Matrix
#'
#' This function converts a dataframe to a binary matrix based on specified rows and columns.
#'
#' @param anno.tbl A dataframe with annotation information.
#' @param rows The name of column to use as rows in the matrix.
#' @param cols The name of column to use as columns in the matrix.
#'
#' @return A binary matrix.
#'
#' @export
df2matrix <- function(anno.tbl, rows="gene", cols="id"){
  genes=unique(anno.tbl[,rows])
  ids=unique(anno.tbl[,cols])
  genexanno=matrix(0,nrow=length(genes), ncol=length(ids), dimnames=list(genes,ids))
  genexanno[cbind(anno.tbl[,rows],anno.tbl[,cols])]=1
  return(genexanno)
}

#' Assess Clustering Results
#'
#' This function assesses the performance of different clustering methods using various metrics.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param manual_modules Optional. A dataframe with true module assignments for comparison.
#' @param return_df_only Boolean. If TRUE, only the dataframe of results is returned.
#' @param internal_db_use Optional. Name of the internal annotation database to use.
#' @param go_cat Optional. GO categories to include if using internal database.
#' @param external_db_tbl Optional. A dataframe with external annotation data.
#' @param verbose Boolean. If TRUE, prints progress messages.
#'
#' @return An object with updated cluster metrics or a dataframe of results.
#'
#' @export
assess_result <- function(object, manual_modules=NULL, return_df_only=F, internal_db_use=NULL, go_cat=NULL, external_db_tbl=NULL, verbose=T){
  exp_sim=object@exp.sim
  exp_sim_names=dimnames(exp_sim)[[3]]
  anno_sim=object@func.sim
  anno_sim_names=dimnames(anno_sim)[[3]]
  clus_all=object@clusters
  ## if correlations between genes not calculated before, do that
  if(!"correlation"%in%dimnames(object@exp.sim[[3]])){
    correlation=cor(t(object@exp.data[object@genes,]))
    object@exp.sim=abind(object@exp.sim,correlation, along=3,make.names=T)
  }
  corM=object@exp.sim[,,'correlation']
  ## Construct a table for all the metrics
  result_list=list()
  metrics=c("n_cluster","n_singlet","max_cluster", "intra_mean_exp_cor", "pct_gene_in_cluster",
            "similarity", "similarity_mode", "anno_sim", "exp_sim", "cluster_method","leiden_res")
  if(!is.null(manual_modules)){
    metrics=c(metrics, "AMI")
  }
  if(!is.null(internal_db_use)){
    metrics=c(metrics, paste0(internal_db_use,"_mean(-logP)"),paste0(internal_db_use,"_median(-logP)"))
    db_in=object@annotations[[internal_db_use]]
    if(!is.null(go_cat)&&("category"%in%colnames(db_in))){
      db_in=db_in[db_in$category%in%go_cat,]
    }
    genexanno_in=df2matrix(db_in, rows="gene", cols="id")
  }
  if(!is.null(external_db_tbl)){
    exdb_name=deparse(substitute(external_db_tbl))
    metrics=c(metrics, paste0(exdb_name,"_mean(-logP)"),paste0(exdb_name,"_median(-logP)"))
    genexanno_ex=df2matrix(external_db_tbl, rows="gene", cols="id")
  }
  assess_tbl<-data.frame(row.names = metrics, stringsAsFactors = F)
  assess_tbl=as.data.frame(t(assess_tbl))
  for(s in names(clus_all)){
    print(paste0("Calculating stats for clustering results using ", s))
    clus_s=clus_all[[s]]
    ngenes=dim(clus_s)[1]
    add_tbl=matrix(NA, nrow=length(metrics), ncol = dim(clus_s)[2])
    rownames(add_tbl)=metrics
    colnames(add_tbl)=colnames(clus_s)
    add_tbl=as.data.frame(t(add_tbl))
    add_tbl['similarity']=s
    ## separate the similarity name field to get which expression and/or functional similarities were used
    ## separator character ":" is hard coded in the integrate_all_exp_anno() function
    if(grepl(":", s,fixed = T)){
      #functional and expression similarities were combined
      scores=unlist(strsplit(s,":"))
      add_tbl['similarity_mode']=scores[2] #this is either "AND", "OR", or "+"
      add_tbl['anno_sim']=scores[1]
      add_tbl['exp_sim']=scores[3]
    }else{
      #only functional or only expression similarities were used
      if(s%in%anno_sim_names){
        add_tbl['similarity_mode']="function only"
        add_tbl['anno_sim']=s
        add_tbl['exp_sim']="none"
      }else if(s%in%exp_sim_names){
        add_tbl['similarity_mode']="expression only"
        add_tbl['anno_sim']="none"
        add_tbl['exp_sim']=s
      }
    }
    if(!is.null(manual_modules)){
      add_tbl["AMI"]=as.vector(comp_par(clus_s, manual_modules))
    }
    clus_m=strsplit(colnames(clus_s),"_")
    add_tbl["cluster_method"]=unlist(lapply(clus_m,function(x) x[1]))
    add_tbl["leiden_res"]=unlist(lapply(clus_m,function(x) as.numeric(x[2]))) #non-leiden methods will have NA leiden_res
    membs=apply(clus_s,2,table)
    add_tbl["n_cluster"]=unlist(lapply(membs, length))
    add_tbl["n_singlet"]=unlist(lapply(membs, function(x){return(sum(x==1))}))
    add_tbl["max_cluster"]=unlist(lapply(membs, max))
    add_tbl['pct_gene_in_cluster']=100*(ngenes-add_tbl['n_singlet'])/ngenes
    ## fill in other metrics
    for(c in colnames(clus_s)){
      clus.list=split(row.names(clus_s), clus_s[,c])
      clus.list=Filter(function(x) length(x)>1, clus.list)
      mean_cor=lapply(clus.list, function(x){return(mean(corM[x,x][lower.tri(corM[x,x])]))})
      add_tbl[c,"intra_mean_exp_cor"]=mean(unlist(mean_cor))
      if(!is.null(internal_db_use)){
        g.enrich=lapply(clus.list, fisher.exact.enrichment,gXm=genexanno_in, p.adj_only=T) #calculate enriched terms in each cluster
        top_enrich=unlist(lapply(g.enrich,function(x){if(!is.null(x)){return(-log(x[1,'p.adj']))}else{return(0)}})) #get -log p-val for the most enriched term in each cluster
        top_enrich[!is.finite(top_enrich)]=40 # an arbitrary large number
        add_tbl[c, paste0(internal_db_use,"_mean(-logP)")]=mean(top_enrich)
        add_tbl[c, paste0(internal_db_use,"_median(-logP)")]=median(top_enrich)
      }
      if(!is.null(external_db_tbl)){
        g.enrich=lapply(clus.list, fisher.exact.enrichment,gXm=genexanno_ex, p.adj_only=T) #calculate enriched terms in each cluster
        top_enrich=unlist(lapply(g.enrich,function(x){if(!is.null(x)){return(-log(x[1,'p.adj']))}else{return(0)}})) #get -log p-val for the most enriched term in each cluster
        top_enrich[!is.finite(top_enrich)]=40 # an arbitrary large number
        add_tbl[c, paste0(exdb_name,"_mean(-logP)")]=mean(top_enrich)
        add_tbl[c, paste0(exdb_name,"_median(-logP)")]=median(top_enrich)
      }
    }
    assess_tbl=rbind(assess_tbl,add_tbl)
  }
  rownames(assess_tbl)=1:dim(assess_tbl)[1]
  if(return_df_only){
    return(assess_tbl)
  }else{
    object@cluster_metrics <- assess_tbl
    return(object)
  }
}

#' Add External Assessment Results
#'
#' This function adds external assessment results to the existing cluster metrics.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param return_df_only Boolean. If TRUE, only the dataframe of results is returned.
#' @param external_db_tbl A dataframe with external annotation data.
#' @param verbose Boolean. If TRUE, prints progress messages.
#'
#' @return An object with updated cluster metrics or a dataframe of results.
#'
#' @export
add_external_assess <- function(object, return_df_only=F, external_db_tbl=NULL, verbose=T){
  if(dim(object@cluster_metrics)[1]==0){
    stop("No cluster_metrics calculated yet. Call assess_result() function first.")
  }
  exdb_name=deparse(substitute(external_db_tbl))
  clus_all=object@clusters
  ## Construct a table for all the metrics
  result_list=list()
  metrics=c('similarity', "cluster_method","leiden_res", paste0(exdb_name,"_mean(-logP)"),paste0(exdb_name,"_median(-logP)"))
  genexanno_ex=df2matrix(external_db_tbl, rows="gene", cols="id")
  
  assess_tbl<-data.frame(row.names = metrics, stringsAsFactors = F)
  assess_tbl=as.data.frame(t(assess_tbl))
  for(s in names(clus_all)){
    print(paste0("Calculating stats for clustering results using ", s))
    clus_s=clus_all[[s]]
    add_tbl=matrix(NA, nrow=length(metrics), ncol = dim(clus_s)[2])
    rownames(add_tbl)=metrics
    colnames(add_tbl)=colnames(clus_s)
    add_tbl=as.data.frame(t(add_tbl))
    add_tbl['similarity']=s
    clus_m=strsplit(colnames(clus_s),"_")
    add_tbl["cluster_method"]=unlist(lapply(clus_m,function(x) x[1]))
    add_tbl["leiden_res"]=unlist(lapply(clus_m,function(x) as.numeric(x[2]))) #non-leiden methods will have NA leiden_res
    ## fill in other metrics
    for(c in colnames(clus_s)){
      clus.list=split(row.names(clus_s), clus_s[,c])
      clus.list=Filter(function(x) length(x)>1, clus.list)
      g.enrich=lapply(clus.list, fisher.exact.enrichment,gXm=genexanno_ex, p.adj_only=T) #calculate enriched terms in each cluster
      top_enrich=unlist(lapply(g.enrich,function(x){if(!is.null(x)){return(-log(x[1,'p.adj']))}else{return(0)}})) #get -log p-val for the most enriched term in each cluster
      top_enrich[!is.finite(top_enrich)]=40 # an arbitrary large number
      add_tbl[c, paste0(exdb_name,"_mean(-logP)")]=mean(top_enrich)
      add_tbl[c, paste0(exdb_name,"_median(-logP)")]=median(top_enrich)
    }
    assess_tbl=rbind(assess_tbl,add_tbl)
  }
  rownames(assess_tbl)=1:dim(assess_tbl)[1]
  if(return_df_only){
    return(assess_tbl)
  }else{
    rownames(assess_tbl)=paste(assess_tbl[,'similarity'], assess_tbl[,"cluster_method"],assess_tbl[,"leiden_res"])
    row_match=paste(object@cluster_metrics[,'similarity'], object@cluster_metrics[,"cluster_method"],object@cluster_metrics[,"leiden_res"])
    object@cluster_metrics <- cbind(object@cluster_metrics, assess_tbl[row_match,c(paste0(exdb_name,"_mean(-logP)"),paste0(exdb_name,"_median(-logP)"))])
    return(object)
  }
}

library(plotly)
#' Add Vertical Line to Plotly Plot
#'
#' @param x X-coordinate for the vertical line.
#' @param color Color of the vertical line. Default is "grey".
#'
#' @return A list specifying the vertical line.
#'
#' @export
vline <- function(x = 0, color = "grey") {
  return(list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash="dot")
  ))
}

#' Add Horizontal Line to Plotly Plot
#'
#' @param y Y-coordinate for the horizontal line.
#' @param color Color of the horizontal line. Default is "grey".
#'
#' @return A list specifying the horizontal line.
#'
#' @export
hline <- function(y = 0, color = "grey") {
  return(list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash="dot")
  ))
}

#' Create a Plotly Scatter Plot
#'
#' This function creates a scatter plot using Plotly with various customization options.
#'
#' @param data A dataframe containing the data to plot.
#' @param x The name of the column to use for the x-axis.
#' @param y The name of the column to use for the y-axis.
#' @param color Optional. The name of the column to use for coloring the points.
#' @param hover_text Optional. The name of the column to use for hover text.
#' @param marker_opacity Opacity of the markers. Default is 0.5.
#' @param marker_size Size of the markers. Default is 5.
#' @param color_pal Color palette to use. Default is "Set1".
#' @param addline_x Optional. X-coordinate for a vertical line to add.
#' @param addline_y Optional. Y-coordinate for a horizontal line to add.
#' @param ylab Optional. Label for the y-axis.
#' @param xlab Optional. Label for the x-axis.
#' @param figsize Optional. A vector specifying the width and height of the plot.
#'
#' @return A Plotly scatter plot.
#'
#' @export
plotly_scatter <- function(data, x, y, color=NULL, hover_text=NULL, marker_opacity=0.5, marker_size=5, 
                           color_pal="Set1", addline_x=NULL, addline_y=NULL, ylab=NULL, xlab=NULL, figsize=NULL){
  if(is.null(color)){
    fig <- plot_ly(x=data[[x]], y=data[[y]], type="scatter", mode="markers",text = data[[hover_text]], 
                   marker = list(opacity = marker_opacity,size=marker_size),colors=color_pal)
  }else{
    fig <- plot_ly(x=data[[x]], y=data[[y]], color = data[[color]], type="scatter", mode="markers", text = data[[hover_text]],
                   marker = list(opacity = marker_opacity,size=marker_size),colors=color_pal)
    tryCatch({fig <- fig  %>% colorbar(title = color)}, error=function(e){NULL}, warning=function(w){NULL})
  }
  line_list=NULL
  if(!is.null(addline_y)){
    line_list=c(line_list,list(hline(addline_y)))
  }
  if(!is.null(addline_x)){
    line_list=c(line_list,list(vline(addline_x)))
  }
  if(!is.null(line_list)){
    fig <- fig %>% layout(shapes = line_list)
  }
  if(is.null(xlab)){
    xlab=x
  }
  if(is.null(ylab)){
    ylab=y
  }
  fig <- fig %>% layout(xaxis = list(title=xlab), 
                        yaxis = list(title=ylab),
                        legend = list(title = list(text = color)))
  if(!is.null(figsize)){
   fig <- fig %>%layout(width = figsize[1],  height = figsize[2])
  }
  return(fig)
}

#' Get Cluster Parameters
#'
#' This function extracts similarity and clustering method parameters from a given cluster.
#'
#' @param cluster_use A dataframe containing cluster information.
#'
#' @return A vector with similarity and clustering method parameters.
#'
#' @export
get_cluster_pars <- function(cluster_use){
  similarity_use=cluster_use$similarity
  leiden_res_use=cluster_use$leiden_res
  if(is.na(leiden_res_use)){
    cluster_method=cluster_use$cluster_method
  }else{
    cluster_method=paste0('leiden_',leiden_res_use, '_')
  }
  return(c(similarity_use, cluster_method))
}

#' Get Clusters from Object
#'
#' This function retrieves clusters from an object based on specified parameters.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param cluster_par A vector with similarity and clustering method parameters.
#'
#' @return A dataframe with cluster assignments.
#'
#' @export
get_clusters <- function(object, cluster_par){
  if(length(cluster_par)>2){
    cluster_par=get_cluster_pars(cluster_par)
  }
  clus.df=object@clusters[[cluster_par[1]]]
  clus=clus.df[,grepl(paste0("^",cluster_par[2]), colnames(clus.df)),drop=F]
  return(clus)
}

#' Add Cluster Results to Object
#'
#' This function adds cluster results to an object based on specified metrics.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param metric_use A dataframe with metrics to use for adding cluster results.
#'
#' @return An updated object with added cluster results.
#'
#' @export
add_cluster_result <- function(object, metric_use){
  for(i in 1:dim(metric_use)[1]){
    cluster_par=get_cluster_pars(metric_use[i,,drop=F])
    cluster_use=get_clusters(object, cluster_par)
    colnames(cluster_use)=paste0(cluster_par[1],";",colnames(cluster_use))
    object@gene.info[rownames(cluster_use),colnames(cluster_use)]=cluster_use
  }
  return(object)
}

plot_to_grob <- function(plot_func) {
  # Create a temporary file
  tmpfile <- tempfile(fileext = ".png")
  
  # Create the plot in the temporary file
  png(tmpfile, width = 400, height = 400)
  plot_func()
  dev.off()
  
  # Read the image back into R
  img <- readPNG(tmpfile)
  img_grob <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  
  return(img_grob)
}

library(patchwork)
library(Vennerable)
library(cowplot)
library(png)
#' Plot Clusters
#'
#' This function generates plots for clusters using expression data and enrichment analysis.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param cluster_use The name of the cluster to plot.
#' @param exp_use The expression data to use for plotting. Default is "smoothed.exp".
#' @param ncols Number of columns in the plot layout. Default is 5.
#' @param use_description Boolean. If TRUE, use descriptions for enrichment analysis. Default is TRUE.
#' @param save_pdf Boolean or string. If TRUE, save the plot as a PDF. If string, use the string as the file name.
#' @param save_width Width of the saved plot. Default is 20.
#' @param save_height Height of the saved plot. Default is 40.
#' @param return_plots Boolean. If TRUE, return the list of plots instead of displaying them.
#' @param show_plot Boolean. If TRUE, display the plot.
#'
#' @return A list of plots if return_plots is TRUE.
#'
#' @export
plotClusters <- function(object, cluster_use, exp_use=c("smoothed.exp","exp.data"),ncols=5,use_description=T,
                         save_pdf=F, save_width=20, save_height=40, return_plots=F,show_plot=T){
  if(is.character(save_pdf)){
    if(!grepl("\\.pdf$", save_pdf, ignore.case = TRUE)){
      save_pdf=paste0(save_pdf,'.pdf')
    }
  }
  if(is.character(cluster_use)){
    cluster_use=object@gene.info[,cluster_use,drop=F]
  }
  subexpr = slot(object, exp_use)
  if(dim(subexpr)[1]==0){ #if input is smoothed.exp but smoothed version of expression has not been calculated, calculate it first
    ## calculate smoothed expression requires pseudotime
    pt=object@pseudotime
    sd=diff(range(pt))/100
    subexpr = gauss.blur.matrix(object@exp.data, pt=pt, sd=sd, return_nn = F)
  }
  #n_memb=table(cluster_use)
  #n_clus=sum(n_memb>1)
  clus.list=split(rownames(cluster_use), cluster_use)
  plot.list=list()
  #turn the wide format expression matrix into a long format dataframe, and add cluster id as a new column (module)
  all_df = data.frame(subexpr) %>%
    mutate(gene_id = row.names(.)) %>%
    tidyr::pivot_longer(-gene_id) %>%
    mutate(module = cluster_use[gene_id,])
  all_df$name=as.numeric(gsub("X","",all_df$name))
  #all_df has columns gene_id, name (pseudotime), value (exp values), module (cluster id)
  
  #prepare annotation tables for functional annotation enrichment calculations
  anno.list=list()
  cols=c("gene","id")
  if(use_description){
    cols=c("gene","id","description")
  }
  for(anno in names(object@annotations)){
    anno.list[[anno]]=object@annotations[[anno]][,cols]
  }
  db_tbl= do.call(rbind,anno.list)
  if(use_description){
    id_desc=unique(db_tbl[,c("id","description")])
    rownames(id_desc)=id_desc$id
  }
  genexanno = df2matrix(db_tbl, rows="gene", cols="id")
  setnames=c("Top1st","Top2nd","Top3rd")
  for(cls in names(clus.list)){
    if(length(clus.list[[cls]])>1){
      cls_df=all_df[all_df$module==cls,]
      cls.genes=clus.list[[cls]]
      p <- ggplot(cls_df, aes(x = name, y = value, group=gene_id)) + geom_line(alpha=0.25) +
        labs(title=paste0('cluster',cls), x = "Time",  y = "Expression") + theme_minimal()
      #plot.list[[paste0(cls,".exp")]]=p
      ## get the 3 top enriched functional terms
      fisher_df=fisher.exact.enrichment(cls.genes, gXm=genexanno)
      nrows=min(3,dim(fisher_df)[1])
      if(nrows==0){
        v <- ggplot()+theme_void()
        anno <- ggplot()+theme_void()
      }else{
        fisher_df=fisher_df[1:nrows,]
        venn.list=lapply(1:nrows, function(x){unlist(strsplit(fisher_df[x,"genes.in.selection"], ", "))})
        names(venn.list)=setnames[1:nrows]
        venn.list[["all"]]=cls.genes
        venn_plot <- function() {#define function for plotting venn diagram
          v <- Venn(venn.list)
          if(nrows==1){
            plot(v[,c(setnames[1],"all")], doWeights = TRUE)
          }else{
            sets_equal=F
            if(setequal(venn.list[[setnames[1]]],venn.list[[setnames[2]]])){
              sets_equal=T
            }else{
              if(nrows==3){
                if(setequal(venn.list[[setnames[1]]],venn.list[[setnames[3]]]) || setequal(venn.list[[setnames[3]]],venn.list[[setnames[2]]])){
                  sets_equal=T
                }
              }
            }
            if(sets_equal){
              plot(v[,setnames[1:nrows]], doWeights = F)
            }else{
              tryCatch({plot(v[,setnames[1:nrows]], doWeights = TRUE)}, error=function(e){plot(v[,setnames[1:nrows]], doWeights = F)})
            }
          }
        }
        vpg <- plot_to_grob(venn_plot)
        v <- qplot(c(1,10), c(1,10), geom="blank") + 
          annotation_custom(vpg, xmin = 0, xmax = 10, ymin = 0, ymax = 10) + 
          theme_void() + coord_fixed(ratio = 1)
        if(use_description){
          text_add = paste(paste0(setnames[1:nrows], " enriched (p.adj=", format(fisher_df[1:nrows,'p.adj'], digits=2), "):\n ",
                                  id_desc[rownames(fisher_df)[1:nrows],'description']),
                           collapse = "\n")
        }else{
          text_add = paste(paste0(setnames[1:nrows], " enriched (p.adj=", format(fisher_df[1:nrows,'p.adj'], digits=2), "):\n ",
                                  rownames(fisher_df)[1:nrows]),
                           collapse = "\n")
        }
        anno <- qplot(c(1,20), c(1,20), geom="blank") + theme_void() +
          annotate("text", x=0, y=20, label = text_add, hjust = 0, vjust = "top") + theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
      }
      plot.list[[paste0("cluster",cls)]] <- (p + v) / anno
    }
  }
  pp=wrap_plots(plot.list, ncol = ncols)
  if(show_plot){
    print(pp)
  }else if(!is.character(save_pdf) && !return_plots){
    print(pp)
  }
  if(is.character(save_pdf)){
    ggsave(save_pdf, pp, width = save_width, height = save_height)
    ggsave("../test.pdf", pp, width = 20, height = 40)
  }
  if(return_plots){
    return(plot.list)
  }
}

#' Get Auto Module from Clusters
#'
#' This function generates modules from clusters, optionally naming them by enriched annotations.
#'
#' @param object An object containing clustering results and other relevant data.
#' @param cluster_use The name of the cluster to use.
#' @param name_by_enriched Boolean. If TRUE, name modules by enriched annotations.
#' @param use_description Boolean. If TRUE, use descriptions for enrichment analysis. Default is TRUE.
#'
#' @return A list of modules.
#'
#' @export
get_auto_module <- function(object, cluster_use, name_by_enriched=T, use_description=T){
  cluster_use=object@gene.info[,cluster_use,drop=F]
  clus.list=split(rownames(cluster_use), cluster_use)
  n_memb=table(cluster_use)
  clus.singl=clus.list[names(n_memb)[n_memb==1]]
  clus.list=clus.list[setdiff(names(clus.list),names(n_memb)[n_memb==1])]
  clus.list[["singlets"]]=as.character(unlist(clus.singl))
  if(!name_by_enriched){
    return(clus.list)
  }
  #prepare annotation tables for functional annotation enrichment calculations
  anno.list=list()
  cols=c("gene","id")
  if(use_description){
    cols=c("gene","id","description")
  }
  for(anno in names(object@annotations)){
    anno.list[[anno]]=object@annotations[[anno]][,cols]
  }
  db_tbl= do.call(rbind,anno.list)
  if(use_description){
    id_desc=unique(db_tbl[,c("id","description")])
    rownames(id_desc)=id_desc$id
  }
  genexanno = df2matrix(db_tbl, rows="gene", cols="id")
  new_name=c()
  for(cls in names(clus.list)){
    if(cls!="singlets"){
      cls.genes=clus.list[[cls]]
      fisher_df=fisher.exact.enrichment(cls.genes, gXm=genexanno)
      if(use_description){
        name_add = id_desc[rownames(fisher_df)[1],'description']
      }else{
        name_add = rownames(fisher_df)[1]
      }
      if(name_add%in%new_name){
        i=2
        name_add=paste0(name_add," (",i,")")
        while(name_add%in%new_name){
          i=i+1
          name_add=paste0(name_add," (",i,")")
        }
      }
      new_name=c(new_name, name_add)
    }else{
      new_name=c(new_name,"singlets")
    }
  }
  names(clus.list)=new_name
  return(clus.list)
}

#' Save Modules to TSV
#'
#' This function saves modules to a TSV file.
#'
#' @param auto_module A list of modules to save.
#' @param save_path The path to save the TSV file.
#'
#' @return None
#'
#' @export
module2tsv <- function(auto_module, save_path, return_df=F){
  save_module=lapply(auto_module, paste, collapse=", ")
  names(save_module)=names(auto_module)
  save_df=as.data.frame(t(data.frame(save_module)))
  rownames(save_df)=names(save_module)
  colnames(save_df)=c("member")
  save_df$name=rownames(save_df)
  save_df=save_df[,c('name','member')]
  if(!is.null(save_path)){
    write.table(save_df, save_path, quote = F, row.names = F, sep="\t")
  }
  if(return_df){
    return(save_df)
  }
}
