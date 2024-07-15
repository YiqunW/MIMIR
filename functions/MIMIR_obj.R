#' MIMIR Class
#'
#' A class to represent MIMIR objects with various attributes related to gene expression data.
#'
#' @slot name A character string representing the name of the project.
#' @slot exp.data A matrix representing the expression data, with genes as rownames and pseudotimes as colnames. Alternatively, pseudotimes can be provided directly into the pseudotime slot.
#' @slot smoothed.exp A matrix representing the expression data smoothed (gaussian blurred) along time.
#' @slot genes A vector of gene names.
#' @slot pseudotime A vector of pseudotime values.
#' @slot exp.dis An array representing expression distances (gene x gene x distance metric).
#' @slot exp.sim An array representing expression similarities (gene x gene x similarity metric).
#' @slot func.sim An array representing functional similarities (gene x gene x similarity metric).
#' @slot integrated.sim An array representing integrated similarities (gene x gene x similarity metric).
#' @slot clusters A list of clustering results.
#' @slot cluster_metrics A data frame of metrics and statistics of clustering results using different parameters. Columns are quality metrics and rows are different clustering results.
#' @slot annotations A list of annotations from different annotation databases.
#' @slot gene.info A data frame containing gene information, such as expression onset times, clustering memberships, etc.
#'
#' @export
MIMIR <- methods::setClass("MIMIR", 
                           slots = c(name="character",
                                     exp.data = "matrix",
                                     smoothed.exp = "matrix",
                                     genes = "vector",
                                     pseudotime = "vector",
                                     exp.dis = "array",
                                     exp.sim = "array",
                                     func.sim = "array",
                                     integrated.sim = "array",
                                     clusters = "list",
                                     cluster_metrics="data.frame",
                                     annotations = "list",
                                     gene.info = "data.frame"
                                     ))

#' Show Method for MIMIR Class
#'
#' Displays a summary of the MIMIR object.
#'
#' @param object An object of class MIMIR.
#'
#' @export
setMethod(f = "show", signature = "MIMIR",
          definition = function(object) {
            if(length(object@pseudotime)==0){
              pt="with no temporal information:"
            }else{
              pt="with temporal information:"
            }
            cat("MIMIR object for", object@name, pt,
                nrow(object@exp.data), "genes x",
                ncol(object@exp.data), "cells.\n")
            invisible(NULL)
            }
          )

#' Create MIMIR Object
#'
#' Creates a new MIMIR object with the specified parameters.
#'
#' @param name A character string representing the name of the project.
#' @param exp.data A matrix representing the expression data, with genes as rownames and pseudotimes as colnames. Alternatively, pseudotimes can be provided directly into the pseudotime slot.
#' @param genes A vector of gene names. Defaults to NULL.
#' @param pseudotime A vector of pseudotime values matching the column order in the exp.data matrix. Defaults to NULL.
#'
#' @return An object of class MIMIR.
#'
#' @export
createMIMIR <- function(name, exp.data, genes = NULL, pseudotime = NULL) {
  if(is.null(genes)){
    genes=rownames(exp.data)
  }else{
    genes=intersect(genes,rownames(exp.data))
  }
  if(is.null(pseudotime)){
    ## test if the expression data column name contain numeric values that can be used as pseudotime
    if(sum(is.na(as.numeric(colnames(exp.data))))==0){
      pseudotime=as.numeric(colnames(exp.data))
    }
  }else{
    if(sum(is.na(as.numeric(pseudotime)))>0){
      message("pseeudotime argument unused. Input contains non-numeric values.")
      pseudotime=NULL
    }else{
      if(length(pseudotime)!=dim(exp.data)[2]){
        message("pseeudotime argument unused. Input doesn't match number of columns in exp.data.")
        pseudotime=NULL
      }else{
        pseudotime=as.numeric(pseudotime)
      }
    }
  }
  object <- methods::new("MIMIR", name=name, exp.data=exp.data, genes=genes)
  if(!is.null(pseudotime)){
    object@pseudotime <- pseudotime
    sums=apply(exp.data,1,sum)
    ave.time=apply(pseudotime*t(exp.data),2,sum)/sums
    object@gene.info <- data.frame("weighted.time"=ave.time, row.names = rownames(exp.data))
  }
  return(object)
} 

#' Anno Class
#'
#' A class to represent annotation objects with various attributes related to gene annotations.
#'
#' @slot db A character string representing the database used.
#' @slot gene.anno A data frame of gene annotations, with columns "gene", "id", "description"..
#' @slot genes A vector of genes for which functional similarities are to be calculated.
#' @slot bg.genes A vector of background genes used to adjust calculated similarities.
#' @slot IC A vector of information content values of annotation terms.
#' @slot cat.IC A list storing the information content of annotation terms for each GO category.
#' @slot hierarchy.df A data frame storing the parent-child relations between annotation terms.
#' @slot hierarchy.l A list representing the hierarchy of annotation terms.
#' @slot descrip A data frame of descriptions.
#' @slot similarity A matrix of gene x gene annotation similarities.
#' @slot cat.similarity A list of similarity matrices for GO categories.
#' @slot prior A numeric value of prior similarity value (mean similarity across the union of bg.genes and genes).
#' @slot cat.prior A list of prior values for GO categories.
#'
#' @export
Anno <- methods::setClass("Anno", 
                           slots = c(db="character",
                                     gene.anno="data.frame",
                                     genes = "vector",
                                     bg.genes="vector",
                                     IC = "vector",
                                     cat.IC = "list",
                                     hierarchy.df = "data.frame",
                                     hierarchy.l = "list",
                                     descrip = "data.frame",
                                     similarity="matrix",
                                     cat.similarity="list",
                                     prior="numeric",
                                     cat.prior="list"
                           ))

#' Show Method for Anno Class
#'
#' Displays a summary of the Anno object.
#'
#' @param object An object of class Anno.
#'
#' @export
setMethod(f = "show", signature = "Anno",
          definition = function(object) {
            cat("Anno object. Database used: ", object@db, ".\n")
            invisible(NULL)
          }
)

#' Create Anno Object
#'
#' Creates a new Anno object with the specified parameters.
#'
#' @param db A character string representing the database used.
#' @param gene.anno A data frame of gene annotations, with columns "gene", "id", "description"..
#' @param gene.col The column name or index for genes in the gene.anno data frame.
#' @param id.col The column name or index for annotation IDs in the gene.anno data frame.
#' @param desp.col The column name or index for annotation descriptions in the gene.anno data frame. Defaults to NULL.
#' @param cat.col The column name or index for categories of annotations in the gene.anno data frame. Defaults to NULL.
#' @param genes A vector of genes for which functional similarities are to be calculated. Defaults to NULL.
#' @param bg.genes A vector of background genes. Defaults to NULL.
#' @param hierarchy.df A data frame storing the parent-child relations between annotation terms. Defaults to NULL.
#' @param parent.col The column name or index for parent nodes in the hierarchy.df data frame. Defaults to 1.
#' @param child.col The column name or index for child nodes in the hierarchy.df data frame. Defaults to 2.
#' @param hierarchy.l A list representing the hierarchy of annotations. Defaults to NULL.
#'
#' @return An object of class Anno.
#'
#' @export
createAnno <- function(db, gene.anno, gene.col='gene', id.col='id', desp.col=NULL, cat.col=NULL, genes = NULL, bg.genes = NULL, 
                       hierarchy.df=NULL, parent.col=1, child.col=2, hierarchy.l=NULL){
  gene.anno.use=data.frame('gene'=gene.anno[,gene.col], 'id'=gene.anno[,id.col])
  if(!is.null(cat.col)){
    gene.anno.use$category=gene.anno[,cat.col]
  }
  if(!is.null(desp.col)){
    gene.anno.use$description=gene.anno[,desp.col]
  }
  if(!is.null(genes)){
    if(length(setdiff(genes, gene.anno.use$gene))>0){
      n=length(setdiff(genes, gene.anno.use$gene))
      message(paste0(n,' input genes are not found in the gene.anno table'))
      genes=intersect(genes, gene.anno.use$gene)
    }
  }
  if(!is.null(bg.genes)){
    if(length(setdiff(bg.genes, gene.anno.use$gene))>0){
      n=length(setdiff(bg.genes, gene.anno.use$gene))
      message(paste0(n,' input background genes are not found in the gene.anno table'))
      bg.genes=intersect(bg.genes, gene.anno.use$gene)
    }
  }
  if(!is.null(genes) && is.null(bg.genes)){
    all.genes=union(genes,bg.genes)
    if(length(setdiff(gene.anno.use$genes, all.genes))>0){
      message("Subsetting gene.anno table to remove genes not in the input gene or bg.genes...")
      gene.anno.use=gene.anno.use[which(gene.anno.use$gene%in%all.genes),]
    }
  }
  object <- methods::new("Anno", db=db, gene.anno=gene.anno.use)
  if(!is.null(genes)){
    object@genes = genes
  }
  if(!is.null(bg.genes)){
    object@bg.genes = bg.genes
  }
  if(!is.null(hierarchy.df)){
    object@hierarchy.df <- data.frame('parent'=hierarchy.df[,parent.col], 'child'=hierarchy.df[, child.col])
    object@hierarchy.l <- offspring.list(object@hierarchy.df)
  }
  if(!is.null(hierarchy.l)){
    if(is.null(hierarchy.df)){
      object@hierarchy.l = hierarchy.l
      p2c=stack(hierarchy.l)
      p2c=p2c[,c(2,1)]
      colnames(p2c)=c("parent","child")
      p2c=as.data.frame(as.matrix(p2c))
      object@hierarchy.df = p2c
    }else{
      message("Argument hierarchy.l unused. hierarchy.l calculated from the supplied hiearchy.df.")
    }
  }
  return(object)
} 