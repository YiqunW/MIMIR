#' Scale data by max and min values calculated with a top and bottom fraction of cells (to confer robustness to outliers)
#' @param exp_data A matrix representing gene expression data, with genes as rows and cells as columns.
#' @param max.scl The proportion of the total number of cells used for calculating the upper expression threshold of a gene. Default is 0.01.
#' @param min.sub The proportion of the total number of cells used for calculating the lower expression threshold of a gene. Default is 0.01.
#' @param unlog A logical value indicating whether to unlog the data before processing. If the data is in log scale, then 'unlog=T' should be used.
#'              If TRUE, the data is unlogged (using a base-2 logarithm) at the start and re-logged at the end.
#'
#' @return A matrix of the same dimensions as `exp_data`, with the expression data scaled.
#'
#' @export
scale_trim <- function(exp_data, max.scl=0.01, min.sub=0.01,unlog=T){
  num.cells=dim(exp_data)[2]
  num.max=max(c(round(num.cells*max.scl),1))
  num.min=max(c(round(num.cells*min.sub),1))
  if(unlog){
    exp_data=2^exp_data-1
  }
  exp_sort=apply(exp_data,1,function(x) sort(x,decreasing = T)) #this produce a transposed matrix of cells by genes
  exp.min=apply(exp_sort[(num.cells-num.min+1):num.cells,,drop=F],2,mean)
  exp_sort=sweep(exp_sort, 2, exp.min, "-")
  exp.max=apply(exp_sort[1:num.max,,drop=F],2,mean)
  exp_scl=sweep(exp_data, 1, exp.min, "-")
  exp_scl=sweep(exp_scl, 1, exp.max, "/")*(mean(exp.max))
  exp_scl=exp_scl*(exp_scl>0)
  if(unlog){
    exp_scl=log2(exp_scl+1)
  }
  return(exp_scl)
}

#' Calculate Cosine Distance Between Rows of a Data Matrix
#'
#' The cosine distance is calculated as one minus the cosine of the angle between them.
#'
#' @param x A numeric matrix with rows representing different items (e.g. genes) and 
#'          columns representing features (e.g. pseudotimes or cells).
#'
#' @return A square, symmetric matrix where each element [i, j] represents the cosine distance 
#'         between the i-th and j-th rows of the input matrix `x`. Values range from 0 (identical)
#'         to 1 (orthogonal).
#'
#' @examples
#' data_matrix <- matrix(rnorm(20), nrow = 5)
#' distance_matrix <- cosineDist(data_matrix)
#' print(distance_matrix)
#'
#' @export
cosineDist <- function(x){
  1-(x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

#' Generate a Similarity Matrix with Gaussian Blur
#'
#' This function creates a similarity matrix for a given set of values (e.g. pseudotime - numeric values reflecting developmental progression).
#' Each value is convolved with a Gaussian function to apply a blur effect,
#' resulting in a matrix where similarities are smoothed across nearby time points.
#'
#' @param pt Numeric vector for which the similarity matrix will be generated.
#'           Values are automatically sorted in ascending order.
#' @param sd Standard deviation for the Gaussian blur. Default is 0.01.
#'           Determines the width of the Gaussian distribution.
#' @param cutoff Positive number. Cutoff for the Gaussian function. Default is 4 standard deviations.
#'           Values outside the cutoff range have no effect.
#'
#' @return A square matrix of size length(pt) with blurred similarities.
#'         Each element represents the similarity between two points,
#'         smoothed by a Gaussian blur.
#' @examples
#' pts <- c(1, 2, 3, 4, 5)
#' sim_matrix <- simM(pts)
#' print(sim_matrix)
#' 
#' @export
simM <- function(pt, sd=0.01, cutoff=4){
  pt=sort(as.numeric(pt), decreasing = F)
  s=matrix(0,nrow=length(pt),ncol=length(pt))
  diag(s)=1
  s.blur=s*0
  for(i in 1:dim(s)[2]){
    pt_center=pt[i]
    gaus.conv=dnorm(pt,mean=pt_center,sd=sd)
    gaus.conv[which(pt<(pt_center-cutoff*sd))]=0
    gaus.conv[which(pt>(pt_center+cutoff*sd))]=0
    s.blur[,i]=gaus.conv
  }
  s.blur=s.blur/s.blur[1,1]
  colnames(s.blur)=pt
  rownames(s.blur)=pt
  return(s.blur)
}

#' Soft Cosine Distance Matrix Calculation
#'
#' This function calculates a distance matrix for a given dataset, using a soft cosine similarity measure.
#' It applies an eigenvalue decomposition on a provided similarity matrix 's' and transforms the input data 'x'
#' accordingly. The function then computes the cosine similarity between pairs of rows in the transformed data.
#'
#' @param x A numeric matrix where rows represent items (e.g. genes) and columns represent features (e.g. pseudotimes).
#'          The function computes the distance between each pair of rows.
#' @param s A square, symmetric matrix representing feature similarities (calculated by simM()).
#'          The columns of 's' should match that of 'x'.
#'
#' @return A square matrix of distances between each pair of rows in 'x'.
#'         Distances are calculated using a soft cosine measure, which accounts for feature similarities.
#'         Values range from 0 (identical) to 1 (maximally distant).
#'
#' @export
soft.cos.dist <- function(x, s){
  simM=matrix(1,nrow=dim(x)[1],ncol=dim(x)[1])
  rownames(simM)=rownames(x)
  colnames(simM)=rownames(x)
  #x=x[,colnames(s)]
  eigen_s=eigen(s)
  x_n=t(eigen_s$vectors) %*% t(x)
  self.cos=apply(x_n,2,function(i) sum(i*i*eigen_s$values))
  num_gen=dim(x)[1]
  for(i in 1:(num_gen-1)){
    vi=x_n[,i]
    for(j in (i+1):(num_gen)){
      vj=x_n[,j]
      de=sqrt(self.cos[i]*self.cos[j])
      simM[i,j]=sum(vi*vj*eigen_s$values)/de
      simM[j,i]=simM[i,j]
    }
  }
  return(1-simM)
}

gauss.blur.matrix <- function(exp_data, pt=NULL, sd=0.1, unlog=T, cutoff=3, return_nn=F){
  #exp_data is a gene by cell matrix with pseudotime as colnames
  if(is.null(pt)){
    pt=as.numeric(colnames(exp_data))
  }
  expv=exp_data
  if(unlog){
    expv=2^expv-1
  }
  expv.blur=expv*0
  num.points=c()
  for(i in 1:dim(expv)[2]){
    pt_center=pt[i]
    gaus.conv=dnorm(pt,mean=pt_center,sd=sd)
    gaus.conv[which(pt<(pt_center-cutoff*sd))]=0
    gaus.conv[which(pt>(pt_center+cutoff*sd))]=0
    gaus.conv=gaus.conv/sum(gaus.conv)
    if(return_nn){
      num.points=c(num.points,sum(gaus.conv>0))
    }
    expv.blur[,i]=apply(expv*matrix(gaus.conv,nrow=dim(exp_data)[1],ncol=length(gaus.conv),byrow=T), 1, sum)
  }
  if(unlog){
    expv.blur=log2(expv.blur+1)
  }
  if(return_nn){
    return(num.points)
  }else{
    return(expv.blur)
  }
}

library(abind) # required for storing results in 3D arrays
library(philentropy) # required for JS divergence calculations
library(gplots) # required for plotting the gaussian blur matrix if plot_blur is TRUE.
#' Calculate the expression distances (dissimilarities) between genes in an expression matrix using several methods and store the results in a 3D data array.
#'
#' @param x A numeric matrix with rows representing different items (e.g. genes) and 
#'          columns representing features (e.g. pseudotimes or cells).
#' @param methods A vector of strings specifying the methods for calculating expression distance. 
#'                Available methods are "cosine", "euclidean", "canberra", "JS" (Jensen-Shannon distance), and "soft_cosine".
#' @param smooth Boolean. Whether to smooth gene expression over pseudotime (through gaussian blur) to reduce noise. 
#'               If TRUE, pt (pseudotime of cells) must be provided. When calculating soft_cosine distance, 
#'               this step will be skipped as soft_cosine inherently blurs the data over pt. 
#' @param pt If method is "soft_cosine" or if smooth=TRUE: numeric vector of pseudotime corresponding to the columns of x. 
#' @param sd Standard deviation for the Gaussian blur used in soft_cosine distance calculation and/or for smoothing the data. 
#'           Default is 0.01. 
#' @param cutoff Used in soft_cosine distance calculation and/or when smooth=TRUE. Positive number. 
#'               Cutoff for the Gaussian distribution used for blurring gene expression along pseudotime. 
#'               Default is 4 standard deviations. Values outside the cutoff range have no effect.
#' @param plot_blur Boolean. Used when "soft_cosine" is in methods and/or if smooth=TRUE. 
#'                  Whether to visualized the similarity matrix used in soft_cosine calculation 
#'                  and/or smoothing the data. The similarity matrix is used to modify a gene's expression 
#'                  in a cell based on its expression in cells with similar pseudotimes. 
#'               
#' @return A 3D data array of gene by gene by method. Each element [i, j, m] represents the expression distance
#'         between gene i and gene j using method m. Values range from 0 (identical) to 1 (orthogonal).
#'
#' @export
exp_dis <- function(x, methods=c("cosine", "euclidean", "canberra", "JS", "soft_cosine"), smooth=F, 
                    pt=NULL, sd=0.01, cutoff=4, plot_blur=F, verbose=T){
  ## if "soft_cosine" is in methods and/or if smooth=T, replace column names of x with pt.
  if("soft_cosine"%in%methods || smooth){
    if(is.null(pt) || length(pt)!=dim(x)[2]){
      print("Need to supply pt that matches the columns of the input data matrix.")
      return(NA)
    }else{
      colnames(x)=pt
      pt_order=order(as.numeric(pt))
      x=x[,pt_order]
      pt=pt[pt_order]
    }
    if(plot_blur || "soft_cosine"%in%methods){
      if(verbose){
        print("Calculating gaussian blur matrix for soft_cosine and/or expression smoothing...")
      }
      s.blur=simM(pt = as.numeric(pt),sd = sd,cutoff = cutoff)
      if(plot_blur){
        heatmap.2(s.blur, Rowv = F, Colv = F, dendrogram = 'none', main='gaussian blur matrix',
                  sepwidth = c(0,0),trace='none',density.info='none',symm=T)
      }
    }
  }
  
  distMs=list() # creat a list to store distance matrices
  
  ## calculate soft_cosine if it is in the methods
  if("soft_cosine"%in%methods){
    print("Calculating soft cosine distance...")
    distMs[["soft_cosine_dist"]]=soft.cos.dist(x, s.blur)
  }
  
  ## calculate a gaussian blurred version of expression data if smooth=T
  if(smooth){
    print("Smoothing expression data...")
    x=gauss.blur.matrix(x, sd=sd, return_nn = F)
  }
  
  if("cosine"%in%methods){
    print("Calculating cosine distance...")
    distMs[["cosine_dist"]]=cosineDist(x)
  }
  
  if("euclidean"%in%methods){
    print("Calculating euclidean distance...")
    distMs[["euclidean_dist"]]=as.matrix(dist(x,method = 'euclidean'))
  }
   
  if("canberra"%in%methods){
    print("Calculating canberra distance...")
    distMs[["canberra_dist"]]=as.matrix(dist(x,method = 'canberra'))
  }  
  
  if("JS"%in%methods){
    print("Calculating Jensen-Shannon distance...")
    jsdiv=JSD(x, est.prob = 'empirical')
    rownames(jsdiv)=rownames(x)
    colnames(jsdiv)=rownames(x)
    distMs[['JS_dist']]=sqrt(jsdiv)
  }
  
  ## return results in form of a 3D array
  return(do.call(abind,c(distMs,along=3)))
}

#' Convert expression distances between genes to expression similarities
#'
#' @param x A distance matrix (2D) of non-negative values with rows and columns representing genes and elements representing their gene expression distances.
#' @param max.score A numeric value <=1. The max value of the final scaled similarity.
#'               
#' @return A matrix of the same dimensions to the input matrix
#'
#' @export
dist_to_sim <- function(x, max.score=0.95){
  expM=x*(x>0)
  expM=1-(expM/max(expM))
  expM=expM*max.score
  return(expM)
}

#' Convert expression distances between genes to expression similarities
#'
#' @param x A 3D array of expression distance (gene x gene x method). Out put from the exp_dis function.
#' @param max.score A numeric value <=1. The max value of the final scaled similarity.
#'               
#' @return A 3D array of the same dimensions to the input.
#'
#' @export
all_dist_to_sim <- function(x, max.score=0.95){
  exp.sim=x*0
  for (exp.dist in dimnames(x)[[3]]){
    exp.sim[,,exp.dist]=dist_to_sim(x[,,exp.dist],max.score=max.score)
  }
  dim3names=gsub("dist", "sim", dimnames(x)[[3]])
  dimnames(exp.sim)[[3]]=dim3names
  return(exp.sim)
}