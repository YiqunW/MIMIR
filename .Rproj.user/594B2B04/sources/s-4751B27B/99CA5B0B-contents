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
#' This function creates a similarity matrix for a given set of points.
#' Each point is convolved with a Gaussian function to apply a blur effect,
#' resulting in a matrix where similarities are smoothed across nearby points.
#'
#' @param pt Numeric vector of points for which the similarity matrix will be generated.
#'           Points are automatically sorted in ascending order.
#' @param sd Standard deviation for the Gaussian blur. Default is 0.01.
#'           Determines the width of the Gaussian distribution.
#' @param cutoff Positive number. Cutoff for the Gaussian function. Default is 4 standard deviations.
#'           Points outside the cutoff range have no effect.
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