#' @title Localization procedure
#'
#' @description The localization procedure to detect change-points.
#'
#' @param data A numeric matrix of observations with each horizontal axis being time, and each column being the multivariate time series
#' @param intervals A numeric matrix of intervals with each row be a vector representing the interval
#' @param l A non-negative integer of order of polynomial (\code{l=0} means piecewise constant)
#' @param scaling A logical scalar representing whether to perform refinement for locally stationary data. Only useful when \code{l=0}
#' @param q A positive integer of norm
#'
#' @return The positions of estimated change-points
#'
#' @examples
#'set.seed(0)
#'data=rep(c(0,2,0),each=40)+rnorm(120)
#'d=rid(data,M=1000,tau="clustering")
#'cpt=localization(data,d$Good_Intervals)
#'print(cpt)
#'
#' @export

localization <- function(data,intervals,l=0,scaling=FALSE,q=1){
  stopifnot(
    is.numeric(data),
    is.vector(data)|is.matrix(data),
    is.null(intervals)|is.matrix(intervals),
    is.numeric(l)
  )
  l=as.integer(l)
  if(length(intervals)==0){return(NULL)}
  if(is.vector(data)){
    data=t(as.matrix(data))
  }
  others=vector("list", nrow(data))
  if(l>0){
    fs=Get_fs_l(data,intervals,l,q,IProj)$Pos
  }
  if(l==0 & scaling==F){
    fs=Get_pos_univariate(data,intervals,g_cusum,as.double(q),others)+1
  }
  if(l==0 & scaling==TRUE){
    nrw=nrow(data)
    ncl=ncol(data)
    lrv=list()
    for(j in 1:nrw){
      lrv[[j]]=Longrundim1_lconst(data[j,])
    }
    others=lrv
    fs=Get_pos_univariate(data,intervals,g_cusum_lrv,as.double(q),others)+1
  }
  return(fs)
}
