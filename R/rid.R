#' @title Random intervals distillation procedure
#'
#' @description Distilled intervals that cover change-points are constructed.
#'
#' @param data A numeric matrix of observations with each horizontal axis being time, and each column being the multivariate time series
#' @param M A positive integer of random intervals, used only when \code{intervals=NULL}
#' @param l A non-negative integer of order of polynomial (\code{l=0} means piecewise constant)
#' @param scaling A logical scalar representing whether to perform refinement for locally stationary data. Only useful when \code{l=0}
#' @param q A positive integer of norm
#' @param intervals A numeric matrix of intervals with each row be a vector representing the interval. If \code{intervals=NULL}, random intervals are sampled
#' @param tau A non-negative number representing the threshold of detection. If \code{tau="clustering"} (only useful when \code{l=0}), a clustering-based adaptive approach is applied. If \code{tau="ref"}, a method based on simulated reference values is applied
#' @param bw A parameter passed into function \code{\link{density}}. Only useful when \code{tau="clustering"}
#' @param adjust A parameter passed into function \code{\link{density}}. Only useful when \code{tau="clustering"}
#' @param k.max A positive integer representing the maximum value of clusters in threshold determining. Only useful when \code{tau="clustering"}
#' @param adj A positive number used to multiply onto the threshold \code{tau}, providing threshold adjustments for small sample size
#'
#' @return A list containing:
#' \item{Good_Intervals}{A numeric matrix with each row being an interval that covers a change-point}
#' \item{Threshold}{A positive number of the threshold}
#'
#' @seealso{\code{\link{localization}}}
#'
#' @docType package
#' @useDynLib rid, .registration = TRUE
#'
#' @examples
#'## An example for the univariate case
#'set.seed(0)
#'data=rep(c(0,2,0),each=40)+rnorm(120)
#'d=rid(data,M=1000,tau="clustering")
#'cpt=localization(data,d$Good_Intervals)
#'print(cpt)
#'
#'## An example for the multivariate case
#'set.seed(0)
#'data1=rep(c(0,2,0),each=40)+rt(120,8)
#'data2=rep(c(0,2,0),each=40)+rnorm(120)
#'data=rbind(data1,data2)
#'d=rid(data,M=1000,tau="clustering")
#'cpt=localization(data,d$Good_Intervals)
#'print(cpt)
#'
#'## An example for the piecewise polynomial case
#'set.seed(0)
#'n=300
#'cp=c(0,round(n/3),round(2*n/3),n)
#'mu=matrix(c(0.004,-0.1,0,-0.01,0.02,0,0.01,-0.04,0),nrow=3,byrow = TRUE)
#'mu1=mu[1,1]*(1:n)^2+mu[1,2]*(1:n)+mu[1,3]
#'for(j in 2:3){
#'  index=which((1:n)-cp[j]>0)
#'  tmp1=(1:n)-cp[j]
#'  tmp=mu[j,1]*tmp1^2+mu[j,2]*tmp1+mu[j,3]
#'  tmp[1:(index[1]-1)]=0
#'  mu1=mu1+tmp
#'}
#'data=mu1+runif(n,-6,6)
#'plot(data,type="l")
#'d=rid(data,M=500,tau="ref",l=2)
#'cpt=localization(data,d$Good_Intervals,l=2)
#'print(cpt)
#'
#'## An example for refinement in the locally stationary time series
#'set.seed(0)
#'n=1000
#'cp=c(0,round(n/4),round(3*n/4),n)
#'epsilon=rnorm(500+n,0,1)
#'ei=rep(0,500+n)
#'for(j in 2:(500+n)){ei[j]=0.5*ei[j-1]+epsilon[j]}
#'x=ei[(501):(500+n)]
#'lrv=purrr::map_dbl(1:n,function(j){sqrt(max(1,2000*j/n))})
#'x=x*lrv
#'x=x+c(rep(0,round(n/4)),rep(20,round(n/2)),rep(-20,n-round(n/4)-round(n/2)))
#'data=x
#'plot(data,type="l")
#'d=rid(data,M=1000,scaling = TRUE,tau="clustering")
#'cpt=localization(data,d$Good_Intervals)
#'print(cpt)
#'
#' @export

rid <- function(data,M=1000,l=0,scaling=FALSE,q=1,
                intervals=NULL,tau=c("clustering","ref"),
                bw="nrd0",adjust=0.5,k.max=4,adj=1.3){
  stopifnot(
    is.numeric(data),
    is.vector(data)|is.matrix(data),
    is.numeric(tau)|tau%in%c("clustering","ref"),
    is.null(intervals)|is.matrix(intervals),
    is.numeric(l)
  )
  l=as.integer(l)
  if(M<100){stop("M is too small!")}
  if(is.vector(data)){
    data=t(as.matrix(data))
  }
  if(ncol(data)<10){stop("Sample size is too small!")}

  if(is.null(intervals)){
    intervals=Random_Intervals(ncol(data),M)
    M=nrow(intervals)
  }else{
    M=nrow(intervals)
  }
  others=list()
  if(l==0 & scaling==F){
    fs=Get_fs(data,intervals,f_cusum,as.double(q),list())
  }else if(l==0 & scaling==TRUE){
    nrw=nrow(data)
    ncl=ncol(data)
    lrv=data
    for(j in 1:nrw){
      lrv[j,]=Longrundim1_lconst(data[j,])
    }
    if(min(lrv)<=0){lrv=lrv-min(lrv)+min(lrv[lrv>0])}
    others=list(lrv)
    fs=Get_fs(data,intervals,f_cusum_lrv,as.double(q),others)
  }
  if(l>0){
    if(!is.numeric(tau)){tau="ref"}
    fs_l=pp(data,M=M,q=q,l=l,intervals=intervals,tau=tau)
    fs=fs_l$Statistic_values
  }

  if(l>0 & tau=="ref"){
    tau=fs_l$Threshold*log(log(ncol(data)))
  }
  if(l==0 & tau=="clustering"){
    tau_can=Clustering_tau(fs,bw=bw,k.max=k.max,adjust=adjust)
    tau_can=tau_can[[1]]
    tau_sim=Sim_tau(data,scaling,q,others)*log(log(ncol(data)))
    ind=(tau_can>tau_sim/10)&(tau_can<tau_sim*10)
    if(sum(ind)==0){tau=tau_sim
    }else{
      tau_can=tau_can[ind]
      tau=tau_can[which.min(abs(tau_can-tau_sim))]
    }
  }
  if(l==0 & tau=="ref"){
    tau=Sim_tau(data,scaling,q,others)*log(log(ncol(data)))
  }
  tau=tau*adj
  inn=which(fs>tau)
  if(length(inn)==0){
    good_intervals=NULL
  }else{
    good_intervals=Get_good_intervals(intervals[inn,,drop=F])
    colnames(good_intervals)=c("Left","Right")
  }
  return(list(Good_Intervals=good_intervals,
              Threshold=tau))
}
