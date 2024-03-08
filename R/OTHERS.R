#'
Get_pos_univariate<-function(data,intervals,g,q,others){
  M=nrow(intervals)
  p=nrow(data)
  result=rep(0,M)
  for(j in 1:M){
    s=intervals[j,1]
    e=intervals[j,2]
    cusum=rep(0,p*(e-s+1))
    dim(cusum)=c(p,e-s+1)
    for (k in 1:p){
      cusum[k,]=g(data[k,],as.integer(s),as.integer(e),others[[k]])
    }
    result[j]=which.max(apply(cusum,2,function(t)sum(abs(t)^q)))[1]-1+s
  }
  result
}


#'
wapply1 <- function(x, width, by = 1, FUN = mean, ...){
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width

  lenX <- length(x)
  SEQ1 <- seq(1, lenX - width + 1, by = by)
  SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))

  OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
  OUT <- base::simplify2array(OUT, higher = TRUE)
  return(OUT)
}

#'
fastlocal_constant=function(yy,bb,T=1/length(yy),endj=nn){
  nn=length(yy)

  a1=rep(0,nn)
  a2=rep(0,nn)
  a3=rep(0,nn)
  a4=rep(0,nn)
  a5=rep(0,nn)
  a6=rep(0,nn)

  b1=rep(0,nn)
  b2=rep(0,nn)
  b3=rep(0,nn)
  b4=rep(0,nn)
  b5=rep(0,nn)
  b6=rep(0,nn)

  TT=rep(0,nn)
  Tstart=nn*T

  TT[1]=T

  a1[1]=sum(yy*((1:nn)>=1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a2[1]=sum((1:nn)*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a3[1]=sum((1:nn)^2*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  a4[1]=sum((1:nn)^3*yy*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )


  b1[1]=sum(((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b2[1]=sum((1:nn)*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b3[1]=sum((1:nn)^2*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b4[1]=sum((1:nn)^3*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  b5[1]=sum((1:nn)^4*((1:nn)>=  1)*((1:nn)<=floor(Tstart+nn*bb)) )
  Lend=-1*floor((nn*bb)-Tstart)
  Rend=floor(Tstart+(nn*bb))+1
  Endj=endj
  for(i in 2: Endj)
  {lend1=(Lend>=1)
  rend1= (Rend<=Endj)
  lend=Lend*lend1+(Lend<1)
  rend=Rend*rend1+(Rend>Endj)
  a1[i]=a1[i-1]-yy[lend]*(1)*lend1+yy[rend]*(1)*rend1
  a2[i]=a2[i-1]-yy[lend]*(lend)*lend1+yy[rend]*(rend)*rend1
  a3[i]=a3[i-1]-yy[lend]*(lend^2)*lend1+yy[rend]*(rend^2)*rend1
  a4[i]=a4[i-1]-yy[lend]*(lend^3)*lend1+yy[rend]*(rend^3)*rend1

  b1[i]=b1[i-1]-(1)*lend1+(1)*rend1
  b2[i]=b2[i-1]-(lend)*lend1+(rend)*rend1
  b3[i]=b3[i-1]-(lend^2)*lend1+(rend^2)*rend1
  b4[i]=b4[i-1]-(lend^3)*lend1+(rend^3)*rend1
  b5[i]=b5[i-1]-(lend^4)*lend1+(rend^4)*rend1

  T=T+1/nn
  Rend=Rend+1
  Lend=Lend+1
  TT[i]=T
  }
  Y_1=a1
  Y_2=bb^(-1)*(a2/nn-TT*a1)
  Y_3=bb^(-2)*(a3/nn^2-2*a2*TT/nn+TT^2*a1)
  Y_4=bb^(-3)*(a4/nn^3-3*a3*TT/nn^2+3*a2*TT^2/nn-TT^3*a1)

  R_1=b1
  R_2=bb^(-1)*(b2/nn-TT*b1)
  R_3=bb^(-2)*(b3/nn^2-2*b2*TT/nn+TT^2*b1)
  R_4=bb^(-3)*(b4/nn^3-3*b3*TT/nn^2+3*b2*TT^2/nn-TT^3*b1)
  R_5=bb^(-4)*(b5/nn^4-4*b4*TT/nn^3+6*b3*TT^2/nn^2-4*b2/nn*TT^3+TT^4*b1)
K1=3/4*(R_1-R_3)/nn/bb
K2=3/4*(R_2-R_4)/nn/bb

  K3=3/4*(R_3-R_5)/nn/bb
Y1=3/4*(Y_1-Y_3)/nn/bb

  Y2=3/4*(Y_2-Y_4)/nn/bb

  coeff1=Y1/(K1)
  return(list(mean=coeff1[1:endj],time=TT[1:endj]))

}

#'
Longrundim1_lconst=function(xx,m=1,tau=NA,nn2=length(xx)){
  if(is.na(nn2)){nn2=length(xx)}
  if(is.na(tau)){tau=0.5/length(xx)^{1/7}}
  nn=length(xx)
  rolmean=wapply1(x=xx,width=m,by=1)
  Diff=rolmean[1:(nn-2*m)]-rolmean[(m+1):(nn-m)]
  tt=(1:nn2)/nn2
  Q=rep(0,nn)

  Q[(m+1):(nn-m)]=Diff^2*m/2
  Longvar=fastlocal_constant(yy=Q,bb=tau,T=1/nn,endj=nn)$mean

  ss=m
  Longvar[1:ss]=Longvar[(ss+1)]
  Longvar[(nn-ss+1):nn]=Longvar[(nn-ss)]
  return(Longvar)
}


#' @importFrom CircularSilhouette fast.sil
Clustering_tau<-function(fs,bw,k.max,adjust){
  x=fs[fs>0]
  si=rep(0,k.max-1)
  names(si)=1
  density_list=densitypeak_cluster_tmp(x,bw,adjust)
  for(k in 2:k.max){
    res=densitypeak_cluster(density_list,k)
    x1=res[[1]]
    flag=res[[2]]
    si[k-1]=CircularSilhouette::fast.sil(x1,flag)
    names(si)[k-1]=res[[3]]
  }
  kopt=as.numeric(names(which.max(si)))
  res=densitypeak_cluster(density_list,kopt)
  x1=res[[1]]
  res=res[[2]]
  tau_can=sort(as.vector(unlist(tapply(x1,res,range))))
  tau_can=sort(tau_can)
  tau_can=tau_can[2:(length(tau_can)-1)]
  return(list(tau_can))
}

#' @importFrom stats quantile
Sim_tau<-function(data,scaling,q,others){
  ncl=ncol(data)
  gap=ceiling(max(10,3*log(ncl)))
  res=c()
  for(j in 1:(ncl-gap)){res=rbind(res,c(j,j+gap))}
  intervals=res
  if(scaling==F){
    fs=Get_fs(data,intervals,f_cusum,as.double(q),others)
  }else{
    fs=Get_fs(data,intervals,f_cusum_lrv,as.double(q),others)
  }

  bz=quantile(fs,0.999)
  return(bz)
}

#'
Get_good_intervals <- function(bs){
  NN=max(bs)
  bs=cbind(bs,1)
  rend=c()
  while(sum(bs[,3])>0){
    s=min(bs[,2]/bs[,3])
    rend=c(rend,s)
    bs[(bs[,3]==1) & bs[,1]<s,3]=0
  }
  bs[,3]=1
  lend=c()
  while(sum(bs[,3])>0){
    s=max(bs[,1]*bs[,3])
    lend=c(lend,s)
    bs[(bs[,3]==1) & bs[,2]>s,3]=0
  }
  lend=rev(lend)
  if(length(rend)!=length(lend)){
    return(NULL)
  }
  for( j in 1:length(rend)){
    if((rend-lend)[j]<5){
      rend[j]=min(c(round(0.5*(rend[j]+lend[j]))+3,NN))
      lend[j]=max(c(round(0.5*(rend[j]+lend[j]))-3,0))
    }
  }
  return(matrix(c(lend,rend),ncol=2))
}


#' @importFrom purrr map_dbl
#' @importFrom stats density
densitypeak_cluster_tmp<-function(x,bw,adjust){
  x=as.vector(x)
  x=round(x,8)
  x=unique(x)
  a=density(x,bw=bw,adjust = adjust)

  x_rho=suppressWarnings(purrr::map_dbl(x,function(x1){
    y1=which(a$x>x1)[1]
    y2=y1-1
    return((a$y[y1]-a$y[y2])*((x1-a$x[y2])/(a$x[y1]-a$x[y2]))+a$y[y2])
  }))

  x_delta=suppressWarnings(purrr::map_dbl(1:length(x),function(j){
    min(abs(x[x_rho>x_rho[j]]-x[j]))
  }))

  x_delta[is.infinite(x_delta)]=max(x_delta[is.finite(x_delta)],na.rm=T)
  x=x[rev(order(x_delta))]
  x_rho=x_rho[rev(order(x_delta))]
  x_delta=sort(x_delta,decreasing = T)
  return(list(x=x,x_rho=x_rho,x_delta=x_delta))
}

#'
densitypeak_cluster<-function(density_list,k){
  x=density_list$x
  x_rho=density_list$x_rho
  x_delta=density_list$x_delta
  k=min(k,sum(x_delta>(min(x_delta)+max(x_delta))/4))
  k=max(k,2)
  res=list()
  for(i in 1:k){
    res[[i]]=x[i]
  }
  flag=1:k
  yy=x[1:k]
  x=x[rev(order(x_rho))]
  x=c(yy,x)
  for(j in (k+1):length(x)){
    r= which.min(abs(x[1:(j-1)]-x[j]))
    res[[flag[r]]]=c(res[[flag[r]]],x[j])
    flag=c(flag,flag[r])
  }
  x=x[3:length(x)]
  cluster=flag[3:length(flag)]
  x=x[order(cluster)]
  cluster=cluster[order(cluster)]
  return(list(x=x,cluster=cluster,k=k))
}

#'
Random_Intervals<-function (n, M,min.length=6){
  n <- as.integer(n)
  M <- as.integer(M)
  jj=0
  intervals_res=t(matrix(c(1,n)))
  while(jj<M){
    intervals <- matrix(sample(1:n,2*(M-jj),replace = T),ncol=2)
    intervals <- matrix(c(apply(intervals,1,min),apply(intervals,1,max)),ncol=2)
    ind=which(intervals[,2]-intervals[,1]<min.length)
    if(length(ind)>0){
      intervals=intervals[-ind,]
    }
    intervals_res=rbind(intervals_res,intervals)
    jj=nrow(intervals_res)-1
  }
  intervals_res=intervals_res[-1,]
  intervals_res
}

#'
pp<-function(data,M=1000,q=1,l=0,intervals=NULL,tau="ref"){
  M=nrow(intervals)
  fs=Get_fs_l(data,intervals,l,q,IProj)
  fs$index=1:M
  if(tau=="ref"){
    ncl=ncol(data)
    gap=ceiling(min(10,3*log(ncl)))
    interval_sim=t(rbind(1:(ncl-gap),(1+gap):ncl))
    sim=Get_fs_l(data,interval_sim,l,q,IProj)
    tau=quantile(sim$Max,0.999)
  }
  return(list(Threshold=tau,
              Statistic_values=fs$Max))
}

#' @importFrom pracma gramSchmidt
IProj<-function(s,e,l){
  vec=rep(1,(e-s)*(l+1))
  if(l>0){
    for(j in 1:l){
      vec[((e-s)*j+1):((e-s)*(j+1))]=((s+1):e)^j
    }
  }
  x=matrix(vec,ncol=l+1)
  y=gramSchmidt(x)$Q
  diag(e-s)-y%*%t(y)
}



