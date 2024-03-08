#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix pascal_triangle(int n) {
  NumericMatrix triangle(n, n);
  if (n >= 1) triangle(0, 0) = 1;
  if (n >= 2) {
    triangle(1, 0) = 1;
    triangle(1, 1) = 1;
  }
  for (int i = 2; i < n; ++i) {
    triangle(i, 0) = 1;
    triangle(i, i) = 1;
    for (int j = 1; j < i; ++j) {
      triangle(i, j) = triangle(i - 1, j - 1) + triangle(i - 1, j);
    }
  }
  return triangle;
}
// [[Rcpp::export]]
NumericMatrix matrixMultiplication(NumericMatrix mat1, NumericMatrix mat2) {
  Eigen::Map<Eigen::MatrixXd> eigenMat1(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(mat1));
  Eigen::Map<Eigen::MatrixXd> eigenMat2(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(mat2));
  Eigen::MatrixXd result = eigenMat1 * eigenMat2;
  Rcpp::NumericMatrix resultMat(result.rows(), result.cols());
  for (int i = 0; i < result.rows(); ++i) {
    for (int j = 0; j < result.cols(); ++j) {
      resultMat(i, j) = result(i, j);
    }
  }
  return resultMat;
}
// [[Rcpp::export]]
NumericVector fast_computation(NumericMatrix XIPu, int s, int e, int l, int q) {
  NumericMatrix coef = pascal_triangle(l + 1);
  int p = XIPu.nrow();
  int dim = e - s;
  NumericMatrix result(p,dim);
  for (int pp = 0; pp < p; ++pp) {
    NumericVector XIPu_sub = XIPu(pp,_);
    NumericMatrix res(l + 2, dim);
    for(int i=0;i<dim;++i){
      res(0, i) = XIPu_sub(i)*1.0 / pow(e - s, l);
    }
    res(0, 0) = 0;
    for (int ll = 0; ll <= l; ++ll) {
      res(ll + 1, dim - 1) = res(ll, dim - 1);
        for (int t = dim - 2; t >= 1; --t) {
          double sum = 0.0;
          for (int k = 1; k <= ll + 1; ++k) {
            sum += res(k, t + 1) * coef(ll, k - 1);
          }
          res(ll + 1, t) = sum + res(0, t);
        }
    }
    for(int i=0;i<dim;++i){
      result(pp,i)=res(l+1,i)*1.0*pow((e - s), l);
    }
  }
  NumericVector norm(e-s);
  for (int i = 0; i < e-s; ++i) {
    double sum = 0.0;
    for (int j = 0; j < p; ++j) {
      sum += pow(1.0*std::abs(result(j,i)), q);
    }
    norm(i) = pow(sum, 1.0 / q);
  }
  return norm;
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
DataFrame Get_fs_l(NumericMatrix data,IntegerMatrix intervals, int l, int q,Function IProj){
  int M = intervals.nrow();
  NumericVector ff(M);
  NumericVector pos(M);
  NumericVector mid(M);
  for(int j = 0; j < M; ++j) {
    int s = intervals(j, 0) - 1;
    int e = intervals(j, 1);
    if(e-s<=l+2){
      ff[j]=0;
      pos[j]=s+floor((e-s)*0.5);
      mid[j]=0.5*(s+e);
    }else{
    NumericMatrix IPu = IProj(s, e, l);
    NumericMatrix XIPu=matrixMultiplication(data(_,Range(s,e-1)),IPu);
    NumericVector Numerator=fast_computation(XIPu,s,e,l,q);
    NumericVector Denominator=fast_computation(IPu,s,e,l,2);
    NumericVector di = Numerator[Range(2, Numerator.size() - 1)] / Denominator[Range(2, Denominator.size() - 1)];
    ff[j] = max(di);
    pos[j] = which_max(di) + s+2;
    mid[j]=0.5*(s+e);
    }
  }
  return DataFrame::create(_["Left"] = intervals(_, 0),
                           _["Right"] = intervals(_, 1),
                           _["Max"] = ff,
                           _["Pos"] = pos,
                           _["mid"]=mid);
}
// [[Rcpp::export]]
NumericVector Get_fs(NumericMatrix x,NumericMatrix intervals,Function f, double q,Rcpp::Nullable<Rcpp::List> others = R_NilValue) {
  List customList;
  if (others.isNotNull()) {
    customList = as<List>(others);
  }else{
    customList = Rcpp::List::create();
  }
  int M=intervals.nrow();
  int s;
  int e;
  int j;
  NumericVector result(M);
  NumericVector column;
  double norm;
  for (j = 0; j < M; ++j) {
    s=intervals(j,0);
    e=intervals(j,1);
    column = f(x,s,e,q,customList);
    norm=sum(pow(abs(column),q));
    result[j] = pow(norm, 1.0 / q);
  }
  return result;
}
// [[Rcpp::export]]
NumericVector Get_fs_univariate(NumericMatrix x,NumericMatrix intervals,Function f, double q,Rcpp::Nullable<Rcpp::List> others = R_NilValue) {
  List customList;
  if (others.isNotNull()) {
    customList = as<List>(others);
  } else {
    customList = Rcpp::List::create();
  }
  int M=intervals.nrow();
  int s;
  int e;
  int j;
  int k;
  int p = x.nrow();
  NumericVector result(M);
  NumericVector column;
  double norm;
  for (j = 0; j < M; ++j) {
    s=intervals(j,0);
    e=intervals(j,1);
    norm=0.0;
    for(k=0;k<p;++k){
      column = f(x(k,_),s,e,q,customList);
      norm+=pow(abs(column[0]),q);
    }
    result[j] = pow(norm, 1 / q);
  }
  return result;
}
// [[Rcpp::export]]
NumericVector g_cusum(NumericVector x, int s, int e,Rcpp::Nullable<Rcpp::NumericVector> others = R_NilValue) {
  if (e <= s) {
    return NumericVector::create(0.0);
  }
  s = s - 1;
  int n = e - s;
  NumericVector tmp_cusum_value(n,0.0);
  double aa=sum(1.0*x[Range(s , e-1)]);
  double bb=0.0;
  for (int t = s; t < e-1; ++t) {
    bb+=x[t];
    aa-=x[t];
    tmp_cusum_value[t - s] = sqrt(1.0*(e-1 - t) / (e - s) / (t - s+1)) * bb -
      sqrt(1.0*(t - s+1) / (e - s) / (e-1 - t)) * aa;
  }
  return tmp_cusum_value;
}
// [[Rcpp::export]]
NumericVector f_cusum(NumericMatrix x, int s, int e, double q,List others) {
  int n = x.nrow();
  NumericMatrix nn(n, e - s + 1);
  NumericVector result(n);
  int j;
  for ( j = 0; j < n; ++j) {
    nn(j, _) = g_cusum(x(j, _), s, e);
  }
  double u;
  int index;
  double temp;
  u=sum(pow(abs(nn(_, 0)), q));
  index=0;
  for (j = 1; j < e - s + 1; ++j) {
    temp = sum(pow(abs(nn(_, j)), q));
    if(temp>u){
      u=temp;
      index=j;
    }
  }
  for (j=0;j<n;++j){
    result[j]=nn(j, index);
  }
  return result;
}
// [[Rcpp::export]]
NumericVector g_cusum_lrv(NumericVector x, int s, int e,NumericVector lrv) {
  if (e <= s) {
    return NumericVector::create(0.0);
  }
  s = s - 1;
  int n = e - s;
  NumericVector tmp_cusum_value(n,0.0);
  double sigma;
  double aa=sum(1.0*x[Range(s , e-1)]);
  double bb=0.0;
  double aa_lrv=sum(1.0*lrv[Range(s , e-1)]);
  double bb_lrv=0.0;
  for (int t = s; t < e-1; ++t) {
    bb+=x[t];
    aa-=x[t];
    bb_lrv+=lrv[t];
    aa_lrv-=lrv[t];
    tmp_cusum_value[t - s] = sqrt(1.0*(e-1 - t) / (e - s) / (t - s+1)) * bb -
      sqrt(1.0*(t - s+1) / (e - s) / (e-1 - t)) * aa;
    sigma=(1.0*(e-1-t)/(e-s)/(t-s+1))*bb_lrv+(1.0*(t-s+1)/(e-s)/(e-1-t))*aa_lrv;
    tmp_cusum_value[t - s]=tmp_cusum_value[t - s]/sqrt(sigma);
  }
  return tmp_cusum_value;
}
// [[Rcpp::export]]
NumericVector f_cusum_lrv(NumericMatrix x, int s, int e, double q,List others) {
  NumericMatrix lrv;
  int n = x.nrow();
  NumericMatrix nn(n, e - s + 1);
  NumericVector result(n);
  int j;
  lrv=as<NumericMatrix>(others[0]);
  for ( j = 0; j < n; ++j) {
    nn(j, _) = g_cusum_lrv(x(j, _), s, e,lrv(j,_));
  }
  double u;
  int index;
  double temp;
  u=sum(pow(abs(nn(_, 0)), q));
  index=0;
  for (j = 1; j < e - s + 1; ++j) {
    temp = sum(pow(abs(nn(_, j)), q));
    if(temp>u){
      u=temp;
      index=j;
    }
  }
  for (j=0;j<n;++j){
    result[j]=nn(j, index);
  }
  return result;
}
