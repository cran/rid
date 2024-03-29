# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

pascal_triangle <- function(n) {
    .Call(`_rid_pascal_triangle`, n)
}

matrixMultiplication <- function(mat1, mat2) {
    .Call(`_rid_matrixMultiplication`, mat1, mat2)
}

fast_computation <- function(XIPu, s, e, l, q) {
    .Call(`_rid_fast_computation`, XIPu, s, e, l, q)
}

Get_fs_l <- function(data, intervals, l, q, IProj) {
    .Call(`_rid_Get_fs_l`, data, intervals, l, q, IProj)
}

Get_fs <- function(x, intervals, f, q, others = NULL) {
    .Call(`_rid_Get_fs`, x, intervals, f, q, others)
}

Get_fs_univariate <- function(x, intervals, f, q, others = NULL) {
    .Call(`_rid_Get_fs_univariate`, x, intervals, f, q, others)
}

g_cusum <- function(x, s, e, others = NULL) {
    .Call(`_rid_g_cusum`, x, s, e, others)
}

f_cusum <- function(x, s, e, q, others) {
    .Call(`_rid_f_cusum`, x, s, e, q, others)
}

g_cusum_lrv <- function(x, s, e, lrv) {
    .Call(`_rid_g_cusum_lrv`, x, s, e, lrv)
}

f_cusum_lrv <- function(x, s, e, q, others) {
    .Call(`_rid_f_cusum_lrv`, x, s, e, q, others)
}

