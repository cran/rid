#' @title Convert a list of matrices into a single large matrix
#'
#' @param x A numeric list with each entity be a numeric matrix with p rows and q columns
#' @param sym A logical scalar representing whether each matrix is symmetric. If true, the duplicated half is removed
#'
#' @return A numeric matrix with pq rows and T columns
#'
#' @importFrom purrr map
#'
#' @examples
#' x=list(x1=1:3,x2=4:6)
#' List2Matrix(x)
#'
#' y=list(y1=matrix(1:4,2),y2=matrix(5:8,2))
#' List2Matrix(y)
#'
#' @export

List2Matrix <- function(x,sym=FALSE){
  stopifnot(is.list(x),is.numeric(x[[1]]))
  data=x
  if(sym){
    d <- purrr::map(data,function(a)a[lower.tri(a,diag = T)])
  }else{
    d <- purrr::map(data,function(a)as.vector(a))
  }
  return(t(matrix(unlist(d), ncol=length(d[[1]]), byrow=TRUE)))
}
