#' @title Fit a Fast Scalable Non-negative Matrix Factorization using L-0 penalization
#'
#' @description Fits a grade of membership model using non-negative matrix factorization using L-0 penalization
#'              that can handle large scale input data and is fast and scalable. Also allows for unsupervised
#'              and semi-supervised set ups for factors/clusters.
#'              X=LF
#'
#' @param X The input data matrix with rows being the features and columns the samples.
#' @param K the number of clusters to fit. Must be an integer.
#' @param knownF The number of known factors to be used. Default is missing in which case the unsupervised NMF0
#'               method is used.
#' @param lambda1  The tuning parameter for the L-0 penalty on the loading matrix L
#' @param lambda2  The tuning parameter for the L-2 penalty on the loading matrix L
#' @param lambda3  The tuning parameter for the L-2 penalty on the factor matrix F
#' @param tol The relative tolerance which when met calls for stoppage in the optimization run
#' @param maxiter The maximum number of iterations for which to run the iterative updates in nmf0
#' @param verb If TRUE, prints the progress of the model fit
#' @param init_method The method for initializing the NMF method. Can be one of two values - random and svd.
#'                    when init_method=random, the initial updates to L and F are decided randomly
#' @param semi If FALSE, then enforces that the factor matrix is non-negative, else does not make that assumption,
#'             Default is FALSE.
#' @param hard_keep The maximum number of clusters that are representative of each sample.
#'
#' @return Outputs the best NMF0 fitted model for cluster K that includes the estimated L and F matrices, and the
#'         model likelihood.
#'
#' @references  Hussein Hazimeh and Rahul Mazumder.2018. Fast Best Subset Selection:
#'              Coordinate Descent and Local Combinatorial Optimization Algorithms.
#'              arXiv preprint arXiv:1803.01454.
#'
#' @examples
#'
#' data("ex.counts")
#' out <- nmf0(t(ex.counts), K=4, tol=0.001, maxiter=1000, verb=FALSE)
#' @export


nmf0 = function(Xmat, K=1L, knownF, lambda1 = 0.1, lambda2 = 0.1, lambda3 = 0.1,
                tol = 1e-04, maxiter = 1000, verb= TRUE,
                init_method = "random", semi=FALSE, hard_keep = NULL){

  N = ncol(Xmat)
  G = nrow(Xmat)

  if(!is.numeric(K) | length(K) > 1){
    stop("K must be an integer scalar entry")
  }

  if(K != floor(K)){
    stop("The value of K, the number of clusters to fit, must be an integer")
  }

  if(!is.matrix(Xmat) & !is.data.frame(Xmat)){
    stop("The input data must either be a matrix or a data frame")
  }

  if(ncol(Xmat) <=1 | nrow(Xmat) <= 1){
    stop("The number of rows or columns in input data must be greater than 1")
  }

  if (K > N | K > G){
    stop("The number of clusters K must be smaller than both the number of rows and number of columns of input data")
  }

  ###################      missing knownF, unsupervised NMF0 method        #########################


  if(missing(knownF)){
    message("knownF not provided, running unsupervised NMF0 method.")
    K2 = K
  }

  ###################   not missing knownF, supervised/semi-supervised NMF0 method        #########################

  if(!missing(knownF)){
    if(nrow(Xmat) !=  nrow(knownF)){
      stop("the number of features must match between input data and knownF known factor matrix.")
    }
    if(ncol(knownF) > 0){
      K2 = K + ncol(knownF)
    }else{
      K2 = K
    }
  }

  ######################  constraints the on the number of hard_kept clusters per sample    ######################

  if(!is.null(hard_keep)){
    if(hard_keep > nrow(Xmat) | hard_keep > ncol(Xmat) | hard_keep > K){
      stop("infeasible value of hard_keep, must be smaller than K, and also N and P.")
    }
  }


  ##############  Initializing the factor loading and the factor distribution matrices W and Q  ######################

  if(init_method == "svd"){
    res = nndsvd(Xmat, ncomponents = K)
    Wold = t(res$F1)
    if(missing(knownF)){
      Qold = res$L1
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K >= 1){
      Qold = cbind(res$L1[,1:K], knownF)
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K <= 0){
      Qold = knownF
    }else{
      stop("knownF if present must have more than 1 column, if not present K must be greater than 0")
    }
  }else if (init_method == "random"){
    avg = sqrt(mean(abs(Xmat))/K2)
    Wold = avg * matrix( abs(rnorm(N*K2,mean=0,sd=1)), nrow = K2, ncol = N)
    if(missing(knownF)){
      if(semi){
        Qold = avg * matrix(rnorm(G*K2,mean=0,sd=1), nrow = G, ncol = K2)
      }else{
        Qold = avg * matrix( abs(rnorm(G*K2,mean=0,sd=1)), nrow = G, ncol = K2)
      }
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K >= 1){
      if(semi){
        Qold = cbind(avg * matrix(rnorm(G*K,mean=0,sd=1), nrow = G, ncol = K), knownF)
      }else{
        Qold = cbind(avg * matrix( abs(rnorm(G*K,mean=0,sd=1)), nrow = G, ncol = K), knownF)
      }
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K <= 0){
      Qold = knownF
    }else{
      stop("knownF if present must have more than 1 column, if not present K must be greater than 0")
    }
  }else{
    stop("The init_method must be one of two options - svd or random")
  }


  ###############   Initial Frobenius norm for the initialization  ###################################

  init_normNMF = (matrixcalc::frobenius.norm(Xmat - Qold%*%Wold))/(ncol(Xmat)*nrow(Xmat))

  normNMFvec = c(init_normNMF)


  #################     Start of the iterative updates of the NMF0 method   ##########################


  for(num in 1:maxiter){

    ###################  updating W when Q is known  ############################

    svdQ = svd(Qold)
    QtQ = svdQ$v %*% diag(svdQ$d^2) %*% t(svdQ$v)
    Lipf = max(svdQ$d^2) + 2*lambda2
    coef =  -(1/Lipf)* QtQ + (Lipf - 2*lambda2)/(Lipf) * diag(1, nrow(QtQ))
    coef2 = (coef %*% Wold) + (1/Lipf)*t(Qold)%*%Xmat

    if(is.null(hard_keep)){
      thresh = sqrt(2 * lambda1/ Lipf)
      coef2[coef2 < thresh] = 0
      Wnew = coef2
    }else{
      coef2[coef2 < 0] = 0
      coef22 = apply(coef2, 2, function(x) {
        id = order(x, decreasing = TRUE)[1:hard_keep]
        z = rep(0, length(x))
        z[id] = x[id]
        return(z)
      })
      Wnew = coef22
    }


    ###################  updating Q when W is known   ##########################

    svdW = svd(Wnew)
    dd = svdW$d/(svdW$d^2 + lambda3)
    coef3 = t(svdW$u %*% diag(dd) %*% t(svdW$v) %*% t(Xmat))
    if(!semi){
      coef3[coef3 < 0] = 0
    }
    if(missing(knownF)){
      Qnew = coef3
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K >= 1){
      Qnew = cbind(coef3[,1:K], knownF)
    }else if (!missing(knownF) & ncol(knownF) >= 1 & K <= 0){
      Qnew = knownF
    }


    #################  calculating norm and checking if to continue the updating  ##################

    normNMF = matrixcalc::frobenius.norm(Xmat - Qnew%*%Wnew)/(ncol(Xmat)*nrow(Xmat))
    Qold = Qnew
    Wold = Wnew
    if(verb){cat("The Euclidean distance is :", normNMF, "\n")}
    tol_iter =  abs(tail(normNMFvec, 1) - normNMF)
    normNMFvec = rbind(normNMFvec, normNMF)
    if(tol_iter < tol){
      break
    }
  }
  ll = list("W" = Wnew, "Q" = Qnew, "norm_updates" = normNMFvec)
  return(ll)
}




nndsvd = function(Xmat, ncomponents){
  svd_out = rsvd::rsvd(Xmat, k=ncomponents)
  W = matrix(0, dim(svd_out$u)[1], dim(svd_out$u)[2])
  H = matrix(0, dim(svd_out$v)[1], dim(svd_out$v)[2])
  W[, 1] = sqrt(svd_out$d[1]) * abs(svd_out$u[,1])
  H[1, ] = sqrt(svd_out$d[1]) * abs(svd_out$v[1,])

  for(j in 2:ncomponents){
    x = svd_out$u[,j]
    y = svd_out$v[j,]
    x_p = pmax(x, 0)
    y_p = pmax(y, 0)
    x_n = abs(pmin(x, 0))
    y_n = abs(pmin(y, 0))
    x_p_nrm = sqrt(sum(x_p^2))
    y_p_nrm = sqrt(sum(y_p^2))
    x_n_nrm = sqrt(sum(x_n^2))
    y_n_nrm = sqrt(sum(y_n^2))
    m_p = x_p_nrm*y_p_nrm
    m_n = x_n_nrm*y_n_nrm
    if(m_p > m_n){
      u = x_p/x_p_nrm
      v = y_p/y_p_nrm
      sigma = m_p
    }else{
      u = x_n/x_n_nrm
      v = y_n/y_n_nrm
      sigma = m_n
    }
    lbd = sqrt(svd_out$d[j]*sigma)
    W[,j] = lbd * u
    H[j,] = lbd * v
  }
  ll=list("L1" = W, "F1" = H)
  return(ll)
}




check.matrix <- function(A, dm = NULL, mode = 'numeric', check.na = FALSE, input.name = '', check.negative = FALSE) {
  if (is.null(A)) return(invisible(NULL));
  if (!is.null(dm) && any(dim(A) != dm, na.rm = TRUE))
    stop(sprintf("Dimension of matrix %s is expected to be (%d, %d), but got (%d, %d)", input.name, nrow(A), ncol(A), dm[1], dm[2]));
  if (mode(A) != mode) stop(sprintf("Matrix %s must be %s.", input.name, mode));
  if (check.negative && any(A[!is.na(A)] < 0)) stop(sprintf("Matrix %s must be non-negative.", input.name));
  if (check.na && any(is.na(A))) stop(sprintf("Matrix %s contains missing values.", input.name));
}
