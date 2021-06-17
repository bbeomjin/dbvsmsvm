svm_compact = function(K = NULL, y, lambda, epsilon = 1e-6, epsilon_D = 1e-8)
{
  out = list()
  n = length(y)
  In = diag(1, n)
  K = K * (y %*% t(y))
  diag(K) = diag(K) + epsilon_D
  
  dvec = rep(1, n)
  Amat = cbind(y, In, -In)
  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind
  bvec = c(0, rep(0, n), rep(-1 / (n * lambda), n))
  alpha = solve.QP.compact(K, dvec, Amat, Aind, bvec, meq = 1, factorized = FALSE)$solution
  
  S = 1:n
  sv_index = S[alpha > epsilon]
  sv_number = length(sv_index)
  
  alpha[-sv_index] = 0
  svii = S[(alpha > epsilon) & (alpha < (1 / (n * lambda)))]
  if (length(svii) == 0) {
    svii = S[(alpha > epsilon) & (alpha <= (1 / (n * lambda)))]
  }
  bias = sum(y[svii] - K[svii, sv_index] %*% alpha[sv_index] * y[svii]) / length(svii)
  fit = bias + y * (K %*% alpha)
  fit[fit == 0] = 1e-10
  fit_class = sign(fit)
  
  out$alpha = alpha
  out$bias = bias
  out$fit = fit
  out$fit_class = fit_class
  out$n = n
  return(out)
}

svm = function(x = NULL, y, lambda, kernel, kparam, scale = FALSE, epsilon = 1e-6, epsilon_D = 1e-8)
{
  out = list()
  n = NROW(x)
  p = NCOL(x)
  
  center = rep(0, p)
  scaled = rep(1, p)
  
  if (scale) {
    x = scale(x)
    center = attr(x, "scaled:center")
    scaled = attr(x, "scaled:scale")
  }
  
  K = KernelMat(x, x, kernel = kernel, kparam = kparam)
  solutions = svm_compact(K = K, y = y, lambda = lambda, epsilon = epsilon, epsilon_D = epsilon_D)
  
  out$x = x
  out$y = y
  out$lambda = lambda
  out$kparam = kparam
  out$kernel = kernel
  out$alpha = solutions$alpha
  out$bias = solutions$bias
  out$epsilon = epsilon
  out$epsilon_D = epsilon_D
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "svm"
  return(out)
}


predict.svm_compact = function(object, newK = NULL)
{
  pred_y = bias + drop(t(object$y * object$alpha) %*% K)
  pred_y[pred_y == 0] = 1e-10
  pred_class = sign(pred_y)
  return(list(class = pred_class, pred_value = pred_y))
}


predict.svm = function(object, newx = NULL, newK = NULL)
{
  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }
  
  if (is.null(newK)) {
    newK = kernelMat(newx, object$x, kernel = object$kernel, kparam = object$kparam)
  }
  
  pred_y = bias + drop(t(object$y * object$alpha) %*% K)
  pred_y[pred_y == 0] = 1e-10
  pred_class = sign(pred_y)
  return(list(class = pred_class, pred_value = pred_y))
}



cstep.ssvm = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                      lambda_seq = 2^seq(-10, 10, length.out = 100), 
                      kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = 1,
                      theta = NULL, scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("ERROR: Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  
  out = list()
  p = NCOL(x)
  
  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)
  
  if (is.null(theta)) {
    theta = rep(1, p)
  }
  
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    kernel_list = list(type = kernel, par = kparam)
    
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq))
    
    for (i_cv in 1:nfolds) {
      
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      subanova_K = make_anovaKernel(train_x, train_x, kernel_list)
      if (is.null(theta)) {
        theta = rep(1, subanova_K$numK)
      }
      subK = combine_kernel(subanova_K, theta)
      
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel_list)
      subK_valid = combine_kernel(subanova_K_valid, theta)
      
      fold_err = mclapply(1:length(lambda_seq),
                          function(j) {
                            error = try({
                              svm_fit = svm_compact(K = subK, y = train_y, lambda = lambda_seq[j], ...)
                            })
                            
                            if (!inherits(error, "try-error")) {
                              pred_val = predict.svm_compact(svm_fit, subK_valid)
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                
                              }
                            } else {
                              svm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = svm_fit))
                          }, mc.cores = nCores)
      
      valid_err = sapply(fold_err, "[[", "error")
      valid_err_mat[i_cv, ] = valid_err
    }
    valid_err = colMeans(valid_err_mat)
    opt_ind = which(valid_err == min(valid_err))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
    valid_x = NULL
    valid_y = NULL
  }
  
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$theta = theta
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$scale = scale
  out$criterion = criterion
  out$nfolds = nfolds
  if (optModel) {
    kernel_list = list(type = kernel, par = kparam)
    anova_K = make_anovaKernel(x, x, kernel_list)
    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }
    K = combine_kernel(anova_K, theta)
    opt_model = svm_compact(K = K, y = y, lambda = opt_param["lambda"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "ssvm"
  return(out)
}


find.c.vec <- function(Y.vec, K.anova, L.mat, gamma.A, gamma.I, theta.vec, n.labeled, n.unlabeled, epsilon = 1e-7)
{
  # set matrix
  J.mat <- cbind(diag(1, n.labeled), matrix(0, n.labeled, n.unlabeled))
  Y.mat <- diag(as.vector(Y.vec), nrow(Y.vec), nrow(Y.vec))
  m.mat = 0
  temp.mat = 0 
  Q.mat = 0
  for(i in 1:K.anova$num.K)
  {
    temp.mat = temp.mat + theta.vec[i]*Y.mat%*%J.mat%*%K.anova$K[[i]]
    m.mat = m.mat + 2*(theta.vec[i]*gamma.A*K.anova$K[[i]] + gamma.I/(n.labeled+n.unlabeled)^2*theta.vec[i]^2*K.anova$K[[i]]%*%L.mat%*%K.anova$K[[i]])
  }
  
  if(gamma.I != 0 ) {
    m.mat = fixit(m.mat)
  }
  
  diag(m.mat) = diag(m.mat) + 1e-6
  #   m.mat = Ginv(m.mat)$Ginv
  #    m.mat = solvecov(m.mat)$inv
  m.mat = solve(m.mat)
  Q.mat = temp.mat%*%m.mat%*%t(temp.mat)
  #    K.mat = kernel.combine(K.anova, theta.vec)
  #    m.mat <- solve(2*gamma.A*diag(1, n.labeled+n.unlabeled)+2*gamma.I/(n.labeled+n.unlabeled)^2*L.mat%*%K.mat)
  #    Q.mat <- (Y.mat%*%J.mat%*%K.mat)%*%(m.mat)%*%t(J.mat)%*%Y.mat
  
  # solve QP
  Dmat = Q.mat
  Dmat = fixit(Dmat)
  diag(Dmat) = diag(Dmat) + 1e-6
  dvec = matrix(1, n.labeled, 1)
  Amat = cbind(matrix(Y.vec), diag(1, n.labeled), diag(-1, n.labeled))
  bvec = c(0, rep(0, n.labeled), rep(-1/n.labeled, n.labeled))
  beta = solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE)$solution
  
  # find c.vec
  # temp.mat = 0
  # for(i in 1:K.anova$num.K)
  # {
  #   temp.mat = temp.mat + theta.vec[i]*t(Y.mat%*%J.mat%*%K.anova$K[[i]])%*%matrix(beta, length(beta), 1)
  # }
  
  temp.mat = t(Y.mat %*% J.mat %*% kernel.combine(K.anova, theta.vec)) %*% matrix(beta, length(beta), 1)
  
  c.vec = m.mat%*%temp.mat
  
  # find bias
  K = kernel.combine(K.anova, theta.vec)
  S = 1:length(Y.vec)
  sv.index = S[beta > epsilon]
  sv.number = length(sv.index)
  #    sv = list(index = sv.index, number = sv.number)
  # let alpha = 0 if it is too small
  beta[-sv.index] = 0
  # compute bias
  svii = S[(beta > epsilon) & (beta < 1/n.labeled - epsilon)]
  
  if(length(svii) == 0) 
  {
    svii = S[(beta > epsilon)]
  }
  bias = sum(Y.vec[svii] - K[svii,] %*% c.vec)/length(svii)
  #    bias = sum(Y.vec[svii] - K[svii,] %*% c.vec * Y.vec[svii])/length(svii)
  #    bias = sum(Y.vec[svii] - K[svii, sv.index] %*% beta[sv.index] * Y.vec[svii]) / length(svii)
  #    bias = sum(Y.vec[svii] - K[svii, sv.index] %*% beta[sv.index]) / length(svii)
  #    bias = -(max(K[svii, sv.index] %*% beta[sv.index]) + min(K[svii, sv.index] %*% beta[sv.index]))/2
  return(list(c.vec = c.vec, bias = bias))
}
