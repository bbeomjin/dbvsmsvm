svm_compact = function(K = NULL, y, lambda, epsilon = 1e-6, epsilon_D = 1e-8)
{
  out = list()
  
  if (sum(K) == 0) {
    diag(K) = 1
  }
  
  y_temp = factor(y)
  classname = levels(y_temp)
  n_class = length(classname)
  
  y_int = integer(length(y))
  for (j in 1:n_class) y_int[which(y_temp %in% classname[j])] = j
  if (is(y, "numeric")) {classname = as.numeric(classname)}
  
  y_int = ifelse(y_int == 1, 1, -1)
  
  n = length(y_int)
  In = diag(1, n)
  Q = K * (y_int %*% t(y_int))
  diag(Q) = diag(Q) + max(Q) * epsilon_D
  
  dvec = rep(1, n)
  Amat = cbind(y_int, In, -In)
  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind
  bvec = c(0, rep(0, n), rep(-1 / (n * lambda), n))
  alpha = solve.QP.compact(Q, dvec, Amat, Aind, bvec, meq = 1, factorized = FALSE)$solution
  
  S = 1:n
  sv_index = S[alpha > epsilon]
  sv_number = length(sv_index)
  
  alpha[-sv_index] = 0
  svii = S[(alpha > epsilon) & (alpha < (1 / (n * lambda)))]
  if (length(svii) == 0) {
    svii = S[(alpha > epsilon) & (alpha <= (1 / (n * lambda)))]
  }
  
  # compute the bias term
  
  # bias = sum(y_int[svii] - K[svii, sv_index] %*% alpha[sv_index] * y_int[svii]) / length(svii)
  # one = rep(1, n)
  # dalpha = diag(n * lambda * alpha)
  # bias = as.vector((t(one) %*% dalpha %*% (diag(1, n) - dalpha) %*% (y_int - K %*% (y_int * alpha))) / (t(n * lambda * alpha) %*% (one - n * lambda * alpha)))
  
  alp = c(rep(1 / n, n), 0, 0)
  # Alp1 = c(rep(0, n), c(1, -1))
  Alp1 = diag(n)
  Alp2 = cbind(y_int, -y_int)
  Alp = cbind(Alp1, Alp2)
  
  Kmat = as.vector(K %*% (y_int * alpha))
  
  blp = 1 - y_int * Kmat
  const_dir = rep(">=", n)
  cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir, const.rhs = blp)$solution
  bias = -(cposneg[n + 1] - cposneg[n + 2])
  
  fit = bias + Kmat
  fit[fit == 0] = 1e-10
  fit_class = sign(fit)
  
  fit_class = ifelse(fit_class == 1, classname[1], classname[2])
  if (is(y, "factor")) {fit_class = factor(fit_class, classname)}
  
  out$y = y
  out$y_int = y_int
  out$classname = classname
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
  
  K = kernelMatrix(x, x, kernel = kernel, kparam = kparam)
  solutions = svm_compact(K = K, y = y, lambda = lambda, epsilon = epsilon, epsilon_D = epsilon_D)
  
  out$x = x
  out$y = y
  out$y_int = solutions$y_int
  out$classname = solutions$classname
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
  pred_y = object$bias + as.vector(newK %*% (object$y_int * object$alpha))
  pred_y[pred_y == 0] = 1e-10
  pred_class = sign(pred_y)
  pred_class = ifelse(pred_class == 1, object$classname[1], object$classname[2])
  if (is(object$y, "factor")) {pred_class = factor(pred_class, object$classname)}
  if (is(object$y, "numeric")) {pred_class = as.numeric(pred_class)}
  return(list(class = pred_class, pred_value = pred_y))
}


predict.svm = function(object, newx = NULL, newK = NULL)
{
  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }
  
  if (is.null(newK)) {
    newK = kernelMatrix(newx, object$x, kernel = object$kernel, kparam = object$kparam)
  }
  
  pred_y = object$bias + as.vector(newK %*% (object$y * object$alpha))
  pred_y[pred_y == 0] = 1e-10
  pred_class = sign(pred_y)
  pred_class = ifelse(pred_class == 1, object$classname[1], object$classname[2])
  if (is(object$y, "factor")) {pred_class = factor(pred_class, object$classname)}
  if (is(object$y, "numeric")) {pred_class = as.numeric(pred_class)}
  return(list(class = pred_class, pred_value = pred_y))
}
