# dyn.load("./alpha_update2.dll")
ramsvm_solver = function(K = NULL, y, gamma = 0.5, lambda, 
                         weight = NULL, epsilon = 1e-4 * length(y) * length(unique(y)), maxit = 300)
{
  out = list()
  
  K = K + 1
  
  y_temp = as.factor(y)
  y_name = levels(y_temp)
  n_class = length(y_name)
  
  y_int = integer(length(y))
  for (j in 1:n_class) y_int[which(y_temp %in% y_name[j])] = j
  if (is(y, "numeric")) {y_name = as.numeric(y_name)}
  
  n = length(y_int)
  if (is.null(weight)) weight = numeric(n) + 1.0
  
  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  W = XI_gen(k = n_class)
  yyi = Y_matrix_gen(k = n_class,
                     n = n,
                     y = y_int)

  alpha_ij = matrix(data = 0.0, nrow = n, ncol = n_class)
  alpha_yi = numeric(n)

  erci = -diag(K) / 2 / as.double(n) / lambda

  aa = .C("alpha_update",
          as.vector(alpha_ij),
          as.vector(alpha_yi),
          as.vector(t(W)),
          as.vector(yyi),
          as.vector(K),
          as.double(lambda),
          as.vector(weight),
          as.integer(n),
          as.double(n),
          as.integer(n_class),
          as.double(n_class),
          as.vector(erci),
          as.double(gamma),
          as.vector(y_int),
          as.double(epsilon),
          outalpha_ij = as.vector(numeric(n * n_class)),
          maxiter = as.integer(maxit),
          PACKAGE = "GBFSMSVM")

  alpha = matrix(data = aa$outalpha_ij, nrow = n, ncol = n_class)

  cmat_list = coef_kernel(y = y_int,
                          k = n_class,
                          W = W,
                          alpha = alpha,
                          lambda = lambda)


  cmat = cmat_list$cmat
  c0vec = cmat_list$c0vec
  
  Kcmat = (K %*% cmat) %*% W
  W_c0vec = drop(t(c0vec) %*% W)
  # compute the fitted values
  fit = (matrix(W_c0vec, nrow = n, ncol = n_class, byrow = T) + Kcmat)
  fit_class = apply(fit, 1, which.max)
  
  out$y = y
  out$y_name = y_name
  out$alpha = alpha
  out$cmat = cmat
  out$c0vec = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n = n
  out$n_class = n_class
  out$gamma = gamma
  out$weight = weight
  out$lambda = lambda
  return(out)
}

ramsvm_compact = function(K, y, gamma = 0.5, lambda, epsilon = 1e-6, eig_tol_D = 0, epsilon_D = 1e-6)
{
  out = list()
  
  y_temp = as.factor(y)
  y_name = levels(y_temp)
  n_class = length(y_name)
  
  y_int = integer(length(y))
  for (j in 1:n_class) {y_int[which(y_temp %in% y_name[j])] = j}
  if (is(y, "numeric")) {y_name = as.numeric(y_name)}
  
  n = length(y_int)
  qp_dim = n * n_class
  
  code_mat = code_ramsvm(y_int)
  yyi = code_mat$yyi
  W = code_mat$W
  y_index = code_mat$y_index
  Hmatj = code_mat$Hmatj
  Lmatj = code_mat$Lmatj
  
  D = 0
  Amat = matrix(0, n * n_class, n_class - 1)
  for (j in 1:(n_class - 1)) {
    D = D + Hmatj[[j]] %*% K %*% t(Hmatj[[j]])
    Amat[, j] = -Lmatj[[j]]
  }
  
  D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(D))
  diag(D) = diag(D) + max_D * epsilon_D
  
  g_temp = matrix(-1, n, n_class)
  g_temp[y_index] = 1 - n_class
  g = as.vector(g_temp)
  
  dvec = -g * n * lambda
  
  Amat = cbind(Amat, diag(-1, n * n_class), diag(1, n * n_class))
  
  bvec_temp = matrix(gamma - 1, nrow = n, ncol = n_class)
  bvec_temp[y_index] = -gamma
  if (gamma == 0 | gamma == 1) {
    bvec_temp = bvec_temp - epsilon
  }
  # bvec = c(rep(0, qp_dim + n_class), as.vector(bvec_temp))
  bvec = c(rep(0, n_class - 1), as.vector(bvec_temp), rep(0, n * n_class))
  
  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind
  
  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, Amat, bvec, meq = n_class - 1)
  alpha = dual$solution
  alpha[alpha < 0] = 0
  
  alpha_mat = matrix(alpha, nrow = n, ncol = n_class)
  
  cmat = matrix(0, n, n_class - 1)
  for (j in 1:(n_class - 1)) {
    cmat[, j] = t(Hmatj[[j]]) %*% alpha / (n * lambda)
  }
  
  Kcmat = (K %*% cmat) %*% W
  
  alp_temp = matrix(1 - gamma, nrow = n, ncol = n_class)
  alp_temp[y_index] = gamma
  
  alp = c(as.vector(alp_temp), rep(0, 2 * (n_class - 1)))
  
  Alp1 = diag(qp_dim)
  Alp2 = matrix(0, nrow = qp_dim, ncol = 2 * (n_class - 1))
  
  for (i in 1:(n_class - 1)) {
    Alp2[, (2 * i - 1)] = Lmatj[[i]]
    Alp2[, (2 * i)] = -Lmatj[[i]]
  }
  
  Alp = cbind(Alp1, Alp2)
  
  blp_temp = Kcmat + 1
  blp_temp[y_index] = (n_class - 1) - Kcmat[y_index]
  blp = as.vector(blp_temp)
  
  const_dir = rep(">=", qp_dim)
  # const_dir[1] = "="
  cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir,const.rhs = blp)$solution[(qp_dim + 1):(qp_dim + 2 * (n_class - 1))]
  c0vec = rep(0, n_class - 1)
  for(j in 1:(n_class - 1)) {
    c0vec[j] = cposneg[(2 * j - 1)] - cposneg[(2 * j)]
  }
  
  W_c0vec = drop(t(c0vec) %*% W)
  # compute the fitted values
  fit = (matrix(W_c0vec, nrow = n, ncol = n_class, byrow = T) + Kcmat)
  fit_class = apply(fit, 1, which.max)
  # table(y, fit_class)
  
  # Return the output
  out$y = y
  out$y_name = y_name
  out$alpha = alpha_mat
  out$cmat = cmat
  out$c0vec = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n = n
  out$n_class = n_class
  out$gamma = gamma
  out$weight = NULL
  out$lambda = lambda
  return(out)
}


ramsvm = function(x = NULL, y, gamma = 0.5, lambda, kernel, kparam, scale = FALSE, type = c("type1", "type2"), ...)
{
  out = list()
  type = match.arg(type)
  # y_temp = as.factor(y)
  # y_name = levels(y_temp)
  n_class = length(unique(y))
  # 
  # y_int = integer(length(y))
  # for(j in 1:n_class) y_int[which(y_temp %in% y_name[j])] = j
  # if (is(y, "numeric")) {y_name = as.numeric(y_name)}
  
  n = NROW(x)
  p = ncol(x)
  
  center = rep(0, p)
  scaled = rep(1, p)
  
  if (scale) {
    x = scale(x)
    center = attr(x, "scaled:center")
    scaled = attr(x, "scaled:scale")
  }
  
  K = kernelMat(x, x, kernel = kernel, kparam = kparam)
  if (type == "type1") {
    solutions = ramsvm_solver(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  } else {
    solutions = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  }
  
  out$x = x
  out$y = y
  out$y_name = solutions$y_name
  out$gamma = gamma
  out$n_class = n_class
  out$lambda = lambda
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$c0vec = solutions$c0vec
  out$alpha = solutions$alpha
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  out$epsilon = list(...)
  class(out) = "ramsvm"
  return(out)
}


predict.ramsvm_core = function(object, newK = NULL) {
  
  cmat = object$cmat
  c0vec = object$c0vec
  n_class = object$n_class
  W = XI_gen(n_class)
  
  W_c0 = drop(t(c0vec) %*% W)
  
  fit = matrix(W_c0, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_y = apply(fit, 1, which.max)
  
  for(i in 1:object$n_class) {
    pred_y[pred_y == i] = object$y_name[i]
  }
  
  # return(list(class = pred_y, inner_prod = inner_prod))
  return(list(class = pred_y, pred_value = fit))
}


predict.ramsvm = function(object, newx = NULL, newK = NULL, ...) {

  if (is.null(newK)) {
    newK = kernelMat(newx, object$x, kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }
  
  cmat = object$cmat
  c0vec = object$c0vec
  n_class = object$n_class
  W = XI_gen(n_class)
  
  W_c0 = drop(t(c0vec) %*% W)
  fit = matrix(W_c0, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_y = apply(fit, 1, which.max)
  
  for(i in 1:object$n_class) {
    pred_y[pred_y == i] = object$y_name[i]
  }

  return(list(class = pred_y, pred_value = fit))
}


coef_kernel = function(y, k, W, alpha, lambda){
  
  n = length(y)
  
  cmat = matrix(data = 0, nrow = n, ncol = (k - 1))
  c0vec = numeric(k - 1)
  
  for (q in 1:(k - 1)) {
    temp = numeric(n)
    temp0 = 0
    for (i in 1:n) {
      for (j in 1:k) {
        t1 = alpha[i, j] * W[q, j]
        t2 = (2 * {y[i] == j} - 1)
        
        temp[i] = temp[i] + t2 * t1
        temp0 = temp0 + t2 * t1
      }
    }
    cmat[, q] = temp / n / lambda
    c0vec[q] = temp0 / n / lambda
  }
  return(list(cmat = cmat, c0vec = c0vec))
}

Kfold_ramsvm = function(x, y, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
                      gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                      scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  
  if (scale) {
    x = scale(x)
    if (!is.null(valid_x)) {
      means = attr(x, "scaled:center")
      stds = attr(x, "scaled:scale")
      valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE)) / matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
    }
  }
  
  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  kparam = sort(kparam, decreasing = TRUE)
  
  # The number of classes
  k = length(unique(y))
  
  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, kparam = kparam)
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL
    
    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          msvm_fit = ramsvm(x = x, y = y, gamma = gamma, lambda = params$lambda[j],
                                            kernel = kernel, kparam = params$kparam[j], ...)
                          pred_val = predict.ramsvm(msvm_fit, newx = valid_x)
                          
                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val$class) / length(valid_y)
                            err = 1 - acc
                          } else {
                            err = ramsvm_hinge(valid_y, pred_val$pred_value, k = k, gamma = gamma)
                          }
                          return(list(error = err, fit_model = msvm_fit))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = min(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  } else {
    # set.seed(y[1])
    # fold_list = createFolds(y, k = nfolds, list = TRUE)
    fold_list = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
    model_list = vector("list", nfolds)
    
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      
      #  Parallel computation on the combination of hyper-parameters
      fold_err = mclapply(1:nrow(params),
                          function(j) {
                            msvm_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = params$lambda[j], 
                                              kernel = kernel, kparam = params$kparam[j], ...)
                            pred_val = predict(msvm_fit, newx = x_valid)
                            
                            if (criterion == "0-1") {
                              acc = sum(y_valid == pred_val$class) / length(y_valid)
                              err = 1 - acc
                            } else {
                              err = ramsvm_hinge(y_valid, pred_val$pred_value, k = k, gamma = gamma)
                            }
                            return(list(error = err, fit_model = msvm_fit))
                          }, mc.cores = nCores)
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
      model_list[[i]] = lapply(fold_err, "[[", "fit_model")
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    # opt_ind = max(which(valid_err == min(valid_err)))
	opt_ind = min(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }
  
  out = list()
  out$opt_param = c(lambda = opt_param$lambda, kparam = opt_param$kparam)
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$fold_models = lapply(model_list, "[[", opt_ind)
  out$fold_ind = fold_list
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$gamma = gamma
  out$scale = scale
  if (optModel) {
    opt_model = ramsvm(x = x, y = y, gamma = gamma, lambda = opt_param$lambda,
                       kernel = kernel, kparam = opt_param$kparam, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "ramsvm"
  return(out)
}



