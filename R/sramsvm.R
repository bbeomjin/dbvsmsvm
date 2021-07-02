# SMSVM using ramsvm
# library(quadprog)
# library(lpSolve)

cstep_ram_cv = function(x, y, lambda, kernel, kernel_par = 1,
                        theta = NULL, criterion = "hinge", x_test = NULL, y_test = NULL,
                        cv = FALSE, fold = 5)
{
  fit_list = vector("list", length(kernel_par))
  err_vec = rep(NA, length(kernel_par))
  for (i in 1:length(kernel_par)) {
    cstep_fit = cstep_ram(x = x, y = y, lambda = lambda, kernel = kernel, cv = cv, fold = folds,
                          criterion = criterion, kernel_par = kernel_par[i])
    fit_list[[i]] = cstep_fit
    err_vec[i] = min(cstep_fit$error)
  }
  min_err_ind = which.min(err_vec)
  final_model = fit_list[[min_err_ind]]
  return(list(model = final_model, opt_param = kernel_par[min_err_ind]))
}

sramsvm = function(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                   lambda_seq = 2^{seq(-10, 10, length.out = 100)},
                   lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                   kernel, kparam, scale = FALSE, criterion = c("0-1", "loss"),
                   isCombined = FALSE, cv_type = "original", type = c("type1", "type2"), nCores = 1, verbose = 1, ...)
{
  # initialize
  out = list()
  cat("Fit c-step \n")
  
  cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                            lambda_seq = lambda_seq, theta = NULL, kernel = kernel, kparam = kparam,
                            scale = scale, criterion = criterion, optModel = FALSE, type = type, nCores = nCores, ...)
  
  cat("Fit theta-step \n")
  
  theta_step_fit = thetastep.sramsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, cv_type = cv_type, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                            lambda_seq = lambda_seq, theta = theta_step_fit$opt_theta, kernel = kernel, kparam = kparam,
                            scale = scale, criterion = criterion, optModel = TRUE, type = type, nCores = nCores, ...)  
  
  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
    cat("CV-error(theta-step):", theta_step_fit$opt_valid_err, "\n")
    cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
  }
  
  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$cstep_valid_err = opt_cstep_fit$valid_err
  out$theta_valid_err = theta_step_fit$valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param["kparam"]
  out$opt_theta = theta_step_fit$opt_theta
  out$theta = theta_step_fit$theta_seq
  out$x = x
  out$y = y
  out$n_class = opt_cstep_fit$n_class
  class(out) = "sramsvm"
  return(out)  
}

cstep.sramsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                         lambda_seq = 2^{seq(-10, 10, length.out = 100)}, theta = NULL,
                         kernel = c("linear", "radial", "poly", "spline", "anova_radial", "radial2"), kparam = c(1),
                         scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1,
                         type = c("type1", "type2"), ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  type = match.arg(type)
  
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("ERROR: Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  
  if (type == "type1") {
    ramsvm_fun = ramsvm_solver
  } else {
    ramsvm_fun = ramsvm_compact
  }
  
  out = list()
  p = ncol(x)
  n_class = length(unique(y))
  
  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)
  
  # if (is.null(theta)) {
  #   theta = rep(1, p)
  # }
  
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  kparam = sort(kparam, decreasing = TRUE)
  
    
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    valid_err_by_kparam = vector("list", length(kparam))
    
    for (i in 1:length(kparam)) {
      par = kparam[i]
      # anova_kernel = make_anovaKernel(x, x, kernel, par)
      # 
      # if (is.null(theta)) {
      #   theta = rep(1, anova_kernel$numK)
      # }
      
      # K = combine_kernel(anova_kernel, theta)
      
      ran = data_split(y, nfolds)
      valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq))
      
      for (i_cv in 1:nfolds) {
        
        omit = ran == i_cv
        train_x = x[!omit, ]
        train_y = y[!omit]
        valid_x = x[omit, ]
        valid_y = y[omit]
        
        subanova_K = make_anovaKernel(train_x, train_x, kernel, par)
        if (is.null(theta)) {
          theta = rep(1, subanova_K$numK)
        }
        subK = combine_kernel(subanova_K, theta)
        
        subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, par)
        subK_valid = combine_kernel(subanova_K_valid, theta)
        
        fold_err = mclapply(1:length(lambda_seq),
                            function(j) {
                              error = try({
                                msvm_fit = ramsvm_fun(K = subK, y = train_y, gamma = gamma, lambda = lambda_seq[j], ...)
                              })
                              
                              if (!inherits(error, "try-error")) {
                                pred_val = predict.ramsvm_core(msvm_fit, subK_valid)
                                if (criterion == "0-1") {
                                  acc = sum(valid_y == pred_val$class) / length(valid_y)
                                  err = 1 - acc
                                } else {
                                  # 수정 필요 valid_y가 factor나 character일 경우
                                  err = ramsvm_hinge(valid_y, pred_val$pred_value, k = k, gamma = gamma)
                                }
                              } else {
                                msvm_fit = NULL
                                err = Inf
                              }
                              return(list(error = err, fit_model = msvm_fit))
                            }, mc.cores = nCores)
        
        valid_err = sapply(fold_err, "[[", "error")
        valid_err_mat[i_cv, ] = valid_err
      }
      valid_err_by_kparam[[i]] = colMeans(valid_err_mat)
    }
    valid_err_list = do.call(rbind, valid_err_by_kparam)
    opt_ind = which(valid_err_list == min(valid_err_list), arr.ind = TRUE)
    opt_ind = opt_ind[nrow(opt_ind), ]
    opt_param = c(lambda = lambda_seq[opt_ind[2]], kparam = kparam[opt_ind[1]])
    opt_valid_err = min(valid_err_list)
    valid_x = NULL
    valid_y = NULL
  }
  
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err_list
  out$x = x
  out$y = y
  out$gamma = gamma
  out$theta = theta
  out$n_class = n_class
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = opt_param["kparam"]
  out$scale = scale
  out$criterion = criterion
  out$nfolds = nfolds
  out$classifier = ramsvm_fun
  if (optModel) {
    anova_K = make_anovaKernel(x, x, kernel, opt_param["kparam"])
    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }
    K = combine_kernel(anova_K, theta)
    opt_model = ramsvm_fun(K = K, y = y, gamma = gamma, lambda = opt_param["lambda"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "sramsvm"
  return(out)
}

thetastep.sramsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, ...)
{
  out = list()
  call = match.call()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  kernel = object$kernel
  kparam = object$opt_param["kparam"]
  n_class = object$n_class
  x = object$x
  y = object$y
  gamma = object$gamma
  
  valid_x = object$valid_x
  valid_y = object$valid_y
  
  ramsvm_fun = object$classifier
  nfolds = object$nfolds
  
  # y_temp = as.factor(y)
  # y_name = levels(y_temp)
  # n_class = length(y_name)
  
  # y_int = integer(length(y))
  # for(j in 1:n_class) y_int[which(y_temp %in% y_name[j])] = j
  # if (is(y, "numeric")) {y_name = as.numeric(y_name)}
  
  anova_K = make_anovaKernel(x, x, kernel, kparam)
  # valid_anova_K = make_anovaKernel(valid_x, x, kernel, kparam)
  
  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, object$theta)
    opt_model = ramsvm_fun(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  } else {
    opt_model = object$opt_model
  }
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq))
    ran = data_split(y, nfolds)
    
    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      subK = combine_kernel(subanova_K, object$theta)
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)
      
      init_model = ramsvm_fun(K = subK, y = train_y, gamma = gamma, lambda = lambda, ...)
      cmat = init_model$cmat
      c0vec = init_model$c0vec
      
      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = findtheta.sramsvm(y = train_y, anova_kernel = subanova_K, gamma = gamma, cmat = cmat, c0vec = c0vec,
                                                        lambda = lambda, lambda_theta = lambda_theta_seq[j])
                              if (isCombined) {
                                subK = combine_kernel(subanova_K, theta)
                                init_model = ramsvm_fun(K = subK, y = train_y, gamma = gamma, lambda = lambda, ...)
                              }
                            })
                            
                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.ramsvm_core(init_model, newK = subK_valid)
                              
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                # 수정 필요 valid_y가 factor나 character일 경우
                                err = ramsvm_hinge(valid_y, pred_val$pred_value, k = k, gamma = gamma)
                              }
                            } else {
                              err = Inf
                              theta = rep(0, anova_K$numK)
                            }
                            return(list(error = err, theta = theta))
                          }, mc.cores = nCores)
      valid_err_mat[i_cv, ] = sapply(fold_err, "[[", "error")
    }
    valid_err = colMeans(valid_err_mat)
    
    # if the optimal values are not unique, choose the largest value
    # assuming that lambda_theta is in increasing order.
    if (cv_type == "original") {
      opt_ind = max(which(valid_err == min(valid_err)))
    } else {
      cv_se = (apply(valid_err_mat, 2, sd) / sqrt(nfolds))
      opt_ind = max(which(valid_err == min(valid_err)))
      opt_ind = max(which(valid_err <= (min(valid_err) + cv_se[opt_ind])))
    }
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)
    # opt_theta = findtheta.sramsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
    #                               lambda = lambda, lambda_theta = opt_lambda_theta)
    
    # generate a sequence of theta vectors
    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                       function(j) {
                         error = try({
                           theta = findtheta.sramsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
                                                     lambda = lambda, lambda_theta = lambda_theta_seq[j])
                         })
                         if (inherits(error, "try-error")) {
                           theta = rep(0, anova_K$numK)
                         }
                         return(theta)
                       }, mc.cores = nCores)
    theta_seq = do.call(cbind, theta_seq_list)
    opt_theta = theta_seq[, opt_ind]
  }
  
  out$opt_lambda_theta = opt_lambda_theta
  out$opt_ind = opt_ind
  out$opt_theta = opt_theta
  out$theta_seq = theta_seq
  out$opt_valid_err = opt_valid_err
  out$valid_err = valid_err
  return(out)
}

findtheta.sramsvm = function(y, anova_kernel, gamma = 0.5, cmat, c0vec, lambda, lambda_theta) 
{
  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }
  
  # standard LP form :
  # min a^T x , subject to A1x <= a1
  y_temp = as.factor(y)
  y_name = levels(y_temp)
  n_class = length(y_name)
  
  y_int = integer(length(y))
  for (j in 1:n_class) y_int[which(y_temp %in% y_name[j])] = j
  if (is(y, "numeric")) {y_name = as.numeric(y_name)}
  
  n = length(y_int)
  y_index = cbind(1:n, y_int)
  
  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, nobs = n, y = y_int)
  
  # calculate the 'a' matrix 
  a_tmp = matrix(gamma / n, nrow = n, ncol = n_class)
  a_tmp[y_index] = (1 - gamma) / n
  a = matrix(a_tmp, ncol = 1)
  
  # initialize M 
  M = matrix(rep(0, anova_kernel$numK), ncol = 1)
  
  # calculate M
  for (d in 1:anova_kernel$numK) {
    for (j in 1:(n_class - 1)) {
      M[d] = (M[d] + t(cmat[, j]) %*% anova_kernel$K[[d]] %*% cmat[, j])
    }
    M[d] = (lambda / 2 * M[d] + (lambda_theta))
  }
  a = rbind(a, M)
  
  Y_code = XI_gen(n_class)
  
  # calculate N matrix
  for (d in 1:anova_kernel$numK) {
    K = anova_kernel$K[[d]]
    for (j in 1:n_class) {
      if (j == 1) {
        temp_N = matrix(rowSums(K %*% cmat * matrix(Y_code[, j], nrow = n, ncol = ncol(cmat), byrow = TRUE)))
      } else {
        temp_N = rbind(temp_N, matrix(rowSums(K %*% cmat * matrix(Y_code[, j], nrow = n, ncol = ncol(cmat), byrow = TRUE))))
      }
    }
    if (d == 1) {
      N = temp_N
    } else {
      N = cbind(N, temp_N)
    }
  }
  
  sign_mat = matrix(1, n, n_class)
  sign_mat[y_index] = -1
  
  # constraints
  # A matrix
  I_nonzeroIndex = diag(1, n * n_class)
  N_nonzeroIndex = as.matrix(N) * as.vector(sign_mat)
  A_theta = cbind(matrix(0, anova_kernel$numK, n * n_class),
                  diag(-1, anova_kernel$numK))
  A_ineq = rbind(cbind(I_nonzeroIndex, -N_nonzeroIndex), A_theta)
  
  bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * c0vec[x]))
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(rowSums(t(Y_code) * matrix(c0vec, nrow = ncol(Y_code), ncol = nrow(Y_code), byrow = TRUE)), nrow = n, ncol = n_class, byrow = TRUE)
  bb_j[y_index] = bb_yi
  
  b_nonzeroIndex = as.vector(bb_j)
  b_ineq = rbind(as.matrix(b_nonzeroIndex), matrix(-1, anova_kernel$numK, 1))
  
  # constraint directions
  const_dir = matrix(rep(">=", nrow(matrix(b_ineq))))
  
  # find solution by LP
  lp = lp("min", objective.in = a, const.mat = A_ineq, const.dir = const_dir,
          const.rhs = b_ineq)$solution
  
  # find the theta vector only from the solution
  theta = cbind(matrix(0, anova_kernel$numK, n * n_class),
                diag(1, anova_kernel$numK)) %*% matrix(lp, ncol = 1)
  return(as.vector(theta))
}




