ssvm_m = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                  lambda_seq = 2^seq(-10, 10, length.out = 100), 
                  lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE,
                  kernel = c("linear", "radial", "poly", "spline", "anova_radial", "radial2"), kparam = 1, type = c("OVO", "OVR"),
                  scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, verbose = 1, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep_m.ssvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                           lambda_seq = lambda_seq, kernel = kernel, kparam = kparam, type = type, theta_mat = NULL, 
                           scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  
  cat("Fit theta-step \n")
  thetastep_fit = thetastep_m.ssvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, cv_type = cv_type, isCombined = isCombined, nCores = nCores, ...)
  
  cat("Fit c-step \n")
  opt_cstep_fit = cstep_m.ssvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                           lambda_seq = lambda_seq, kernel = kernel, kparam = kparam, type = type, theta_mat = thetastep_fit$opt_theta_mat, 
                           scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
  
  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
    cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
    cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
  }
  
  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$cstep_valid_err = opt_cstep_fit$valid_err
  out$theta_valid_err = thetastep_fit$valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param["kparam"]
  out$opt_theta = thetastep_fit$opt_theta_mat
  out$theta = thetastep_fit$theta_mat_seq
  out$x = x
  out$y = y
  out$n_class = opt_cstep_fit$n_class
  return(out)  
}


cstep_m_core.ssvm = function(x = NULL, y = NULL, lambda, theta_mat = NULL, kernel, kparam, type = c("OVO", "OVR"), ...)
{
  type = match.arg(type)
  n_class = length(unique(y))
  n = length(y)
  p = NCOL(x)
  out = list()
  
  if (kernel %in% c("radial2")) {
    numK = p + choose(p, 2)
  } else {
    numK = p
  }
  
  classname = levels(factor(y))
  if (is(y, "numeric")) {classname = as.numeric(classname)}
  
  if (type == "OVO") {
    # y_int = as.integer(y)
    comb = combn(n_class, 2)
    
    if (is.null(theta_mat)) {
      theta_mat = matrix(1, nrow = numK, ncol = choose(n_class, 2))
    }
    
    fitted_mat = matrix(0, nrow = n, ncol = n_class)
    colnames(fitted_mat) = classname
    rownames(fitted_mat) = as.character(1:n)
    model_list = list()
    
    for (j in 1:ncol(comb)) {
      theta = theta_mat[, j]
      index = y %in% classname[comb[, j]]
      yy = y[index]
      xx = x[index, ]
      
      subanova_K = make_anovaKernel(xx, xx, kernel, kparam)
      subK = combine_kernel(subanova_K, theta)
      svm_fit = svm_compact(K = subK, y = yy, lambda = lambda, ...)
      # svm_fit = svm_compact(K = subK, y = yy, lambda = lambda)
      model_list[[j]] = svm_fit
      fitted_mat[cbind(which(index), as.character(svm_fit$fit_class))] = fitted_mat[cbind(which(index), as.character(svm_fit$fit_class))] + 1
    }
    fit_class = classname[apply(fitted_mat, 1, which.max)]
    if (is(y, "factor")) {fit_class = factor(fit_class, classname)}
    
    out$theta_mat = theta_mat
    out$n_class = n_class
    out$x = x
    out$y = y
    out$classname = classname
    out$comb = comb
    out$kernel = kernel
    out$kparam = kparam
    out$models = model_list
    out$type = type
    out$fit_class = fit_class
  }
  
  if (type == "OVR") {
    
    if (is.null(theta_mat)) {
      theta_mat = matrix(1, nrow = numK, ncol = n_class)
    }
    
    fitted_mat = matrix(0, nrow = n, ncol = n_class)
    colnames(fitted_mat) = classname
    rownames(fitted_mat) = as.character(1:n)
    model_list = list()
    
    for (j in 1:n_class) {
      theta = theta_mat[, j]
      index = y %in% classname[j]
      yy = ifelse(index, 1, -1)
      subanova_K = make_anovaKernel(x, x, kernel, kparam)
      subK = combine_kernel(subanova_K, theta)
      svm_fit = svm_compact(K = subK, y = yy, lambda = lambda, ...)
      # svm_fit = svm_compact(K = subK, y = yy, lambda = lambda)
      model_list[[j]] = svm_fit
      nj = which(svm_fit$fit_class == 1)
      if (length(nj) == 0) {
        next
      } else {
        fitted_mat[cbind(nj, classname[j])] = 1
      }
    }
    fit_class = classname[apply(fitted_mat, 1, which.max)]
    if (is(y, "factor")) {fit_class = factor(fit_class, classname)}
    
    out$theta_mat = theta_mat
    out$n_class = n_class
    out$x = x
    out$y = y
    out$classname = classname
    out$kernel = kernel
    out$kparam = kparam
    out$models = model_list
    out$type = type
    out$fit_class = fit_class
  }
  return(out)
}

predict.cstep_m_core = function(object, newx, theta_mat = NULL)
{
  n = nrow(newx)
  pred_mat = matrix(0, nrow = n, ncol = object$n_class)
  classname = object$classname
  colnames(pred_mat) = classname
  rownames(pred_mat) = as.character(1:n)
  
  if (object$type == "OVO") {
    comb = object$comb
    if (is.null(theta_mat)) {
      theta_mat = object$theta_mat
    }
    
    for (j in 1:length(object$models)) {
      model = object$models[[j]]
      index = object$y %in% classname[comb[, j]]
      xx = object$x[index, ]
      theta = theta_mat[, j]
      
      K_valid = combine_kernel(make_anovaKernel(newx, xx, object$kernel, object$kparam), theta)
      pred_val = predict.svm_compact(model, newK = K_valid)$class
      pred_mat[cbind(1:n, as.character(pred_val))] = pred_mat[cbind(1:n, as.character(pred_val))] + 1
    }
    pred_class = classname[apply(pred_mat, 1, which.max)]
    if (is(y, "factor")) {pred_class = factor(pred_class, classname)}
  }
  
  if (object$type == "OVR") {
    if (is.null(theta_mat)) {
      theta_mat = object$theta_mat
    }
    
    for (j in 1:object$n_class) {
      model = object$models[[j]]
      theta = theta_mat[, j]
      
      K_valid = combine_kernel(make_anovaKernel(newx, object$x, object$kernel, object$kparam), theta)
      pred_val = predict.svm_compact(model, newK = K_valid)$class
      nj = which(pred_val == 1)
      if (length(nj) == 0) {
        next
      } else {
        pred_mat[cbind(as.character(nj), classname[j])] = 1
      }
    }
    pred_class = classname[apply(pred_mat, 1, which.max)]
    if (is(y, "factor")) {pred_class = factor(pred_class, classname)}
  }
  return(pred_class)
}


cstep_m.ssvm = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                        lambda_seq = 2^seq(-10, 10, length.out = 100), 
                        kernel = c("linear", "radial", "poly", "spline", "anova_radial", "radial2"), kparam = 1, type = c("OVO", "OVR"),
                        theta_mat = NULL, scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  type = match.arg(type)
  
  out = list()
  p = ncol(x)
  
  # The number of classes
  n_class = length(unique(y))
  
  lambda_seq = as.numeric(lambda_seq)
  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  kparam = as.numeric(kparam)
  
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq))
    
    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      fold_err = mclapply(1:length(lambda_seq), function(j) {
        error = try({
          cstep_fit = cstep_m_core.ssvm(x = train_x, y = train_y, lambda = lambda_seq[j], theta_mat = theta_mat,
                                        kernel = kernel, kparam = kparam, type = type, ...)
        })
        
        if (!inherits(error, "try-error")) {
          pred_val = predict.cstep_m_core(cstep_fit, newx = valid_x, theta_mat = theta_mat)
          if (criterion == "0-1") {
            acc = sum(valid_y == pred_val) / length(valid_y)
            err = 1 - acc
          } else {
            
          }
        } else {
          cstep_fit = NULL
          err = Inf
        }
        return(list(error = err, fit_model = cstep_fit))
      }, mc.cores = nCores)
      
      valid_err = sapply(fold_err, "[[", "error")
      valid_err_mat[i_cv, ] = valid_err
    }
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
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
  out$n_class = n_class
  out$theta_mat = theta_mat
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$type = type
  out$scale = scale
  out$criterion = criterion
  out$nfolds = nfolds
  if (optModel) {
    if (!is.null(valid_x) & !is.null(valid_y)) {
      out$opt_model = model_list[[opt_ind]]
    } else {
      opt_model = cstep_m_core.ssvm(x = x, y = y, lambda = opt_param["lambda"], theta_mat = theta_mat,
                                    kernel = kernel, kparam = kparam, type = type, ...)
      out$opt_model = opt_model
    }
  }
  out$call = call
  return(out)
}


thetastep_m.ssvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, cv_type = "original", isCombined = TRUE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$kparam
  type = object$type
  nfolds = object$nfolds
  
  x = object$x
  y = object$y
  
  valid_x = object$valid_x
  valid_y = object$valid_y
  
  if (is.null(object$opt_model)) {
    opt_model = cstep_m_core.ssvm(x = x, y = y, lambda = lambda, theta_mat = NULL, kernel = kernel, kparam = kparam, type = type, ...)
  } else {
    opt_model = object$opt_model
  }
  
  if (!is.null(valid_x) & !is.null(valid_y)) {
    
  } else {
    
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq))
    
    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      init_model = cstep_m_core.ssvm(x = train_x, y = train_y, lambda = lambda, theta_mat = object$theta_mat, 
                                     kernel = kernel, kparam = kparam, type = type, ...)
      
      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta_mat = findtheta_m.ssvm(y = train_y, x = train_x, models = init_model$models, 
                                                           lambda = lambda, lambda_theta = lambda_theta_seq[j], kernel = kernel, kparam = kparam, 
                                                           type = type)
                              
                              if (isCombined) {
                                init_model = cstep_m_core.ssvm(x = train_x, y = train_y, lambda = lambda, theta_mat = theta_mat, 
                                                               kernel = kernel, kparam = kparam, type = type, ...)
                              }
                            })
                            if (!inherits(error, "try-error")) {
                              pred_val = predict.cstep_m_core(init_model, valid_x, theta_mat = theta_mat)
                              acc = sum(pred_val == valid_y) / length(valid_y)
                              err = 1 - acc
                            } else {
                              err = Inf
                              theta_mat[] = 0
                            }
                            return(list(error = err, theta_mat = theta_mat))
                          }, mc.cores = nCores)
      
      valid_err = sapply(fold_err, "[[", "error")
      valid_err_mat[i_cv, ] = valid_err
    }
    valid_err = colMeans(valid_err_mat)
    
    if (cv_type == "original") {
      opt_ind = max(which(valid_err == min(valid_err)))
    } else {
      cv_se = (apply(valid_err_mat, 2, sd) / sqrt(nfolds))
      opt_ind = max(which(valid_err == min(valid_err)))
      opt_ind = max(which(valid_err <= (min(valid_err) + cv_se[opt_ind])))
    }
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)
    
    theta_mat_seq_list = mclapply(1:length(lambda_theta_seq),
                                  function(j) {
                                    theta_mat = findtheta_m.ssvm(y = y, x = x, models = opt_model$models,
                                                                 lambda = lambda, lambda_theta = lambda_theta_seq[j], type = type,
                                                                 kernel = kernel, kparam = kparam)
                                    return(theta_mat)
                                  })
    opt_theta_mat = theta_mat_seq_list[[opt_ind]]
  }
  
  out$opt_lambda_theta = opt_lambda_theta
  out$opt_ind = opt_ind
  out$opt_theta_mat = opt_theta_mat
  out$theta_mat_seq = theta_mat_seq_list
  out$opt_valid_err = opt_valid_err
  out$valid_err = valid_err
  return(out)
}



findtheta_m.ssvm = function(y, x, models, lambda, lambda_theta, kernel, kparam, type = c("OVO", "OVR"))
{
  call = match.call()
  type = match.arg(type)
  
  classname = levels(factor(y))
  if (is(y, "numeric")) {classname = as.numeric(classname)}
  n_class = length(classname)
  
  p = NCOL(x)
  if (kernel %in% c("radial2")) {
    numK = p + choose(p, 2)
  } else {
    numK = p
  }
  
  if (type == "OVO") {
    
    comb = combn(n_class, 2)
    theta_mat = matrix(NA, nrow = numK, ncol = ncol(comb))
    
    for (j in 1:ncol(comb)) {
      index = y %in% classname[comb[, j]]
      yy = y[index]
      xx = x[index, ]
      
      subanova_K = make_anovaKernel(xx, xx, kernel, kparam)
      
      alpha = models[[j]]$alpha
      bias = models[[j]]$bias
      
      theta = findtheta.ssvm(y = yy, anova_kernel = subanova_K, alpha, bias, lambda, lambda_theta)
      theta_mat[, j] = theta
    }
  }
  
  if (type == "OVR") {
    
    theta_mat = matrix(NA, nrow = numK, ncol = n_class)
    
    for (j in 1:n_class) {
      index = y %in% classname[j]
      yy = ifelse(index, 1, -1)
      subanova_K = make_anovaKernel(x, x, kernel, kparam)
      
      alpha = models[[j]]$alpha
      bias = models[[j]]$bias
      
      theta = findtheta.ssvm(y = yy, anova_kernel = subanova_K, alpha, bias, lambda, lambda_theta)
      theta_mat[, j] = theta
    }
  }
  return(theta_mat)
}




cstep.ssvm = function(x = NULL, y = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                      lambda_seq = 2^seq(-10, 10, length.out = 100), 
                      kernel = c("linear", "radial", "poly", "spline", "anova_radial", "radial2"), kparam = 1,
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
    
    anova_K = make_anovaKernel(x, x, kernel, kparam)
    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }
    K = combine_kernel(anova_K, theta)
    
    anova_K_valid = make_anovaKernel(valid_x, x, kernel, kparam)
    K_valid = combine_kernel(anova_K_valid, theta)
    
    fold_err = mclapply(1:length(lambda_seq),
                        function(j) {
                          error = try({
                            svm_fit = svm_compact(K = K, y = y, lambda = lambda_seq[j], ...)
                          })
                          
                          if (!inherits(error, "try-error")) {
                            pred_val = predict.svm_compact(svm_fit, K_valid)
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
    model_list = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
    
  } else {
    
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq))
    
    for (i_cv in 1:nfolds) {
      
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]
      
      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      if (is.null(theta)) {
        theta = rep(1, subanova_K$numK)
      }
      subK = combine_kernel(subanova_K, theta)
      
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)
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
    opt_ind = max(which(valid_err == min(valid_err)))
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
    
    if (!is.null(valid_x) & !is.null(valid_y)) {
      out$opt_model = model_list[[opt_ind]]
    } else {
      anova_K = make_anovaKernel(x, x, kernel, kparam)
      if (is.null(theta)) {
        theta = rep(1, anova_K$numK)
      }
      K = combine_kernel(anova_K, theta)
      opt_model = svm_compact(K = K, y = y, lambda = opt_param["lambda"], ...)
      out$opt_model = opt_model
    }
  }
  out$call = call
  class(out) = "ssvm"
  return(out)
}


thetastep.ssvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, cv_type = "original", nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  kernel = object$kernel
  kparam = object$kparam
  x = object$x
  y = object$y
  # theta = object$theta
  valid_x = object$valid_x
  valid_y = object$valid_y
  
  nfolds = object$nfolds
  
  anova_K = make_anovaKernel(x, x, kernel, kparam)
  
  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, object$theta)
    opt_model = svm_compact(K = K, y = y, lambda = lambda, ...)
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
      
      init_model = svm_compact(K = subK, y = train_y, lambda = lambda, ...)
      alpha = init_model$alpha
      bias = init_model$bias
      
      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = findtheta.ssvm(y = train_y, anova_kernel = subanova_K,
                                                     alpha = alpha, bias = bias, lambda = lambda, lambda_theta = lambda_theta_seq[j])
                              if (isCombined) {
                                subK = combine_kernel(subanova_K, theta)
                                init_model = svm_compact(K = subK, y = train_y, lambda = lambda, ...)
                              }
                            })
                            
                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.svm_compact(init_model, newK = subK_valid)
                              
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                
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
    
    if (cv_type == "original") {
      opt_ind = max(which(valid_err == min(valid_err)))
    } else {
      cv_se = (apply(valid_err_mat, 2, sd) / sqrt(nfolds))
      opt_ind = max(which(valid_err == min(valid_err)))
      opt_ind = max(which(valid_err <= (min(valid_err) + cv_se[opt_ind])))
    }
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)
    
    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = findtheta.ssvm(y = y, anova_kernel = anova_K, 
                                                         alpha = opt_model$alpha, bias = opt_model$bias, 
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
  out$vliad_err = valid_err
  return(out)
}


findtheta.ssvm = function(y, anova_kernel, alpha, bias, lambda, lambda_theta)
{
  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }
  
  y_temp = factor(y)
  classname = levels(y_temp)
  n_class = length(classname)
  
  y_int = integer(length(y))
  for (j in 1:n_class) {y_int[which(y_temp %in% classname[j])] = j}
  if (is(y, "numeric")) {classname = as.numeric(classname)}
  
  y_int = ifelse(y_int == 1, 1, -1)
  
  n = length(y_int)
  
  cvec = alpha * y_int
  
  dvec = NULL
  A_mat = NULL
  
  for (i in 1:anova_kernel$numK) {
    temp_d = (1 / 2) * lambda * t(cvec) %*% anova_kernel$K[[i]] %*% cvec + lambda_theta
    temp_A = y_int * as.vector(anova_kernel$K[[i]] %*% cvec)
    
    if (i == 1) {
      dvec = temp_d
      A_mat = temp_A
    } else {
      dvec = c(dvec, temp_d)
      A_mat = cbind(A_mat, temp_A)
    }
  }
  
  dvec = c(dvec, rep(1 / n, n))
  A_mat = cbind(A_mat, diag(1, n))
  # A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)
  
  bvec = c(1 - bias * y_int, rep(-1, anova_kernel$numK))
  
  # constraint directions
  const_dir = matrix(rep(">=", length(bvec)))
  
  # find solution by LP
  lp = lp("min", objective.in = dvec, const.mat = A_mat, const.dir = const_dir,
          const.rhs = bvec)$solution
  
  # find the theta vector only from the solution
  theta = lp[1:anova_kernel$numK]
  return(round(theta, 6))
}


