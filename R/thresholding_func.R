###########################################################################
# GBFSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

threshold_fun = function(x, ...)
{
  UseMethod("threshold_fun")
}


threshold_fun.default = function(x, y, valid_x = NULL, valid_y = NULL, lambda = 1, nfolds = 10, thresh_Ngrid = 10,
                                 gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                                 scale = FALSE, cv_type = c("original", "osr"), criterion = c("0-1", "loss"),
                                 interaction = FALSE, nCores = 1, ...)
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

  # The number of classes
  k = length(unique(y))

  lambda = as.numeric(lambda)
  kparam = as.numeric(kparam)
  # if (!is.numeric(lambda)) {
  #   lambda = as.numeric(lambda)
  # }
  # 
  # if (!is.numeric(kparam)) {
  #   kparam = as.numeric(kparam)
  # }

  # Initial fitting
  fit = ramsvm(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, ...)

  # Compute the gradient with respect to x
  gd = gradient(alpha = fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam)

  # Compute thresholding path
  gd_vec = seq(0, max(gd), length.out = thresh_Ngrid)
  gd_vec = gd_vec[-c(length(gd_vec))]

  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(gd_vec,
                       function(thresh) {
                         msvm_fit = ramsvm(x = x[, gd > thresh, drop = FALSE], y = y, gamma = gamma,
                                                 lambda = lambda, kernel = kernel, kparam = kparam, ...)

                         pred_val = predict.ramsvm(msvm_fit, newx = valid_x[, gd > thresh, drop = FALSE])

                         if (criterion == "0-1") {
                           acc = sum(valid_y == pred_val$class) / length(valid_y)
                           err = 1 - acc
                         } else {
                           err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                         }
                         return(err)
                       }, mc.cores = nCores)
    valid_err = unlist(fold_err)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_thresh = gd_vec[opt_ind]
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)

  } else {

    fold_list = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec))

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      
      # Initial fitting RAMSVM for computing the gradients
      init_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                        kernel = kernel, kparam = kparam, ...)
      
      init_gd = gradient(alpha = init_fit$cmat, x = x_fold, y = y_fold, kernel = kernel, kparam = kparam)
      
      fold_err = mclapply(gd_vec,
                         function(thresh) {
                           error = try({
                             msvm_fit = ramsvm(x = x_fold[, init_gd > thresh, drop = FALSE], y = y_fold, gamma = gamma,
                                               lambda = lambda, kernel = kernel, kparam = kparam, ...) 
                           })
                           
                           if (!inherits(error, "try-error")) {
                             pred_val = predict.ramsvm(msvm_fit, newx = x_valid[, init_gd > thresh, drop = FALSE])
                             
                             if (criterion == "0-1") {
                               acc = sum(y_valid == pred_val$class) / length(y_valid)
                               err = 1 - acc
                             } else {
                               err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                             }
                           } else {
                             err = Inf
                           }
                           return(err)
                         }, mc.cores = nCores)
      valid_err_mat[i, ] = unlist(fold_err)
    }
    valid_err = colMeans(valid_err_mat)
    valid_se = apply(valid_err_mat, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(valid_err == min(valid_err)))
    
    if (cv_type == "osr") {
      opt_ind_osr = max(which(valid_err <= (min(valid_err) + valid_se[opt_ind])))
      opt_thresh = gd_vec[opt_ind_osr]
    } else {
      opt_thresh = gd_vec[opt_ind]
    }
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)
  }

  out = list()
  out$selected = selected
  cat("The number of selected features out of ", length(selected), ":", sum(selected), "\r", "\n")
  out$gd = gd
  out$thresh_path = gd_vec
  out$opt_thresh = opt_thresh
  out$opt_valid_err = opt_valid_err
  if (is.null(valid_x) | is.null(valid_y)) {
    out$opt_valid_err_se = valid_se[opt_ind]
  }
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  if (interaction) {
    out$inter_selected_var = inter_selected_var
    out$gd_inter = gd_inter
    out$opt_thresh = inter_opt_thresh
  }
  out$cv_type = cv_type
  out$call = call
  return(out)
}


threshold_fun.dbvsmsvm = function(object, thresh_Ngrid = 10, cv_type = c("original", "osr"), criterion = c("0-1", "loss"),
                                  interaction = FALSE, nCores = 1, ...)
{
  call = match.call()
  cv_type = match.arg(cv_type)
  criterion = match.arg(criterion)

  x = object$x
  y = object$y
  valid_x = object$valid_x
  valid_y = object$valid_y
  lambda = object$opt_param["lambda"]
  kparam = object$opt_param["kparam"]
  gamma = object$gamma
  kernel = object$kernel

  if (object$scale & !is.null(valid_x)) {
    means = attr(x, "scaled:center")
    stds = attr(x, "scaled:scale")
    valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE)) / matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
  }

  # The number of classes
  k = length(unique(y))
  p = NCOL(x)
  
  # Initial fitting
  fit = ramsvm(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

  # Compute the gradient with respect to x
  gd = gradient(alpha = fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam)

  # Compute thresholding path
  gd_vec = seq(0, max(gd), length.out = thresh_Ngrid)
  gd_vec = gd_vec[-c(length(gd_vec))]

  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(gd_vec,
                        function(thresh) {
                          msvm_fit = ramsvm(x = x[, gd > thresh, drop = FALSE], y = y, gamma = gamma,
                                                   lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

                          pred_val = predict.ramsvm(msvm_fit, newx = valid_x[, gd > thresh, drop = FALSE])

                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val$class) / length(valid_y)
                            err = 1 - acc
                          } else {
                            err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                          }
                          return(err)
                        }, mc.cores = nCores)
    valid_err = unlist(fold_err)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_thresh = gd_vec[opt_ind]
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)
    
  } else {

    # fold_list = object$fold_ind
    nfolds = object$nfolds
    fold_list = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec))
    
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      
      # Initial fitting RAMSVM for computing the gradients
      init_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                        kernel = kernel, kparam = kparam, scale = FALSE, ...)
      
      init_gd = gradient(alpha = init_fit$cmat, x = x_fold, y = y_fold, kernel = kernel, kparam = kparam)
      
      fold_err = mclapply(gd_vec,
                         function(thresh) {
                           # Fit model under the fold set
                           error = try({
                             msvm_fit = ramsvm(x = x_fold[, init_gd > thresh, drop = FALSE], y = y_fold, gamma = gamma,
                                               lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...) 
                           })
                           
                           if (!inherits(error, "try-error")) {
                             pred_val = predict.ramsvm(msvm_fit, newx = x_valid[, init_gd > thresh, drop = FALSE])
                             
                             if (criterion == "0-1") {
                               acc = sum(y_valid == pred_val$class) / length(y_valid)
                               err = 1 - acc
                             } else {
                               err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                             } 
                           } else {
                             err = Inf 
                           }
                           return(err)
                         }, mc.cores = nCores)
      valid_err_mat[i, ] = unlist(fold_err)
    }
    valid_err = colMeans(valid_err_mat)
    valid_se = apply(valid_err_mat, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(valid_err == min(valid_err)))
    
    if (cv_type == "osr") {
      opt_ind_osr = max(which(valid_err <= (min(valid_err) + valid_se[opt_ind])))
      opt_thresh = gd_vec[opt_ind_osr]
    } else {
      opt_thresh = gd_vec[opt_ind]
    }
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)
  }

  out = list()
  out$selected = selected
  cat("The number of selected features out of ", length(selected), ":", sum(selected), "\r", "\n")
  out$gd = gd
  out$threshold_path = gd_vec
  out$opt_threshold = opt_thresh
  out$opt_valid_err = opt_valid_err
  if (is.null(valid_x) | is.null(valid_y)) {
    out$opt_valid_err_se = valid_se[opt_ind]
  }
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  
  if (interaction) {
    active_set = which(selected == 1)
    # comb_set = combn(1:NCOL(x), 2)
    if (length(active_set) == 1 | length(active_set) == 0) {
      int_selected = rep(0, choose(ncol(x), 2))
      gd_interaction = rep(0, choose(ncol(x), 2))
      opt_thresh_int = NULL
      int_opt_valid_err = NULL
      int_valid_err = NULL
    } else {
      gd_interaction = gradient_interaction(alpha = fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam, active_set = active_set)
      temp = combn(active_set, 2)
      gd_vec_int = seq(0, max(gd_interaction), length.out = thresh_Ngrid)
      # gd_vec_int = gd_vec_int[-c(length(gd_vec_int))]
      int_valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec_int))
      
      for (i in 1:nfolds) {
        cat(nfolds, "- fold CV (interaction) :", i / nfolds * 100, "%", "\r")
        fold = which(fold_list == i)
        y_fold = y[-fold]
        x_fold = x[-fold, , drop = FALSE]
        y_valid = y[fold]
        x_valid = x[fold, , drop = FALSE]
        
        # Initial fitting RAMSVM for computing the gradients
        init_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                          kernel = kernel, kparam = kparam, scale = FALSE, ...)
        fold_gd_int = gradient_interaction(alpha = init_fit$cmat, x = x_fold, y = y_fold, 
                                           kernel = kernel, kparam = kparam, active_set = active_set)
        
        fold_err_int = mclapply(gd_vec_int,
                                function(thresh) {
                                  KK = interaction_kernel(x_fold, x_fold, kernel = kernel, kparam = kparam), 
                                                          active_set, temp[, fold_gd_int > thresh, drop = FALSE])
                                  
                                  # Fit model under the fold set
                                  msvm_fit = ramsvm_solver(K = KK, y = y_fold, gamma = gamma,
                                                           lambda = lambda, kernel = kernel, kparam = kparam, ...)
                                  
                                  valid_KK = interaction_kernel(x_valid, x_fold, kernel = kernel, kparam = kparam, 
                                                                active_set, temp[, fold_gd_int > thresh, drop = FALSE])
                                  pred_val = predict.ramsvm(msvm_fit, newK = valid_KK)
                                  
                                  if (criterion == "0-1") {
                                    acc = sum(y_valid == pred_val$class) / length(y_valid)
                                    err = 1 - acc
                                  } else {
                                    err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                                  }
                                  return(err)
                                }, mc.cores = nCores)
        int_valid_err_mat[i, ] = unlist(fold_err_int)
      }
      int_valid_err = colMeans(int_valid_err_mat, na.rm = TRUE)
      int_valid_se = apply(int_valid_err_mat, 2, sd) / sqrt(nfolds)
      int_opt_ind = max(which(int_valid_err == min(int_valid_err)))
      int_opt_ind_osr = max(which(int_valid_err <= (int_valid_err[int_opt_ind] + int_valid_se[int_opt_ind])))
      
      if (cv_type == "osr") {
        opt_thresh_int = gd_vec_int[int_opt_ind_osr]
      } else {
        opt_thresh_int = gd_vec_int[int_opt_ind]
      }
      int_selected = as.integer(gd_interaction > opt_thresh_int)
    }
    out$int_selected = int_selected
    out$gd_interaction = gd_interaction
    out$opt_threshold_int = opt_thresh_int
    out$threshold_path_int = gd_vec_int
    out$int_opt_valid_err = min(int_valid_err)
    out$int_valid_err = int_valid_err
    out$se_int = int_valid_se
  }
  out$cv_type = cv_type
  out$call = call
  return(out)
}







