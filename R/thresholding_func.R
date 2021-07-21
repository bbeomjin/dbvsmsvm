###########################################################################
# GBFSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

threshold_fun = function(x, ...)
{
  UseMethod("threshold_fun")
}


threshold_fun.default = function(x, y, valid_x = NULL, valid_y = NULL, lambda = 1, nfolds = 10,
                                 v_seq = NULL, Nofv = 100, u_seq = NULL, Nofu = 100,
                                 gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                                 scale = FALSE, cv_type = c("original", "osr"), criterion = c("0-1", "loss"),
                                 interaction = FALSE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel)
  cv_type = match.arg(cv_type)
  criterion = match.arg(criterion)

  if (scale) {
    x = scale(x)
    if (!is.null(valid_x)) {
      means = attr(x, "scaled:center")
      stds = attr(x, "scaled:scale")
      valid_x = (valid_x - matrix(means, NROW(valid_x), NCOL(valid_x), byrow = TRUE)) / matrix(stds, NROW(valid_x), NCOL(valid_x), byrow = TRUE)
    }
  }

  # The number of classes
  k = length(unique(y))
  p = NCOL(x)
  
  lambda = as.numeric(lambda)
  kparam = as.numeric(kparam)
  
  # Initial fitting
  init_fit = ramsvm(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

  # Compute the partial derivatives with respect to x
  pderiv_vec = pderiv(alpha = fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam)

  # Compute threshold path
  if (is.null(v_seq)) {
    v_seq = seq(0, max(pderiv_vec), length.out = Nofv)
    v_seq = v_seq[-c(length(v_seq))]
  } 
  
  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(v_seq,
                       function(v) {
                         msvm_fit = ramsvm(x = x[, pderiv_vec > v, drop = FALSE], y = y, gamma = gamma,
                                                 lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

                         pred_val = predict.ramsvm(msvm_fit, newx = valid_x[, pderiv_vec > v, drop = FALSE])

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
    opt_v = v_seq[opt_ind]
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_v)

  } else {

    fold_list = data_split(y, nfolds)
    model_list = vector("list", nfolds)
    valid_err = matrix(NA, nrow = nfolds, ncol = length(v_seq), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      
      # Initial fitting RAMSVM for computing the partial derivatives
      fold_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                        kernel = kernel, kparam = kparam, scale = FALSE, ...)
      
      # Save the fitted model
      model_list[[i]] = fold_fit
      
      fold_pderiv_vec = pderiv(alpha = fold_fit$cmat, x = x_fold, y = y_fold, kernel = kernel, kparam = kparam)
      
      fold_err = mclapply(v_seq,
                         function(v) {
                           error = try({
                             msvm_fit = ramsvm(x = x_fold[, fold_pderiv_vec > v, drop = FALSE], y = y_fold, gamma = gamma,
                                               lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...) 
                           })
                           
                           if (!inherits(error, "try-error")) {
                             pred_val = predict.ramsvm(msvm_fit, newx = x_valid[, fold_pderiv_vec > v, drop = FALSE])
                             
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
      valid_err[i, ] = unlist(fold_err)
    }
    mean_valid_err = colMeans(valid_err)
    valid_se = apply(valid_err, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
    
    if (cv_type == "osr") {
      opt_ind_osr = max(which(mean_valid_err <= (min(mean_valid_err) + valid_se[opt_ind])))
      opt_v = v_seq[opt_ind_osr]
    } else {
      opt_v = v_seq[opt_ind]
    }
    opt_valid_err = min(mean_valid_err)
    selected = as.integer(pderiv_vec > opt_v)
  }

  out$selected = selected
  out$selection_inform = list(partial_deriv = pderiv_vec,
                              v_path = v_seq,
                              opt_v = opt_v,
                              opt_valid_err = opt_valid_err,
                              opt_ind = opt_ind,
                              valid_err = valid_err)
  
  if (interaction) {
    
    active_set = which(selected == 1)
    # comb_set = combn(1:NCOL(x), 2)
    if (length(active_set) == 1 | length(active_set) == 0) {
      interaction_selected = rep(0, choose(p, 2))
      so_pderiv_vec = rep(0, choose(p, 2))
      opt_u = NULL
      valid_err = NULL
      u_path = NULL
    } else {
      
      so_pderiv_vec = pderiv_so(alpha = init_fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam, active_set = active_set)
      
      if (is.null(u_seq)) {
        u_seq = seq(0, max(so_pderiv_vec), length.out = Nofu)
        u_seq = u_seq[-c(length(u_seq))]
      }
      
      temp = combn(active_set, 2)
      valid_err = matrix(NA, nrow = nfolds, ncol = length(u_seq), dimnames = list(paste0("Fold", 1:nfolds)))
      
      for (i in 1:nfolds) {
        cat(nfolds, "- fold CV (interaction) :", i / nfolds * 100, "%", "\r")
        fold = which(fold_list == i)
        y_fold = y[-fold]
        x_fold = x[-fold, , drop = FALSE]
        y_valid = y[fold]
        x_valid = x[fold, , drop = FALSE]
        
        # Initial fitting RAMSVM for computing the partial derivatives
        # init_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
        #                   kernel = kernel, kparam = kparam, scale = FALSE, ...)
        fold_fit = model_list[[i]]
        fold_so_pderiv_vec = pderiv_so(alpha = fold_fit$cmat, x = x_fold, y = y_fold, 
                                       kernel = kernel, kparam = kparam, active_set = active_set)
        
        fold_err = mclapply(u_seq,
                            function(u) {
                              KK = interaction_kernel(x_fold, x_fold, kernel = kernel, kparam = kparam, 
                                                      active_set, temp[, fold_so_pderiv_vec > u, drop = FALSE])
                              
                              # Fit model under the fold set
                              msvm_fit = ramsvm(K = KK, y = y_fold, gamma = gamma, lambda = lambda, ...)
                              
                              valid_KK = interaction_kernel(x_valid, x_fold, kernel = kernel, kparam = kparam, 
                                                            active_set, temp[, fold_so_pderiv_vec > u, drop = FALSE])
                              pred_val = predict.ramsvm(msvm_fit, newK = valid_KK)
                              
                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {
                                err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                              }
                              return(err)
                            }, mc.cores = nCores)
        valid_err[i, ] = unlist(fold_err)
      }
      mean_valid_err = colMeans(valid_err)
      valid_se = apply(valid_err, 2, sd) / sqrt(nfolds)
      opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
      opt_ind_osr = max(which(mean_valid_err <= (mean_valid_err[opt_ind] + valid_se[opt_ind])))
      
      if (cv_type == "osr") {
        opt_u = u_seq[opt_ind_osr]
      } else {
        opt_u = u_seq[opt_ind]
      }
      opt_valid_err = min(mean_valid_err)
      interaction_selected = as.integer(so_pderiv_vec > opt_u)
      comb_f = combn(1:p, 2)
      int_comb = temp[, interaction_selected == 1, drop = FALSE]
      interaction_selected = as.integer(paste0(comb_f[1, ], comb_f[2, ]) %in% paste0(int_comb[1, ], int_comb[2, ]))
    }
    out$selected = c(selected, interaction_selected)
    out$interaction_selection_inform = list(second_order_partial_deriv = so_pderiv_vec,
                                            u_path = u_seq,
                                            opt_u = opt_u,
                                            opt_valid_err = opt_valid_err,
                                            opt_ind = opt_ind,
                                            valid_err = valid_err)
  }
  out$cv_type = cv_type
  out$call = call
  cat(" The number of selected features out of ", length(selected), ":", sum(selected), "\r", "\n")
  return(out)
}


threshold_fun.dbvsmsvm = function(object, v_seq = NULL, Nofv = 100, u_seq = NULL, Nofu = 100,
                                  cv_type = c("original", "osr"), criterion = c("0-1", "loss"), 
                                  interaction = FALSE, nCores = 1, ...)
{
  out = list()
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
  
  # The number of classes
  k = length(unique(y))
  p = NCOL(x)
  
  # Initial fitting
  init_fit = ramsvm(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

  # Compute the partial derivatives with respect to x
  pderiv_vec = pderiv(alpha = init_fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam)

  # Compute thresholding path
  if (is.null(v_seq)) {
    v_seq = seq(0, max(pderiv_vec), length.out = Nofv)
    v_seq = v_seq[-c(length(v_seq))]
  } 

  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(v_seq,
                        function(v) {
                          msvm_fit = ramsvm(x = x[, pderiv_vec > v, drop = FALSE], y = y, gamma = gamma,
                                                   lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...)

                          pred_val = predict.ramsvm(msvm_fit, newx = valid_x[, pderiv_vec > v, drop = FALSE])

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
    opt_v = v_seq[opt_ind]
    opt_valid_err = min(valid_err)
    selected = as.integer(pderiv_vec > opt_v)
    
  } else {

    # fold_list = object$fold_ind
    nfolds = object$nfolds
    fold_list = data_split(y, nfolds)
    model_list = vector("list", nfolds)
    valid_err = matrix(NA, nrow = nfolds, ncol = length(v_seq), dimnames = list(paste0("Fold", 1:nfolds)))
    
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      
      # Initial fitting RAMSVM for computing the partial derivatives
      fold_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                        kernel = kernel, kparam = kparam, scale = FALSE, ...)
      
      # Save the fitted model
      model_list[[i]] = fold_fit
      
      # Compute the partial derivatives
      fold_pderiv_vec = pderiv(alpha = fold_fit$cmat, x = x_fold, y = y_fold, kernel = kernel, kparam = kparam)
      
      fold_err = mclapply(v_seq,
                         function(v) {
                           # Fit model under the fold set
                           error = try({
                             msvm_fit = ramsvm(x = x_fold[, fold_pderiv_vec > v, drop = FALSE], y = y_fold, gamma = gamma,
                                               lambda = lambda, kernel = kernel, kparam = kparam, scale = FALSE, ...) 
                           })
                           
                           if (!inherits(error, "try-error")) {
                             pred_val = predict.ramsvm(msvm_fit, newx = x_valid[, fold_pderiv_vec > v, drop = FALSE])
                             
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
      valid_err[i, ] = unlist(fold_err)
    }
    mean_valid_err = colMeans(valid_err)
    valid_se = apply(valid_err, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
    
    if (cv_type == "osr") {
      opt_ind_osr = max(which(mean_valid_err <= (min(mean_valid_err) + valid_se[opt_ind])))
      opt_v = v_seq[opt_ind_osr]
    } else {
      opt_v = v_seq[opt_ind]
    }
    opt_valid_err = min(mean_valid_err)
    selected = as.integer(pderiv_vec > opt_v)
  }
  
  out$selected = selected
  out$selection_inform = list(partial_deriv = pderiv_vec,
                              v_path = v_seq,
                              opt_v = opt_v,
                              opt_valid_err = opt_valid_err,
                              opt_ind = opt_ind,
                              valid_err = valid_err)
  
  if (interaction) {
    
    active_set = which(selected == 1)
    # comb_set = combn(1:NCOL(x), 2)
    if (length(active_set) == 1 | length(active_set) == 0) {
      interaction_selected = rep(0, choose(p, 2))
      so_pderiv_vec = rep(0, choose(p, 2))
      opt_u = NULL
      valid_err = NULL
      u_path = NULL
    } else {
      
      so_pderiv_vec = pderiv_so(alpha = init_fit$cmat, x = x, y = y, kernel = kernel, kparam = kparam, active_set = active_set)
      
      if (is.null(u_seq)) {
        u_seq = seq(0, max(so_pderiv_vec), length.out = Nofu)
        u_seq = u_seq[-c(length(u_seq))]
      }
      
      temp = combn(active_set, 2)
      valid_err = matrix(NA, nrow = nfolds, ncol = length(u_seq), dimnames = list(paste0("Fold", 1:nfolds)))
      
      for (i in 1:nfolds) {
        cat(nfolds, "- fold CV (interaction) :", i / nfolds * 100, "%", "\r")
        fold = which(fold_list == i)
        y_fold = y[-fold]
        x_fold = x[-fold, , drop = FALSE]
        y_valid = y[fold]
        x_valid = x[fold, , drop = FALSE]
        
        # Initial fitting RAMSVM for computing the partial derivatives
        # init_fit = ramsvm(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
        #                   kernel = kernel, kparam = kparam, scale = FALSE, ...)
        fold_fit = model_list[[i]]
        fold_so_pderiv_vec = pderiv_so(alpha = fold_fit$cmat, x = x_fold, y = y_fold, 
                                           kernel = kernel, kparam = kparam, active_set = active_set)
        
        fold_err = mclapply(u_seq,
                            function(u) {
                              KK = interaction_kernel(x_fold, x_fold, kernel = kernel, kparam = kparam, 
                                                          active_set, temp[, fold_so_pderiv_vec > u, drop = FALSE])
                                  
                              # Fit model under the fold set
                              msvm_fit = ramsvm(K = KK, y = y_fold, gamma = gamma, lambda = lambda, ...)
                                  
                              valid_KK = interaction_kernel(x_valid, x_fold, kernel = kernel, kparam = kparam, 
                                                            active_set, temp[, fold_so_pderiv_vec > u, drop = FALSE])
                              pred_val = predict.ramsvm(msvm_fit, newK = valid_KK)
                                  
                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {
                                err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                              }
                                return(err)
                              }, mc.cores = nCores)
        valid_err[i, ] = unlist(fold_err)
      }
      mean_valid_err = colMeans(valid_err)
      valid_se = apply(valid_err, 2, sd) / sqrt(nfolds)
      opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
      opt_ind_osr = max(which(mean_valid_err <= (mean_valid_err[opt_ind] + valid_se[opt_ind])))
      
      if (cv_type == "osr") {
        opt_u = u_seq[opt_ind_osr]
      } else {
        opt_u = u_seq[opt_ind]
      }
      opt_valid_err = min(mean_valid_err)
      interaction_selected = as.integer(so_pderiv_vec > opt_u)
      comb_f = combn(1:p, 2)
      int_comb = temp[, interaction_selected == 1, drop = FALSE]
      interaction_selected = as.integer(paste0(comb_f[1, ], comb_f[2, ]) %in% paste0(int_comb[1, ], int_comb[2, ]))
    }
    out$selected = c(selected, interaction_selected)
    out$interaction_selection_inform = list(second_order_partial_deriv = so_pderiv_vec,
                                            u_path = u_seq,
                                            opt_u = opt_u,
                                            opt_valid_err = opt_valid_err,
                                            opt_ind = opt_ind,
                                            valid_err = valid_err)
  }
  out$cv_type = cv_type
  out$call = call
  cat("The number of selected features out of ", length(selected), ":", sum(selected), "\r", "\n")
  return(out)
}







