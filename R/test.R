gbfsmsvm2 = function(x, y, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
                    thresh_Ngrid = 10, kernel = "linear", kparam = 1, scale = FALSE, criterion = "0-1", cv_type = "original", interaction = FALSE,
                    gd_scale = TRUE, gamma = 0.5, optModel = FALSE, nCores = 1, ...)
{
  # Find a optimal lambda in first step
  initial_fit = Kfold_msvm2(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, lambda_seq = lambda_seq,
                           gamma = gamma, kernel = kernel, kparam = kparam, scale = scale, criterion = criterion,
                           gd = TRUE, optModel = FALSE, nCores = nCores, ...)
  
  
  # Find a relevant variable for optimal lambda in second step
  select_fit = threshold_fun2(initial_fit, thresh_Ngrid = thresh_Ngrid, cv_type = cv_type, criterion = criterion, gd_scale = gd_scale,
                             interaction = interaction, nCores = nCores, ...)
  selected = select_fit$selected
  
  # Find a optimal lambda under the selected variable in third step  
  selected_x = x[, selected == 1, drop = FALSE]
  
  out = list()
  out$selected = selected
  out$gd = select_fit$gd
  out$thresh_path = select_fit$thresh_path
  out$opt_thresh = select_fit$opt_thresh
  out$opt_valid_err = select_fit$opt_valid_err
  out$valid_err = select_fit$valid_err
  if (interaction) {
    out$inter_selected_var = select_fit$inter_selected_var
    out$gd_inter = select_fit$gd_inter
    out$opt_thresh = select_fit$inter_opt_thresh
  }
  
  if (optModel) {
    if (!is.null(valid_x)) {
      selected_valid_x = valid_x[, selected == 1, drop = FALSE]
    } else {
      selected_valid_x = valid_x
    }
    
    # temporary estimates sigma
    # opt_sigma = kernlab::sigest(y ~ selected_x, frac = 1, scaled = FALSE)[3]
    final_fit = Kfold_msvm(x = selected_x, y = y, valid_x = selected_valid_x, valid_y = valid_y,
                           nfolds = nfolds, lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = kparam,
                           scale = scale, criterion = criterion, gd = FALSE, optModel = TRUE, nCores = nCores, ...)
    out$opt_model = final_fit$opt_model
    out$valid_err = final_fit$valid_err
    out$opt_valid_err = final_fit$opt_valid_err
    
  }
  
  out$cv_type = select_fit$cv_type
  out$call = call
  class(out) = "GBFSMSVM"  
  return(out)  
}




Kfold_msvm2 = function (x, y, valid_x = NULL, valid_y = NULL, nfolds = 10, 
                        lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+06),
                        gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"),
                        kparam = c(1), scale = FALSE, 
                        criterion = c("0-1", "loss"), gd = TRUE,
                        optModel = FALSE, nCores = 1, ...) 
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)
  if (scale) {
    x = scale(x)
    if (!is.null(valid_x)) {
      means = attr(x, "scaled:center")
      stds = attr(x, "scaled:scale")
      valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE))/matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
    }
  }
  if (!is.numeric(lambda_seq)) {
    lambda = as.numeric(lambda_seq)
  }
  if (!is.numeric(kparam)) {
    kparam = as.numeric(kparam)
  }
  k = length(unique(y))
  params = expand.grid(lambda = lambda_seq, kparam = kparam)
  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL
    fold_err = mclapply(1:nrow(params), function(j) {
      msvm_fit = SRAMSVM_solve(x = x, y = y, gamma = gamma, 
                                          lambda = params$lambda[j], kernel = kernel, kparam = params$kparam[j], 
                                          ...)
      pred_val = predict(msvm_fit, newx = valid_x)
      # if (gd) {
      #   gd_x = gradient(alpha = msvm_fit$beta[[1]], x = x, 
      #                   y = y, scale = gd_scale, kernel = kernel, kparam = list(kparam))
      # }
      # else {
      #   gd_x = NULL
      # }
      if (criterion == "0-1") {
        acc = sum(valid_y == pred_val[[1]][[1]])/length(valid_y)
        err = 1 - acc
      }
      else {
        err = ramsvm_hinge(valid_y, pred_val$inner_prod, 
                                      k = k, gamma = gamma)
      }
      return(list(error = err, model = msvm_fit))
    }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "model")
    opt_ind = min(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }
  else {
    fold_list = data_split(y, nfolds, k = k)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
    model_list = vector("list", nfolds)
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i/nfolds * 100, "%", "\r")
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      fold_err = mclapply(1:nrow(params), function(j) {
        msvm_fit = SRAMSVM_solve(x = x_fold, y = y_fold, 
                                            gamma = gamma, lambda = params$lambda[j], kernel = kernel, 
                                            kparam = params$kparam[j], ...)
        pred_val = predict(msvm_fit, newx = x_valid)
        # if (gd) {
        #   gd_x = gradient(alpha = msvm_fit$beta[[1]], 
        #                   x = x_fold, y = y_fold, scale = gd_scale, 
        #                   kernel = kernel, kparam = list(kparam))
        # }
        # else {
        #   gd_x = NULL
        # }
        if (criterion == "0-1") {
          acc = sum(y_valid == pred_val[[1]][[1]])/length(y_valid)
          err = 1 - acc
        }
        else {
          err = ramsvm_hinge(y_valid, pred_val$inner_prod, 
                                        k = k, gamma = gamma)
        }
        return(list(error = err, model = msvm_fit))
      }, mc.cores = nCores)
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
      model_list[[i]] = lapply(fold_err, "[[", "model")
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }
  out = list()
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$opt_models = lapply(model_list, "[[", opt_ind)
  out$fold_ind = fold_list
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$gamma = gamma
  out$scale = scale
  if (optModel) {
    opt_model = SRAMSVM_solve(x = x, y = y, gamma = gamma, 
                              lambda = opt_param$lambda, kernel = kernel, kparam = opt_param$kparam, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "GBFSMSVM"
  return(out)
}



threshold_fun2 = function (object, thresh_Ngrid = 10, cv_type = c("original", "osr"),
                           criterion = c("0-1", "loss"), gd_scale = TRUE, interaction = FALSE, 
                           nCores = 1, ...) 
{
  call = match.call()
  cv_type = match.arg(cv_type)
  criterion = match.arg(criterion)
  x = object$x
  y = object$y
  valid_x = object$valid_x
  valid_y = object$valid_y
  lambda = object$opt_param$lambda
  kparam = object$opt_param$kparam
  gamma = object$gamma
  kernel = object$kernel
  if (object$scale & !is.null(valid_x)) {
    means = attr(x, "scaled:center")
    stds = attr(x, "scaled:scale")
    valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), 
                                byrow = TRUE))/matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
  }
  k = length(unique(y))
  fit = SRAMSVM_solve(x = x, y = y, gamma = gamma, lambda = lambda, 
                                 kernel = kernel, kparam = kparam, ...)
  gd = gradient(alpha = fit$beta[[1]], x = x, y = y, scale = gd_scale, 
                           kernel = kernel, kparam = list(kparam))
  gd_vec = c(0, seq(min(gd), max(gd), length.out = thresh_Ngrid))
  gd_vec = gd_vec[-c(length(gd_vec))]
  if (!is.null(valid_x) & !is.null(valid_y)) {
    fold_err = mclapply(gd_vec, function(thresh) {
      msvm_fit = SRAMSVM_solve(x = x[, gd > thresh, drop = FALSE], 
                                          y = y, gamma = gamma, lambda = lambda, kernel = kernel, 
                                          kparam = kparam, ...)
      pred_val = predict(msvm_fit, newx = valid_x[, gd > thresh, drop = FALSE])
      if (criterion == "0-1") {
        acc = sum(valid_y == pred_val[[1]][[1]])/length(valid_y)
        err = 1 - acc
      }
      else {
        err = ramsvm_hinge(valid_y, pred_val$inner_prod, 
                           k = k, gamma = gamma)
      }
      return(err)
    }, mc.cores = nCores)
    valid_err = unlist(fold_err)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_thresh = gd_vec[opt_ind]
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)
  }
  else {
    fold_list = object$fold_ind
    nfolds = length(unique(fold_list))
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec))
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i/nfolds * 100, "%", "\r")
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]
      fold_err = mclapply(gd_vec, function(thresh) {
        model = object$opt_models[[i]]
        
        fold_gd = gradient(model$beta[[1]], x_fold, y = y_fold, scale = gd_scale, kernel = kernel, kparam = list(kparam))
        model$x = x_fold[, fold_gd > thresh, drop = FALSE]
        # test_K = kernelMat(x_valid[, fold_gd > thresh, drop = FALSE], x_fold[, fold_gd > thresh, drop = FALSE],
        #                               kernel = kernel, kparam = kparam)
        
        # msvm_fit = SRAMSVM_solve(x = x_fold[, fold_gd > 
        #                                       thresh, drop = FALSE], y = y_fold, gamma = gamma, 
        #                          lambda = lambda, kernel = kernel, kparam = kparam, 
        #                          ...)
        pred_val = predict(model, newx = x_valid[, fold_gd > thresh, drop = FALSE])
        if (criterion == "0-1") {
          acc = sum(y_valid == pred_val[[1]][[1]])/length(y_valid)
          err = 1 - acc
        }
        else {
          err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
        }
        return(err)
      }, mc.cores = nCores)
      valid_err_mat[i, ] = unlist(fold_err)
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    valid_se = apply(valid_err_mat, 2, sd)/sqrt(nfolds)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_ind_osr = max(which(valid_err <= (valid_err[opt_ind] + 
                                            valid_se[opt_ind])))
    if (cv_type == "osr") {
      opt_thresh = gd_vec[opt_ind_osr]
    } else {
      opt_thresh = gd_vec[opt_ind]
    }
    opt_valid_err = min(valid_err)
    selected = as.integer(gd > opt_thresh)
  }
  out = list()
  out$selected = selected
  cat("The number of selected features out of ", length(selected), 
      ":", sum(selected), "\r", "\n")
  out$gd = gd
  out$gd_scale = gd_scale
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
