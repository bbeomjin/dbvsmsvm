###########################################################################
# GBFSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

threshold_fun = function(x, ...)
{
  UseMethod("threshold_fun")
}


threshold_fun.default = function(x, y, valid_x = NULL, valid_y = NULL, lambda = 1, nfolds = 10, thresh_Ngrid = 10,
                                 gamma = 0.5, kernel = c("linear", "radial", "poly", "spline"), kparam = c(1),
                                 scale = FALSE, nCores = 1, cv_type = c("original", "osr"), criterion = c("0-1", "loss"),
                                 interaction = FALSE, gd_scale = TRUE, nCores = 1, ...)
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

  if (!is.numeric(lambda)) {
    lambda = as.numeric(lambda)
  }

  if (!is.numeric(kparam)) {
    kparam = as.numeric(kparam)
  }

  # Initial fitting
  fit = SRAMSVM_solve(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, ...)

  # Compute the gradient with respect to x
  gd = gradient(alpha = fit$beta[[1]], x = x, y = y, scale = gd_scale, kernel = kernel, kparam = list(kparam))

  # Compute thresholding path
  gd_vec = c(0, seq(min(gd), max(gd), length.out = thresh_Ngrid))
  gd_vec = gd_vec[-c(length(gd_vec))]

  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(gd_vec,
                       function(thresh) {
                         msvm_fit = SRAMSVM_solve(x = x[, gd > thresh, drop = FALSE], y = y, gamma = gamma,
                                                 lambda = lambda, kernel = kernel, kparam = kparam, ...)

                         pred_val = predict(msvm_fit, newx = valid_x[, gd > thresh, drop = FALSE])

                         if (criterion == "0-1") {
                           acc = sum(valid_y == pred_val[[1]][[1]]) / length(valid_y)
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

    fold_list = data_split(y, nfolds, k = k)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec))

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]

      fold_err = mclapply(gd_vec,
                         function(thresh) {
                           # Initial fitting RAMSVM for computing the gradients
                           init_fit = SRAMSVM_solve(x = x_fold, y = y_fold, gamma = gamma, lambda = lambda,
                                                    kernel = kernel, kparam = kparam, ...)
                           init_gd = gradient(alpha = init_fit$beta[[1]], x = x_fold, y = y_fold, scale = gd_scale,
                                              kernel = kernel, kparam = list(kparam))

                           msvm_fit = SRAMSVM_solve(x = x_fold[, init_gd > thresh, drop = FALSE], y = y_fold, gamma = gamma,
                                                   lambda = lambda, kernel = kernel, kparam = kparam, ...)

                           pred_val = predict(msvm_fit, newx = x_valid[, init_gd > thresh, drop = FALSE])

                           if (criterion == "0-1") {
                             acc = sum(y_valid == pred_val[[1]][[1]]) / length(y_valid)
                             err = 1 - acc
                           } else {
                             err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                           }
                           return(err)
                         }, mc.cores = nCores)
      valid_err_mat[i, ] = unlist(fold_err)
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    valid_se = apply(valid_err_mat, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_ind_osr = max(which(valid_err <= (valid_err[opt_ind] + valid_se[opt_ind])))

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


threshold_fun.GBFSMSVM = function(object, thresh_Ngrid = 10, cv_type = c("original", "osr"), criterion = c("0-1", "loss"),
                                  interaction = FALSE, nCores = 1, ...)
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
    valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE)) / matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
  }

  # The number of classes
  k = length(unique(y))

  # Initial fitting
  fit = SRAMSVM_solve(x = x, y = y, gamma = gamma, lambda = lambda, kernel = kernel, kparam = kparam, ...)

  # Compute the gradient with respect to x
  gd = gradient(alpha = fit$beta[[1]], x = x, y = y, scale = gd_scale, kernel = kernel, kparam = list(kparam))

  # Compute thresholding path
  gd_vec = c(0, seq(min(gd), max(gd), length.out = thresh_Ngrid))
  gd_vec = gd_vec[-c(length(gd_vec))]

  if (!is.null(valid_x) & !is.null(valid_y)) {

    fold_err = mclapply(gd_vec,
                        function(thresh) {
                          msvm_fit = SRAMSVM_solve(x = x[, gd > thresh, drop = FALSE], y = y, gamma = gamma,
                                                   lambda = lambda, kernel = kernel, kparam = kparam, ...)

                          pred_val = predict(msvm_fit, newx = valid_x[, gd > thresh, drop = FALSE])

                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val[[1]][[1]]) / length(valid_y)
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

    fold_list = object$fold_ind
    nfolds = length(unique(fold_list))
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(gd_vec))
    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      # fold = fold_list[[i]]
      fold = which(fold_list == i)
      y_fold = y[-fold]
      x_fold = x[-fold, , drop = FALSE]
      y_valid = y[fold]
      x_valid = x[fold, , drop = FALSE]

      fold_err = mclapply(gd_vec,
                         function(thresh) {
                           # Pre-computed gradient
                           fold_gd = object$opt_model_gd[[i]]

                           msvm_fit = SRAMSVM_solve(x = x_fold[, fold_gd > thresh, drop = FALSE], y = y_fold, gamma = gamma,
                                                   lambda = lambda, kernel = kernel, kparam = kparam, ...)

                           pred_val = predict(fit_tmp, newx = x_valid[, fold_gd > thresh, drop = FALSE])

                           if (criterion == "0-1") {
                             acc = sum(y_valid == pred_val[[1]][[1]]) / length(y_valid)
                             err = 1 - acc
                           } else {
                             err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                           }
                           return(err)
                         }, mc.cores = nCores)
      valid_err_mat[i, ] = unlist(fold_err)
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    valid_se = apply(valid_err_mat, 2, sd) / sqrt(nfolds)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_ind_osr = max(which(valid_err <= (valid_err[opt_ind] + valid_se[opt_ind])))

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
