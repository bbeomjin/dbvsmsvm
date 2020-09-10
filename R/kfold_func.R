###########################################################################
# GBFSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

Kfold_msvm = function(x, y, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
                      gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                      scale = FALSE, criterion = c("0-1", "loss"), gd = TRUE, gd_scale = TRUE, optModel = FALSE, nCores = 1, ...)
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

  if (!is.numeric(lambda_seq)) {
    lambda = as.numeric(lambda_seq)
  }

  if (!is.numeric(kparam)) {
    kparam = as.numeric(kparam)
  }

  # The number of classes
  k = length(unique(y))

  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, kparam = kparam)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    gd_list = vector("list", 1)
    fold_list = NULL

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          msvm_fit = SRAMSVM_solve(x = x, y = y, gamma = gamma,
                                                   lambda = params$lambda[j], kernel = kernel,
                                                   kparam = params$kparam[j], ...)
                          pred_val = predict(msvm_fit, newx = valid_x)

                          # Compute the gradient with respect to x
                          if (gd) {
                            gd_x = gradient(alpha = msvm_fit$beta[[1]], x = x, y = y, scale = gd_scale,
                                            kernel = kernel, kparam = list(kparam))
                          } else {
                            gd_x = NULL
                          }
                          
                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val[[1]][[1]]) / length(valid_y)
                            err = 1 - acc
                          } else {
                            err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                          }
                          return(list(error = err, gd = gd_x))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    gd_list[[1]] = lapply(fold_err, "[[", "gd")
    opt_ind = min(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  } else {
    # set.seed(y[1])
    # fold_list = createFolds(y, k = nfolds, list = TRUE)
    fold_list = data_split(y, nfolds, k = k)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
    gd_list = vector("list", nfolds)

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
                            msvm_fit = SRAMSVM_solve(x = x_fold, y = y_fold, gamma = gamma,
                                                     lambda = params$lambda[j], kernel = kernel,
                                                     kparam = params$kparam[j], ...)
                            pred_val = predict(msvm_fit, newx = x_valid)

                            # Compute the gradient with respect to x
                            if (gd) {
                              gd_x = gradient(alpha = msvm_fit$beta[[1]], x = x_fold, y = y_fold, scale = gd_scale,
                                              kernel = kernel, kparam = list(kparam))
                            } else {
                              gd_x = NULL
                            }
                            
                            if (criterion == "0-1") {
                              acc = sum(y_valid == pred_val[[1]][[1]]) / length(y_valid)
                              err = 1 - acc
                            } else {
                              err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                            return(list(error = err, gd = gd_x))
                          }, mc.cores = nCores)
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
      gd_list[[i]] = lapply(fold_err, "[[", "gd")
    }
    valid_err = colMeans(valid_err_mat, na.rm = TRUE)
    opt_ind = min(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }

  out = list()
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$opt_model_gd = lapply(gd_list, "[[", opt_ind)
  out$fold_ind = fold_list
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$gamma = gamma
  out$scale = scale
  out$gd_scale = gd_scale
  if (optModel) {
    opt_model = SRAMSVM_solve(x = x, y = y, gamma = gamma,
                              lambda = opt_param$lambda, kernel = kernel, kparam = opt_param$kparam, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "GBFSMSVM"
  return(out)
}



