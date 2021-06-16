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
                   lambda_theta = 2^{seq(-10, 10, length.out = 100)},
                   kernel, kparam, scale = FALSE, criterion = c("0-1", "loss"),
                   isCombined = FALSE, cv_type = "original", nCores = 1, ...)
{
  # initialize
  out = list()
  cat("Fit c-step \n")
  
  cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                            lambda_seq = lambda_seq, theta = NULL, kernel = kernel, kparam = kparam,
                            scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  
  cat("Fit theta-step \n")
  
  theta_step_fit = thetastep.sramsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                            lambda_seq = lambda_seq, theta = theta_step_fit$opt_theta, kernel = kernel, kparam = kparam,
                            scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)  
  
  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$cstep_valid_err = opt_cstep_fit$valid_err
  out$theta_valid_err = theta_step_fit$valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param["kparam"]
  out$opt_theta = theta_step_fit$opt_theta
  out$theta = theta_step_fit$theta
  out$x = x
  out$y = y
  out$n_class = opt_cstep_fit$n_class
  class(out) = "sramsvm"
  return(out)  
}

cstep.sramsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                         lambda_seq = 2^{seq(-10, 10, length.out = 100)}, theta = NULL,
                         kernel = c("linear", "radial", "poly"), kparam = c(1),
                         scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1,
                         type = "type1", ...)
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
      kernel_list = list(type = kernel, par = par)
      # anova_kernel = make_anovaKernel(x, x, kernel_list)
      # 
      # if (is.null(theta)) {
      #   theta = rep(1, anova_kernel$numK)
      # }
      
      # K = combine_kernel(anova_kernel, theta)
      
      ran = data_split(y, nfolds)
      valid_err_mat = matrix(NA, nrow = length(nfolds), ncol = length(lambda_seq))
      
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
        subK = combine_kernel(subanova_kernel, theta)
        
        subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel_list)
        subK_valid = combine_kernel(subanova_kernel_test, theta)
        
        fold_err = mclapply(1:length(lambda_seq),
                            function(j) {
                              error = try({
                                msvm_fit = ramsvm_fun(K = subK, y = train_y, gamma = gamma, lambda = lambda_sq[j], ...)
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
  
  if (optModel) {
    kernel_list = list(type = kernel, par = opt_param["kparam"])
    anova_K = make_anovaKernel(x, x, kernel_list)
    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }
    K = combine_kernel(anova_K, theta)
    opt_model = ramsvm_fun(K = K, y = y, gamma = gamma, lambda = opt_param["lambda"], kernel = kernel, kparam = opt_param["kparam"], ...)
    out$opt_model
  }
  out$call = call
  class(out) = "sramsvm"
  return(out)
}

thetastep.sramsvm = function(x, y, opt_lambda, gamma = 0.5, lambda_theta, kernel,
                         kparam = 1, criterion = "hinge",
                         cv = FALSE, fold = 5, isCombined = TRUE,
                         pretheta = NULL, cv_type = "original", ...)
{
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  y = as.integer(y)
  k = length(unique(y))

  len_lambda_theta = length(lambda_theta)
  ERR = matrix(0, len_lambda_theta, 1)
  HIN = matrix(0, len_lambda_theta, 1)

  kernel_list = list(type = kernel, par = kparam)
  x = as.matrix(x)
  anova_kernel = make_anovaKernel(x, x, kernel_list)

  if (is.null(pretheta)){
    pretheta = matrix(1, anova_kernel$numK, 1)
  }
  K = combine_kernel(anova_kernel, pretheta)
  
  # initial.model = msvm.compact(K, y, exp2.lambda, epsilon, epsilon.H)
  initial_model = SRAMSVM_solve(K = K, y = y, gamma = gamma, lambda = opt_lambda, kernel = kernel, kparam = kparam, ...)
  # ramsvm(x, y, lambda = exp2.lambda, kernel = "linear")@beta0
  theta_seq = matrix(0, len_lambda_theta, anova_kernel$numK)

  if (cv)  # cross-validation
  {
    ran = data_split(y, fold)
    ERR_mat = HIN_mat = matrix(NA, fold, len_lambda_theta)
    for(i.cv in 1:fold)
    {
      # cat("Leaving subset[", i.cv,"] out in",fold,"fold CV:","\n")
      omit = (ran == i.cv)
      x_train = x[!omit,]
      y_train = y[!omit]

      x_test = x[omit,]
      y_test = y[omit]

      subanova_kernel = make_anovaKernel(x_train, x_train, kernel_list)
      subanova_kernel_test = make_anovaKernel(x_test, x_train, kernel_list)
      subK = combine_kernel(subanova_kernel, pretheta)
      
      
      # model.initial = msvm.compact(subK, y_train, exp2.lambda, epsilon, epsilon.H)
      model_initial = SRAMSVM_solve(K = subK, y = y_train, gamma = gamma, lambda = opt_lambda,
                                    kernel = kernel, kparam = kparam, ...)

      row_index = 0
      # cat("lambda_theta of length",len_lambda_theta,"|")

      for(lam_theta in lambda_theta)
      {
        row_index = (row_index + 1)
        model = model_initial
        
        # find the optimal theta vector
        theta = find_theta(y = y_train, anova_kernel = subanova_kernel, gamma = gamma, cmat = model$beta[[1]], bvec = model$beta0[[1]],
                           lambda = opt_lambda, lambda_theta = lam_theta)
        # combine kernels
        subK = combine_kernel(subanova_kernel, theta)
        
        # combined method
        if(isCombined == TRUE) {
          # model = msvm.compact(subK, y_train, exp2.lambda, epsilon, epsilon.H)
          model = SRAMSVM_solve(K = subK, y = y_train, gamma = gamma, lambda = opt_lambda, kernel = kernel, kparam = kparam, ...)
        }
        subK_test = combine_kernel(subanova_kernel_test, theta)
        
        fit_test = predict(object = model, newK = subK_test)
        if (criterion == "0-1") {
          ERR[row_index] = (ERR[row_index] + (1 - (sum(y_test == fit_test[[1]][[1]]) / length(y_test))) / fold)
          ERR_mat[i.cv, row_index] = (1 - (sum(y_test == fit_test[[1]][[1]]) / length(y_test)))
        } else {
          HIN[row_index] = (HIN[row_index] + ramsvm_hinge(y_test, fit_test$inner_prod, k = k, gamma = gamma) / fold)
          HIN_mat[i.cv, row_index] = ramsvm_hinge(y_test, fit_test$inner_prod, k = k, gamma = gamma)
        }

        # cat('*')
      }
      # cat('|\n')
    }
    cat("The minimum of average", fold, "fold cross-validated", criterion,
        "loss:","\n")
    # generate a sequence of theta vectors
    model = initial_model
    row_index = 0
    for(lam_theta in lambda_theta)
    {
      row_index = (row_index + 1)
      theta = find_theta(y, anova_kernel, gamma = gamma, model$beta[[1]], model$beta0[[1]],
                         opt_lambda, lam_theta)
      theta_seq[row_index,] = theta
    }
  }

  # if the optimal values are not unique, choose the largest value
  # assuming that lambda_theta is in increasing order.
  if(criterion == "0-1")
  {
    if (cv_type == "original") {
      optIndex = (len_lambda_theta:1)[which.min(ERR[len_lambda_theta:1])]
    } else {
      err_cv_se = (apply(ERR_mat, 2, sd) / sqrt(fold))

      optIndex = (len_lambda_theta:1)[which.min(ERR[len_lambda_theta:1])]
      optIndex = max(which(ERR <= (min(ERR) + err_cv_se[optIndex])))
    }
    cat(min(ERR),"\n")

  } else if(criterion == "hinge") {
    if (cv_type == "original") {
      optIndex = (len_lambda_theta:1)[which.min(HIN[len_lambda_theta:1])]
    } else {
      hin_cv_se = (apply(HIN_mat, 2, sd) / sqrt(fold))

      optIndex = (len_lambda_theta:1)[which.min(HIN[len_lambda_theta:1])]
      optIndex = max(which(HIN <= (min(HIN) + hin_cv_se[optIndex])))
    }

    cat(min(HIN),"\n")
  }
  opt_lambda_theta = lambda_theta[optIndex]
  cat("The optimal lambda_theta on log2 scale:", opt_lambda_theta,"\n")

  opt_model = initial_model
  opt_theta = find_theta(y, anova_kernel, gamma = gamma, opt_model$beta[[1]], opt_model$beta0[[1]],
                         opt_lambda, opt_lambda_theta)
  nsel = sum(opt_theta > 0)
  shrinkage = mean(opt_theta)
  cat("The number of selected features out of", anova_kernel$numK, ":", nsel, "\n")
  cat("The average shrinkage factor:", shrinkage, "\n\n")

  K = combine_kernel(anova_kernel, opt_theta)
  
  if(isCombined == TRUE) #combined method
  {
    # opt_model = msvm.compact(K, y, exp2.lambda, epsilon, epsilon.H)
    opt_model = SRAMSVM_solve(K = K, y = y, gamma = gamma, lambda = opt_lambda,
                              kernel = kernel, kparam = kparam, ...)
  }
  list(lambda_theta = lambda_theta, opt_lambda_theta = opt_lambda_theta,
       error = ERR,  hinge = HIN, opt_theta = opt_theta, model = opt_model,
       nsel = nsel, shrinkage = shrinkage, theta_seq = theta_seq)
}


