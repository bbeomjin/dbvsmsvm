###########################################################################
# SMSVM v1.2.1: R functions written by SangJun Lee, Yuwon Kim, Yoonkyung
# Lee, and Ja-Yong Koo.
###########################################################################
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

sramsvm = function(x, y, lambda, lambda_theta, kernel, kparam, theta = NULL, isCombined = FALSE, gamma = 0.5, 
                   criterion = "hinge", cv = FALSE, cv_type = "original", fold = 5, ...)
{
  # initialize
  cat('----------\n')
  cat("Initialize \n")
  cat('----------\n')
  cat("c-step...\n")
  
  cstepResult = cstep_ram(x = x, y = y, lambda = lambda, kernel = kernel, kparam = kparam,
                          theta = theta, criterion = criterion, cv = cv, fold = fold, gamma = gamma, ...)
  opt_lambda = cstepResult$opt_lambda

  # first update
  cat('------------\n')
  cat("First Update \n")
  cat('------------\n')
  if(isCombined == TRUE) # combined method
  {
    cat("theta-step and c-step...", "\n")
    thetastepResult = thetastep_ram(x = x, y = y, opt_lambda = opt_lambda, lambda_theta = lambda_theta, kernel = kernel,
                                    kparam = kparam, criterion = criterion, cv = cv, cv_type = cv_type, 
                                    fold = fold, isCombined = isCombined, 
                                    pretheta = theta, gamma = gamma, ...)
    opt_theta = thetastepResult$opt_theta
    final_model = thetastepResult$model
    return(list(cstep0 = cstepResult, thetastep1 = thetastepResult,
                opt_theta = opt_theta, model = final_model))
  }
  else if(isCombined == FALSE) # sequential method
  {
    cat("theta-step...", "\n")
    thetastepResult = thetastep_ram(x = x, y = y, opt_lambda = opt_lambda, lambda_theta = lambda_theta, kernel = kernel,
                                    kparam = kparam, criterion = criterion, cv = cv, cv_type = cv_type, 
                                    fold = fold, isCombined = isCombined, 
                                    pretheta = theta, gamma = gamma, ...)
    opt_theta = thetastepResult$opt_theta

    cat("c-step...", "\n")
    cstep1Result = cstep_ram(x = x, y = y, lambda = lambda, kernel = kernel, kparam = kparam,
                             theta = opt_theta, criterion = criterion, cv = cv, fold = fold, gamma = gamma, ...)
    final_model = cstep1Result$model
    return(list(cstep0 = cstepResult, thetastep1 = thetastepResult,
                cstep1 = cstep1Result, opt_theta = opt_theta, model = final_model))
  }
}

cstep_ram = function(x, y, lambda, gamma = 0.5, kernel, kparam = 1, theta = NULL, criterion = "hinge", cv = FALSE, fold = 5, ...)
{
  if((criterion != "0-1") && (criterion != "hinge"))
  {
    cat("ERROR: Only 0-1 and hinge can be used as criterion!", "\n")
    return(NULL)
  }
  len_lambda = length(lambda)
  y = as.integer(y)
  k = length(unique(y))

  ERR = matrix(0, len_lambda, 1)
  HIN = matrix(0, len_lambda, 1)

  kernel_list = list(type = kernel, par = kparam)
  x = as.matrix(x)
  anova_kernel = make_anovaKernel(x, x, kernel_list)

  if (is.null(theta))
  { theta = matrix(1, anova_kernel$numK, 1)}
  K = combine_kernel(anova_kernel, theta)

  if (cv)  # cross-validation
  {
    ran = data_split(y, fold)
    for(i.cv in 1:fold )
    {
      # cat("Leaving subset[", i.cv,"] out in",fold,"fold CV:","\n")
      omit = (ran == i.cv)
      x_train = x[!omit,]
      y_train = y[!omit]

      x_test = x[omit,]
      y_test = y[omit]

      row_index = 0

      subanova_kernel = make_anovaKernel(x_train, x_train, kernel_list)
      subK = combine_kernel(subanova_kernel, theta)
      
      subanova_kernel_test = make_anovaKernel(x_test, x_train, kernel_list)
      subK_test = combine_kernel(subanova_kernel_test, theta)
      
      # cat("lambda of length",len_lambda,"|")
      for(lam in lambda)
      {
        row_index = row_index + 1
        model = SRAMSVM_solve(y = y_train, gamma = gamma, lambda = lam, kernel = kernel, kparam = kparam, K = subK, ...)
        # fit_test = predict(model, x_test, newK = subK_test)[[1]]
        # model = ramsvm(x_train, y = y_train, gamma = 0.5, lambda = lam, kernel = kernel_type_ramsvm, kparam = 1 / sqrt(kparam))
        # fit_test2 = predict(model, x_test)[[1]]

        fit_test = predict(object = model, newdata = x_test, newK = subK_test)
        if (criterion == "0-1") {
          ERR[row_index] = (ERR[row_index] + (1 - (sum(y_test == fit_test[[1]][[1]]) / length(y_test))) / fold)
        } else {
          HIN[row_index] = (HIN[row_index] + ramsvm_hinge(y_test, fit_test$inner_prod, k = k, gamma = gamma) / fold)
        }

        # cat('*')
      }
      # cat("|\n")
    }
    cat("The minimum of average", fold, "fold cross-validated", criterion,
        "loss:","\n")
  }

  # choose the optimal index for lambda
  # if the optimal values are not unique, choose the largest value
  # assuming that lambda is in increasing order.
  if(criterion == "0-1")
  {
    optIndex = (len_lambda:1)[which.min(ERR[len_lambda:1])]
    cat(min(ERR),"\n")
  }
  else if(criterion == "hinge")
  {
    optIndex = (len_lambda:1)[which.min(HIN[len_lambda:1])]
    cat(min(HIN),"\n")
  }
  # choose the best model
  opt_lambda = lambda[optIndex]
  cat("The optimal lambda on log2 scale:", opt_lambda,"\n\n")

  # opt_model = msvm.compact(K, y, 2^opt_lambda, epsilon, epsilon.H)
  opt_model = SRAMSVM_solve(y = y, gamma = gamma, lambda = opt_lambda, kernel = kernel, kparam = kparam, K = K, ...)
  list(opt_lambda = opt_lambda, error = ERR, hinge = HIN, model = opt_model)
}

thetastep_ram = function(x, y, opt_lambda, gamma = 0.5, lambda_theta, kernel,
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


