###########################################################################
# DBVSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

dbvsmsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
                    thresh_Ngrid = 10, kernel = "linear", kparam = 1, scale = FALSE, criterion = "0-1", cv_type = "original", interaction = FALSE,
                    optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  # Find a optimal lambda in first step
  cat("Step 1 : ")
  initial_fit = Kfold_ramsvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                             lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = kparam,
                             scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  
  cat("\n")
  
  # Find a relevant variable for optimal lambda in second step
  cat("Step 2 : ")
  select_fit = threshold_fun.dbvsmsvm(initial_fit, thresh_Ngrid = thresh_Ngrid, cv_type = cv_type, criterion = criterion,
                                      interaction = interaction, nCores = nCores, ...)
  selected = select_fit$selected[1:NCOL(x)]
  
  # Find a optimal lambda under the selected variable in third step  
  selected_x = x[, selected == 1, drop = FALSE]
  
  out = list()
  out$selected = selected
  out$gd = select_fit$gd
  out$threshold_path = select_fit$threshold_path
  out$opt_threshold = select_fit$opt_threshold
  out$valid_err = select_fit$valid_err
  out$opt_valid_err = select_fit$opt_valid_err
  out$opt_valid_err_se = select_fit$opt_valid_err_se
  if (interaction) {
    out$int_selected = select_fit$int_selected
    out$gd_interaction = select_fit$gd_interaction
    out$threshold_path_int = select_fit$threshold_path_int
    out$opt_threshold_int = select_fit$opt_threshold_int
    out$int_valid_err = select_fit$int_valid_err
    out$int_opt_valid_err = select_fit$int_opt_valid_err
    # out$int_valid_se = select_fit$se_int
  }
  
  if (optModel) {
    if (!is.null(valid_x)) {
      selected_valid_x = valid_x[, selected == 1, drop = FALSE]
    } else {
      selected_valid_x = valid_x
    }
    cat("\n")
    cat("Find optimal model : ")
    # temporary estimates sigma
    # opt_sigma = kernlab::sigest(y ~ selected_x, frac = 1, scaled = FALSE)[3]
    final_fit = Kfold_ramsvm(x = selected_x, y = y, valid_x = selected_valid_x, valid_y = valid_y,
                           nfolds = nfolds, lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = kparam,
                           scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
    out$opt_model = final_fit$opt_model
    out$valid_err = final_fit$valid_err
    out$opt_valid_err = final_fit$opt_valid_err
  }
  out$cv_type = select_fit$cv_type
  out$call = call
  class(out) = "dbvsmsvm"  
  return(out)  
}
  
  
  
  
  
  