###########################################################################
# DBVSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

dbvsmsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10, 
                    lambda_seq = 2^{seq(-10, 15, length.out = 100)}, 
                    v_seq = NULL, Nofv = 100,
                    u_seq = NULL, Nofu = 100, 
                    kernel = c("linear", "gaussian"), kparam = 1, scale = FALSE,
                    criterion = c("0-1", "hinge"), cv_type = c("standard", "osr"), interaction = FALSE, optModel = FALSE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel)
  cv_type = match.arg(cv_type)
  criterion = match.arg(criterion)
  
  # Find a optimal lambda in first step
  cat("Step 1 \n")
  initial_fit = cv.ramsvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, 
                             lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = kparam,
                             scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  
  cat("\n")
  
  # Find a relevant variable for optimal lambda in second step
  cat("Step 2 \n")
  
  if (optModel) {
    threshold_opt = FALSE
  } else {
    threshold_opt = TRUE 
  }
  
  select_fit = threshold_fun.dbvsmsvm(initial_fit, v_seq = v_seq, Nofv = Nofv, u_seq = u_seq, Nofu = Nofu,
                                      cv_type = cv_type, criterion = criterion,
                                      interaction = interaction, optModel = threshold_opt, nCores = nCores, ...)
  
  # Find a optimal lambda under the selected variable in third step  
  selected_x = x[, select_fit$selected[1:NCOL(x)] == 1, drop = FALSE]
  
  out$selected = select_fit$selected
  out$lambda_cv_inform = initial_fit
  out$selection_cv_inform = select_fit$selection_inform
  
  if (interaction) {
    out$interaction_selection_cv_inform = select_fit$interaction_selection_inform
  }
  
  if (optModel) {
    if (!is.null(valid_x)) {
      selected_valid_x = valid_x[, selected == 1, drop = FALSE]
    } else {
      selected_valid_x = valid_x
    }
    cat("\n")
    cat("Find optimal model \n")
    # temporary estimates sigma
    # opt_sigma = kernlab::sigest(y ~ selected_x, frac = 1, scaled = FALSE)[3]
    final_fit = cv.ramsvm(x = selected_x, y = y, valid_x = selected_valid_x, valid_y = valid_y,
                          nfolds = nfolds, lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = kparam,
                          scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
    out$opt_model = final_fit$opt_model
    out$opt_valid_err = final_fit$opt_valid_err
    out$valid_err = final_fit$valid_err
  }
  out$cv_type = select_fit$cv_type
  out$call = call
  class(out) = "dbvsmsvm"  
  return(out)  
}
  
  
  
  
  
  