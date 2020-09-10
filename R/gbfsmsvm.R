###########################################################################
# GBFSMSVM v1.0.0: R functions written by Beomjin Park
###########################################################################

gbfsmsvm = function(x, y, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
                    thresh_Ngrid = 10, kernel = "linear", kparam = 1, scale = FALSE, criterion = "0-1", cv_type = "original", interaction = FALSE,
                    gd_scale = TRUE, gamma = 0.5, optModel = FALSE, nCores = 1, ...)
{
  # Find a optimal lambda in first step
  initial_fit = Kfold_msvm(x = x, y = y, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds, lambda_seq = lambda_seq,
                           gamma = gamma, kernel = kernel, kparam = kparam, scale = scale, criterion = criterion,
                           gd = TRUE, gd_scale = gd_scale, optModel = FALSE, nCores = nCores, ...)
  
  
  # Find a relevant variable for optimal lambda in second step
  select_fit = threshold_fun(initial_fit, thresh_Ngrid = thresh_Ngrid, cv_type = cv_type, criterion = criterion,
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
    opt_sigma = kernlab::sigest(y ~ selected_x, frac = 1, scaled = FALSE)[3]
    final_fit = Kfold_msvm(x = selected_x, y = y, valid_x = selected_valid_x, valid_y = valid_y,
                           nfolds = nfolds, lambda_seq = lambda_seq, gamma = gamma, kernel = kernel, kparam = opt_sigma,
                           scale = scale, criterion = criterion, gd = FALSE, gd_scale = gd_scale, optModel = TRUE, nCores = nCores, ...)
    out$opt_model = final_fit$opt_model
    out$valid_err = final_fit$valid_err
    out$opt_valid_err = final_fit$opt_valid_err
    
  }
  
  out$cv_type = select_fit$cv_type
  out$call = call
  class(out) = "GBFSMSVM"  
  return(out)  
}
  
  
  
  
  
  