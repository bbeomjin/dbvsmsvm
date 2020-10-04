# setwd(r"(C:\Users\Beom\desktop)")
# rm(list = ls())
# gc()
# 
# require(ramsvm)
# require(kernlab)
# require(caret)
# require(dplyr)
# require(GBFSMSVM)
# 
# source("C:\\Users\\Beom\\OneDrive - 서울시립대학교\\feature selection for kernel method\\MSVM_FS\\msvm_fs_20200801/sim_gen.R")
# source("C:\\Users\\Beom\\OneDrive - 서울시립대학교\\feature selection for kernel method\\MSVM_FS\\msvm_fs_20200801/SMSVM_1.2.1.R")
# source("C:\\Users\\Beom\\Documents\\GitHub\\GBFSMSVM\\R\\kfold_func.R")
# source("C:\\Users\\Beom\\Documents\\GitHub\\GBFSMSVM\\R\\SRAMSVM_solve.R")
# source("C:\\Users\\Beom\\Documents\\GitHub\\GBFSMSVM\\R\\sub_funcs.R")
# source("C:\\Users\\Beom\\Documents\\GitHub\\GBFSMSVM\\R\\thresholding_func.R")
# 
# 
# param_set = expand.grid(n = c(100, 200), p = c(10, 20))
# iter = 100
# 
# msvm_list_linear = msvm_list_rbf = msvm_list_spline = 
#   msvm_list_linear_thr = msvm_list_rbf_thr = msvm_list_spline_thr = 
#   msvm_list_linear_thr_hinge = msvm_list_rbf_thr_hinge = msvm_list_spline_thr_hinge = vector(mode = "list", length = nrow(param_set))
# cosso_list_linear = cosso_list_rbf = cosso_list_spline = 
#   cosso_list_linear_hinge = cosso_list_rbf_hinge = cosso_list_spline_hinge = vector(mode = "list", length = nrow(param_set))
# msvm_time_list_linear = msvm_time_list_rbf = msvm_time_list_spline = 
#   msvm_time_list_linear_thr = msvm_time_list_rbf_thr = msvm_time_list_spline_thr = 
#   msvm_time_list_linear_thr_hinge = msvm_time_list_rbf_thr_hinge = msvm_time_list_spline_thr_hinge = vector(mode = "list", length = nrow(param_set))
# cosso_time_list_linear = cosso_time_list_rbf = cosso_time_list_spline = 
#   cosso_time_list_linear_hinge = cosso_time_list_rbf_hinge = cosso_time_list_spline_hinge = vector(mode = "list", length = nrow(param_set))
# 
# sram_list_linear = sram_list_rbf = sram_list_linear_hinge = sram_list_rbf_hinge = vector(mode = "list", length = nrow(param_set))
# sram_time_list_linear = sram_time_list_rbf = vector(mode = "list", length = nrow(param_set))
# 
# sram_perf_linear = sram_perf_linear_hinge = sram_perf_rbf = sram_perf_rbf_hinge = data.frame(matrix(NA, nrow = iter, ncol = 3))   
# sram_time_vec_linear = sram_time_vec_rbf = rep(0, iter)
# cv_err_sram_linear = cv_err_sram_rbf = rep(0, iter)
# cv_err_sram_linear_list = cv_err_sram_rbf_list = vector(mode = "list", length = nrow(param_set))
# 
# msvm_perf_linear = msvm_perf_rbf = msvm_perf_spline = data.frame(matrix(NA, nrow = iter, ncol = 3))
# msvm_perf_linear_thr = msvm_perf_rbf_thr = msvm_perf_spline_thr = data.frame(matrix(NA, nrow = iter, ncol = 3))
# msvm_perf_linear_thr_hinge = msvm_perf_rbf_thr_hinge = msvm_perf_spline_thr_hinge = data.frame(matrix(NA, nrow = iter, ncol = 3))
# cosso_perf_linear = cosso_perf_rbf = cosso_perf_spline = data.frame(matrix(NA, nrow = iter, ncol = 3))
# cosso_perf_linear_hinge = cosso_perf_rbf_hinge = cosso_perf_spline_hinge = data.frame(matrix(NA, nrow = iter, ncol = 3))
# msvm_time_vec_linear = msvm_time_vec_rbf = msvm_time_vec_spline = 
#   msvm_time_vec_linear_thr = msvm_time_vec_rbf_thr = msvm_time_vec_spline_thr = 
#   msvm_time_vec_linear_thr_hinge = msvm_time_vec_rbf_thr_hinge = msvm_time_vec_spline_thr_hinge = rep(0, iter)
# cosso_time_vec_linear = cosso_time_vec_rbf = cosso_time_vec_spline = cosso_time_vec_linear_hinge = cosso_time_vec_rbf_hinge = cosso_time_vec_spline_hinge = rep(0, iter)
# 
# cv_err_thresh_linear_list = cv_err_thresh_rbf_list = cv_err_cosso_linear_list = cv_err_cosso_rbf_list = vector(mode = "list", length = nrow(param_set))
# cv_err_thresh_linear = cv_err_thresh_rbf = cv_err_cosso_linear = cv_err_cosso_rbf = rep(0, iter)
# 
# folds = 5
# i = 1
# j = 1
# 
# cat(j, "th-setting \n")
# n = param_set$n[j]
# # p = param_set$p[j]
# n = 300
# p = 5
# 
# cat(i, "th-iteration \n")
# # cosso -> cosso2
# dat = sim_gen(n = n, p = p, class = 3, seed = i, type = "neuralnet2")
# test_dat = sim_gen(n = n, p = p, class = 3, seed = i + 1, type = "neuralnet2")
# X = dat$x
# X = scale(X)
# test_X = scale(test_dat$x)
# y = as.integer(dat$y)
# test_y = as.integer(dat$y)
# 
# true = factor(dat$true, levels = c(1, 0))
# sigma = sigest(y ~ X, frac = 1, scale = FALSE)[3]
# 
# 
# kfold_fit = Kfold_msvm(X, y, nfold = 5, lambda_seq = 2^seq(-20, 0, length.out = 100), gamma = 0.5, kernel = "radial", kparam = sigma, 
#                        scale = FALSE, criterion = "0-1")
# 
# 
# gd_fit = threshold_fun.GBFSMSVM(kfold_fit, thresh_Ngrid = 100, cv_type = "osr", criterion = "0-1", gd_scale = TRUE, interaction = TRUE)
# gd_fit$selected
# gd_fit$int_selected
# gd_fit$int_valid_err
# gd_fit$gd_interaction
# gd_fit$opt_thresh_int
# 
# K = kernlab::kernelMatrix(rbfdot(sigma), X)
# int_K = interaction_kernel(X, X, kernel = list(type = "radial", par = sigma))
# test_K = kernlab::kernelMatrix(rbfdot(sigma), X, test_X)
# test_int_K = interaction_kernel(X, test_X, kernel = list(type = "radial", par = sigma))
# 
# KK = K - int_K
# test_KK = test_K - test_int_K
# 
# fit1 = GBFSMSVM:::SRAMSVM_solve(y = y, K = KK, lambda = 1e-10, kernel = "radial", kparam = sigma, gamma = 0.5, maxit = 1000, epsilon = 1e-6)
# table(y, GBFSMSVM:::predict.sramsvm(fit1, newK = KK)[[1]][[1]])
# pred_y = GBFSMSVM:::predict.sramsvm(fit1, newK = test_KK)
# tt1 = table(test_y, pred_y[[1]][[1]])
# sum(diag(tt1)) / sum(tt1)
# 
# 
# fit2 = SRAMSVM_solve(y = y, K = K, lambda = 1e-10, kernel = "radial", kparam = sigma, gamma = 0.5, maxit = 1000, epsilon = 1e-6)
# pred_y2 = GBFSMSVM:::predict.sramsvm(fit2, newK = test_K)
# tt2 = table(test_y, pred_y2[[1]][[1]])
# sum(diag(tt2)) / sum(tt2)
# 
# 
# pred_y[[2]][1:10, ]
# pred_y2[[2]][1:10, ]
# 
# table(test_y, apply(pred_y[[2]] - pred_y2[[2]], 1, which.max))
# 
# GBFSMSVM:::gl_fun_interaction(fit2$beta[[1]], X, y, scale = FALSE, kernel = "radial", kernel_par = list(sigma))
# 
# alpha = fit2$beta[[1]]
# 
# gl_fun_interaction = function (alpha, X, y, scale = TRUE, kernel = c("linear", "poly", "radial"), kernel_par = list(), active_set = NULL) 
# {
#   n = length(y)
#   if (kernel == "linear") {
#     if (scale) 
#       scale_const = sum(drop(crossprod(alpha, X))^2)
#   }
#   if (kernel == "poly" | kernel == "poly2") {
#     if (scale) {
#       poly = do.call(polydot, kernel_par)
#       scale_const = drop(crossprod(alpha, kernelMatrix(poly, 
#                                                        X)) %*% alpha)
#     }
#   }
#   if (kernel == "radial" | kernel == "radial2") {
#     names(kernel_par) = "sigma"
#     if (scale) {
#       rbf = do.call(rbfdot, kernel_par)
#       scale_const = drop(crossprod(alpha, kernelMatrix(rbf, 
#                                                        X) + 1) %*% alpha)
#       scale_const = diag(scale_const)
#     }
#   }
#   # ddkernel = switch(kernel, linear = ddlinear, poly = ddpoly, 
#   #                   poly2 = ddpoly, radial = ddrbf, radial2 = ddrbf)
#   ddkernel = GBFSMSVM:::ddrbf
#   comb_set = combn(active_set, 2)
#   gd = lapply(1:ncol(comb_set), function(k) {
#     gd_res = sapply(1:n, FUN = function(i) {
#       return(crossprod(alpha, ddkernel(X, X[i, ], kernel_par, 
#                                        comb_set[1, k], comb_set[2, k])))
#     })
#     return(gd_res)
#   })
#   if (scale) {
#     res = sapply(1:ncol(comb_set), function(x) rowMeans(gd[[x]]^2)/scale_const)
#   }
#   else {
#     res = sapply(1:ncol(comb_set), function(x) rowMeans(gd[[x]]^2))
#   }
#   return(res)
# }
# 
# gd[[1]]
# 
# a = GBFSMSVM:::gradient(fit2$beta[[1]], X, y, scale = FALSE, kernel = "radial", kparam = list(sigma = sigma))
# 
# crossprod(fit2$beta[[1]], GBFSMSVM:::drbf(X, X[1, ], list(sigma = sigma)))
# 
# rbf_cal = function(x, y, sigma) {
#   const = exp(-sigma * (crossprod(x) + crossprod(y)))
#   # kern = 1 + 2 * sigma * crossprod(x, y)
#   kern = exp(2 * sigma * x[1] * y[1]) + exp(2 * sigma * x[2] * y[2])
#   return(const * kern)
# }
# 
# rbf_cal = function(x, y, sigma) {
#   
#   # kern = 1 + 2 * sigma * crossprod(x, y)
#   kern = (exp(2 * sigma * x[1] * y[1]) + exp(2 * sigma * x[2] * y[2]))
#   return(kern)
# }
# 
# rbf_cal2 = function(x, y, sigma) {
#   
#   # kern = 1 + 2 * sigma * crossprod(x, y)
#   kern = exp(-sigma * (crossprod(x) + crossprod(y)) + 2 * sigma * crossprod(x, y))
#   return(kern)
# }
# 
# 
# X = matrix(c(1, 2, 4, 7), nrow = 2)
# aa = kernlab::kernelMatrix(rbfdot(sigma = 1), X)
# kernlab::kernelMatrix(rbfdot(sigma = 1), X[, 1])
# 
# 
# GBFSMSVM:::make_anovaKernel(X[, 1, drop = FALSE], X[, 1, drop = FALSE], kernel = list(type = "radial", par = sigma))
# GBFSMSVM:::make_anovaKernel(X[, 2, drop = FALSE], X[, 2, drop = FALSE], kernel = list(type = "radial", par = sigma))
# rbf_cal(X[1, ], X[2, ], sigma = 1)
# rbf_cal2(X[1, ], X[2, ], sigma = 1)
