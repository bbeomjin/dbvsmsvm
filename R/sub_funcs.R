require(igraph)
main_kernel = function(x, u, kernel)
{
  x = as.matrix(x)
  u = as.matrix(u)
  if (kernel$type == "linear")
    K = (x %*% t(u))
  if (kernel$type == "poly")
    K = (1 + x %*% t(u))^kernel$par
  if (kernel$type == "radial" | kernel$type == "radial2")
  {
    a = as.matrix(rowSums(x^2))
    b = as.matrix(rowSums(u^2))
    one.a = matrix(1, ncol = length(b))   
    one.b = matrix(1, ncol = length(a))
    K1 = one.a %x% a
    K2 = x %*% t(u)
    K3 = t(one.b %x% b)
    K = exp(-(K1 - 2 * K2 + K3) * (kernel$par))
  }
  return(K)
}

spline_kernel = function(x, u)
{
  x = as.matrix(x)
  u = as.matrix(u)
  K1x = (x - 1 / 2)
  K1u = (u - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2u = (K1u^2 - 1 / 12) / 2
  ax = x%x%matrix(1, 1, nrow(u)) 
  au = u%x%matrix(1, 1, nrow(x))
  b = abs(ax - t(au))
  K1 = K1x%x%t(K1u)
  K2 = K2x%x%t(K2u) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  list(K1 = K1, K2 = K2)
}


make_anovaKernel = function(x, u, kernel)
{
  if (!is.matrix(x))  # degenerate case: x is a row vector  
  { x = t(as.matrix(x))}
  else { x = as.matrix(x)}
  
  u = as.matrix(u)
  dimx = ncol(x)
  
  # calculate anova kernels for main effects
  if(kernel$type == "spline")
  {
    # assign the number of anova kernels
    numK = 2*dimx
    # list of kernel matrices
    anova_kernel = vector(mode="list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode="list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }  
  }
  else if (kernel$type == 'spline2')
  {
    numK = (2*dimx) + (2*dimx*(2*dimx-1)/2 - dimx)
    anova_kernel = vector(mode="list", numK)
    kernelCoord = vector(mode="list", numK)
    index = 0
    # main effects
    for(d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }  
    # two-way interactions
    for (i in 1:(dimx-1))
    {
      for (j in (i+1):dimx)
      {
        index = index + 1
        A.linear = as.matrix(anova_kernel[[2*i-1]])
        A.smooth = as.matrix(anova_kernel[[2*i]])
        B.linear = as.matrix(anova_kernel[[2*j-1]])
        B.smooth = as.matrix(anova_kernel[[2*j]])            
        anova_kernel[[index]] = A.linear*B.linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep="")
        index = index + 1
        anova_kernel[[index]] = A.linear*B.smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep="")
        index = index + 1
        anova_kernel[[index]] = A.smooth*B.linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep="")
        index = index + 1
        anova_kernel[[index]] = A.smooth*B.smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep="")
      }
    }
  }
  else if (kernel$type == "spline-t")
  {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
  }
  else if (kernel$type == 'spline-t2')
  {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[, d])
      B = as.matrix(u[, d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
    for (i in 1:(dimx - 1))
    {
      for (j in (i + 1):dimx)
      {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep="")
      }
    }
  } else if (kernel$type == "radial2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      anova_kernel[[index]] = main_kernel(A, B, kernel)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
    for (i in 1:(dimx - 1))
    {
      for (j in (i + 1):dimx)
      {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep="")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx)
    {
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      anova_kernel[[d]] = main_kernel(A, B, kernel)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel)
}

combine_kernel = function(anova_kernel, theta = rep(1, anova_kernel$numK))
{
  K = 0
  for (d in 1:anova_kernel$numK)
  {
    K = (K + theta[d] * anova_kernel$K[[d]])
  }
  return(K)
}


find_theta = function(y, anova_kernel, gamma = 0.5, cmat, bvec, lambda, lambda_theta) 
{
  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }
  # standard LP form :
  # min a^T x , subject to A1x <= a1
  n_data = length(y)
  n_class = length(unique(y))
  bvec = as.matrix(bvec)
  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, as.double(n_class), nobs = n_data, y_train = y)
  
  # calculate the 'a' matrix 
  a_tmp = matrix(gamma / n_data, nrow = n_data, ncol = n_class)
  a_tmp[cbind(1:n_data, y)] = (1 - gamma) / n_data
  a = matrix(a_tmp, ncol=1)
  
  # initialize M 
  M = matrix(rep(0, anova_kernel$numK), ncol=1)
  # calculate M
  for (d in 1:anova_kernel$numK)
  {
    for (j in 1:(n_class - 1))
    {
      M[d] = (M[d] + t(cmat[, j]) %*% anova_kernel$K[[d]] %*% cmat[, j])
    }
    M[d] = (lambda / 2 * M[d] + (lambda_theta))
  }
  a = rbind(a, M)
  
  Y_code = Y_matrix_gen(n_class, as.double(n_class), nobs = n_class, y_train = 1:n_class)
  
  # calculate N matrix
  for (d in 1:anova_kernel$numK)
  {
    K = anova_kernel$K[[d]]
    for (j in 1:n_class)
    {
      if(j == 1)
      {
        temp_N = matrix(rowSums(K %*% cmat * matrix(Y_code[j, ], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE)))
      } 
      else 
      {
        temp_N = rbind(temp_N, matrix(rowSums(K %*% cmat * matrix(Y_code[j, ], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE))))
      }
    }
    if(d == 1)
    {
      N = temp_N
    } 
    else 
    {
      N = cbind(N, temp_N)
    }
  }
  
  sign_mat = matrix(1, n_data, n_class)
  sign_mat[cbind(1:n_data, y)] = -1
  
  # constraints
  
  # A matrix
  I.nonzeroIndex = diag(1, n_data * n_class)
  N.nonzeroIndex = as.matrix(N) * as.vector(sign_mat)
  A.theta = cbind(matrix(0, anova_kernel$numK, n_data * n_class),
                  diag(-1, anova_kernel$numK))
  A_ineq = rbind(cbind(I.nonzeroIndex, -N.nonzeroIndex), A.theta)
  
  bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * bvec[x]))
  
  bb_yi = (n_class - 1) - bb
  
  bb_j = 1 + matrix(rowSums(Y_code * matrix(bvec, nrow = nrow(Y_code), ncol = ncol(Y_code), byrow = TRUE)), nrow = n_data, ncol = n_class, byrow = TRUE)
  bb_j[cbind(1:n_data, y)] = bb_yi
  
  b_nonzeroIndex = as.vector(bb_j)
  b_ineq = rbind(as.matrix(b_nonzeroIndex), matrix(-1, anova_kernel$numK, 1))
  
  # constraint directions
  const_dir = matrix(rep(">=", nrow(matrix(b_ineq))))
  
  # find solution by LP
  
  lp = lp("min", objective.in = a, const.mat = A_ineq, const.dir = const_dir,
          const.rhs = b_ineq)$solution
  
  # find the theta vector only from the solution
  theta = cbind(matrix(0, anova_kernel$numK, n_data * n_class),
                diag(1, anova_kernel$numK)) %*% matrix(lp, ncol=1)
  return(as.vector(theta))
}

find_intercept = function(y, K, gamma = 0.5, cmat, lambda) 
{
  # standard LP form :
  # min a^T x , subject to A1x <= a1
  n_data = length(y)
  n_class = length(unique(y))
  # bvec = as.matrix(bvec)
  
  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, as.double(n_class), nobs = n_data, y_train = y)
  
  # calculate the 'a' matrix 
  # a_tmp = matrix((1 - gamma) / n_data, nrow = n_data, ncol = n_class)
  # a_tmp[cbind(1:n_data, y)] = gamma / n_data
  a_tmp = matrix((1 - gamma), nrow = n_data, ncol = n_class)
  a_tmp[cbind(1:n_data, y)] = gamma
  a = matrix(a_tmp, ncol = 1)
  a = rbind(a, matrix(rep(0, n_class - 1)))
  
  Y_code = Y_matrix_gen(n_class, as.double(n_class), nobs = n_class, y_train = 1:n_class)
  
  sign_mat = matrix(1, n_data, n_class)
  sign_mat[cbind(1:n_data, y)] = -1
  sign_vec = as.vector(sign_mat)
  # constraints
  
  for (j in 1:n_class)
  {
    if(j == 1)
    {
      temp_N = matrix(rowSums(K %*% cmat * matrix(Y_code[j, ], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE)))
    } 
    else 
    {
      temp_N = rbind(temp_N, matrix(rowSums(K %*% cmat * matrix(Y_code[j, ], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE))))
    }
    N = temp_N
  }
  
  # A matrix
  
  A_ineq_h = cbind(diag(1, n_data * n_class), matrix(0, ncol = n_class - 1, nrow = n_data * n_class))
  A_ineq_c = cbind(diag(1, n_data * n_class), -sign_vec *  matrix(rep(Y_code, each = n_data), nrow = n_data * n_class, n_class - 1))
  A_ineq = rbind(A_ineq_h, A_ineq_c)
  
  # I.nonzeroIndex = diag(1, n_data * n_class)
  # N.nonzeroIndex = as.matrix(N) * as.vector(sign_mat)
  # A.theta = cbind(matrix(0, 1, n_data * n_class),
  #                 diag(-1, 1))
  # A_ineq = rbind(cbind(I.nonzeroIndex, -N.nonzeroIndex), A.theta)
  
  # bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * bvec[x]))
  
  # bb_yi = (n_class - 1) - bb
  
  # bb_j = 1 + matrix(rowSums(Y_code * matrix(bvec, nrow = nrow(Y_code), ncol = ncol(Y_code), byrow = TRUE)), nrow = n_data, ncol = n_class, byrow = TRUE)
  # bb_j[cbind(1:n_data, y)] = bb_yi
  
  bb = matrix(1, n_data, n_class)
  bb[cbind(1:n_data, y)] = n_class - 1
  b_ineq_c = as.vector(bb) +  sign_vec * N
  b_ineq_h = matrix(rep(-1e-10, n_data * n_class))
  b_ineq = rbind(b_ineq_h, b_ineq_c)
  
  # constraint directions
  # const_dir = matrix(rep(">=", n_data * n_class + (n_class - 1)))
  const_dir = matrix(rep(">=", 2 * (n_data * n_class)))
  # find solution by LP
  
  lp = lp("min", objective.in = a, const.mat = A_ineq, const.dir = const_dir,
          const.rhs = b_ineq)$solution
  
  # find the theta vector only from the solution
  beta0 = lp[(length(lp) - n_class + 2):length(lp)]
  return(as.vector(beta0))
}



XI_gen = function(k, kd) {
  
  tempA = - (1.0 + sqrt(kd)) / ((kd - 1.0)^(1.5))
  tempB = tempA + sqrt(kd / (kd - 1.0))
  
  XI = matrix(data = tempA, nrow = k - 1L, ncol = k)
  
  XI[, 1L] = 1.0 / sqrt(kd - 1.0)
  
  for( ii in 2L:k ) XI[ii - 1L, ii] = tempB
  
  return(XI)
}

Y_matrix_gen = function(k, kd, nobs, y_train) {
  
  XI = XI_gen(k = k, kd = kd)
  
  Y_matrix = matrix(data = 0.0, nrow = nobs, ncol = k - 1L)
  
  for( ii in 1L:nobs ) Y_matrix[ii, ] = XI[, y_train[ii]]
  
  return( Y_matrix )
  
}


ramsvm_hinge = function(y, fit, k, gamma = 0.5)
{
  n = length(y)
  hinge_j = (k - 1) - fit[cbind(1:length(y), y)]
  hinge_j[hinge_j < 0] = 0
  
  hinge_mj = sapply(1:length(y), function(x) fit[x, -y[x]]) + 1
  hinge_mj[hinge_mj < 0] = 0
  hinge_loss = 1 / n * sum((1 - gamma) * colSums(hinge_mj) + gamma * hinge_j)
  return(hinge_loss)
}


kernelMatrix_spline = function(X, Y) {
  K = 0
  spline_kernel = function(x, u)
  {
    x = as.matrix(x)
    u = as.matrix(u)
    K1x = (x - 1/2)
    K1u = (u - 1/2)
    K2x = (K1x^2 - 1/12) / 2
    K2u = (K1u^2 - 1/12) / 2
    ax = x %x% matrix(1, 1, nrow(u)) 
    au = u %x% matrix(1, 1, nrow(x))
    b = abs(ax - t(au))
    K1 = K1x %x% t(K1u)
    K2 = K2x %x% t(K2u) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
    return(list(K1 = K1, K2 = K2))
  }
  
  for(d in 1:p)
  {
    K_temp = spline_kernel(as.matrix(X[, d]), as.matrix(Y[, d]))
    K = K + K_temp$K1 + K_temp$K2
  }
  return(K)
}


drbf = function(X, xp, kernel_par = list(sigma = NULL))
{
  gamma = kernel_par$sigma
  # diff_mat = sweep(X, xp, MARGIN = 2)
  # diff_mat = t(t(X) - xp)
  np = dim(X)
  diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
  # diff_mat = X - xp[col(X)]
  # tmp = exp(-rowSums(diff_mat^2) / (2 * sigma)) * (diff_mat / sigma)
  tmp = exp(-rowSums(diff_mat^2) * gamma) * (2 * gamma * diff_mat)
  # dgd = drop(crossprod(alpha, tmp))
  return(tmp)
}

# drbf_mul = function(alpha, X, xp, kernel_par = list(sigma = NULL))
# {
#   gamma = kernel_par$sigma
#   np = dim(X)
#   diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
#   dgd = crossprod(alpha, exp(-rowSums(diff_mat^2) * gamma) * 2 * gamma * diff_mat)
#   return(dgd)
# }


ddrbf = function(X, xp, kernel_par = list(sigma = NULL), ind1, ind2)
{
  gamma = kernel_par$sigma
  np = dim(X)
  diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
  tmp = exp(-rowSums(diff_mat^2) * gamma) * (2 * gamma * diff_mat[, ind1]) * (2 * gamma * diff_mat[, ind2])
  # dgd = drop(crossprod(alpha, tmp))
  return(tmp)
}


dlinear = function(X, xp, kernel_par = list())
{
  tmp = X
  return(tmp)
}

# dlinear_mul = function(alpha, X, xp, kernel_par = list())
# {
#   return(drop(crossprod(alpha, X)))
# }

dpoly = function(alpha, X, xp, kernel_par = list(degree = NULL))
{
  # degree = kernel_par$degree
  # scale = kernel_par$scale
  # offset = kernel_par$offset
  degree = kernel_par[[1]]
  scale = 1
  offset = 0
  tmp = degree * (scale * drop(X %*% xp) + offset)^{degree - 1} * scale * X
  return(tmp)
}

# dpoly_mul = function(alpha, X, xp, kernel_par = list(degree = NULL))
# {
#   degree = kernel_par[[1]]
#   scale = 1
#   offset = 0
#   tmp = degree * (scale * drop(X %*% xp) + offset)^{degree - 1} * scale * X
#   return(drop(crossprod(alpha, tmp)))
# }

ddspline = function(x, y)
{
  m1 = (x - (1 / 2)) + (1 / 2) * ((x - (1 / 2))^2 - (1 / 12)) * (y - (1 / 2))
  m2 = ifelse(x > y, -(1 / 6) * (x - y - (1 / 2))^3 + (1 / 24) * (x - y - (1 / 2)), 
              (4 * (y - x - (1 / 2))^3 - (y - x - (1 / 2))) / 24)
  return(m1 - m2)
}

# dspline_mul = function(alpha, X, xp, kernel_par = list())
# {
#   res_mat = matrix(0, nrow = NROW(X), ncol = NCOL(X))
#   for (k in 1:NCOL(X)) {
#     res_mat[, k] = ddspline(X[, k], xp[k])
#   }
#   dgd = drop(crossprod(alpha, res_mat))
#   return(dgd)
# }

gradient = function(alpha, x, y, scale = TRUE, kernel = c("linear", "poly", "radial", "spline", "anova_radial"),
                  kparam = list())
{
  n = length(y)
  k = length(unique(y))
  
  if (kernel == "linear") {
    if (scale) {scale_const = sapply(1:NCOL(alpha), FUN = function(i) sum(crossprod(alpha[, i], x)^2))}
  }
  
  if (kernel == "poly") {
    if (scale) {
      poly = do.call(polydot, kparam)
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix(poly, x)) %*% alpha[, i]))
    }
  }
  
  if ((kernel == "radial") | (kernel == "anova_radial")) {
    names(kparam) = "sigma"
    if (scale) {
      rbf = do.call(rbfdot, kparam)
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix(rbf, x) + 1) %*% alpha[, i]))
      # scale_const = drop(t(alpha) %*% kernelMatrix(rbf, X) %*% alpha)
      # system.time((a = sum(sapply(1:n, FUN = function(i) K_rbf(alpha, X, X[i, ], sigma = sigma)) * alpha)))
    }
  }
  
  if (kernel == "spline") {
    if (scale) {
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix_spline(x, x)) %*% alpha[, i]))
    }
  }
  
  dkernel = switch(kernel,
                   linear = dlinear,
                   poly = dpoly,
                   radial = drbf,
                   spline = dspline,
				   anova_radial = drbf)
  
  
  W_mat = XI_gen(k, as.double(k))
  
  grad_mat = 0
  for (i in 1:n) {
    dK_sq = crossprod(alpha, dkernel(x, x[i, ], kparam))^2
    # gd_mat = gd_mat + crossprod(W_mat, dK)^2 / n
    grad_mat = grad_mat + dK_sq / n
  }
  # browser()
  # print(dim(dK))
  
  # print(dim(gd_mat))
  # print(scaler)
  # gd = colSums(gd_mat) / scaler
  # res = gd
  
  if (scale) {
    scaler = sum(drop(crossprod(W_mat, sqrt(scale_const)))^2)
    # res = colSums(gd_mat) / scaler
    res = colSums(grad_mat) / k / scaler
  } else {
    res = colSums(grad_mat) / k
  }
  
  # if (scale) {res = gd / scale_const} else {res = gd}
  # if (scale) {res = rowMeans(gd^2) / scale_const} else {res = rowMeans(gd^2)}
  return(res)
}


gradient_interaction = function(alpha, x, y, scale = TRUE, kernel = c("linear", "poly", "radial"),
                              kparam = list(), active_set = NULL)
{
  n = length(y)
  k = length(unique(y))
  
  if (kernel == "linear") {
    if (scale) {scale_const = sapply(1:NCOL(alpha), FUN = function(i) sum(crossprod(alpha[, i], x)^2))}
  }
  
  if (kernel == "poly") {
    if (scale) {
      poly = do.call(polydot, kparam)
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix(poly, x)) %*% alpha[, i]))
    }
  }
  
  if ((kernel == "radial") | (kernel == "anova_radial")) {
    names(kparam) = "sigma"
    if (scale) {
      rbf = do.call(rbfdot, kparam)
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix(rbf, x)) %*% alpha[, i]))
      # scale_const = drop(t(alpha) %*% kernelMatrix(rbf, X) %*% alpha)
      # system.time((a = sum(sapply(1:n, FUN = function(i) K_rbf(alpha, X, X[i, ], sigma = sigma)) * alpha)))
    }
  }
  
  if (kernel == "spline") {
    if (scale) {
      scale_const = sapply(1:NCOL(alpha), FUN = function(i) drop(crossprod(alpha[, i], kernelMatrix_spline(x, x)) %*% alpha[, i]))
    }
  }
  
  ddkernel = switch(kernel,
                    linear = ddlinear,
                    poly = ddpoly,
                    poly2 = ddpoly,
                    radial = ddrbf,
                    radial2 = ddrbf)
  
  W_mat = XI_gen(k, as.double(k))
  
  comb_set = combn(active_set, 2)
  
  gd = sapply(1:ncol(comb_set), function(k) {
    gd_res = rowMeans(sapply(1:n, FUN = function(i) {return(crossprod(alpha, ddkernel(x, x[i, ], kparam, comb_set[1, k], comb_set[2, k]))^2)}))
    return(gd_res)
  })
  if (scale) {
    # res = sapply(gd, function(x) return(mean(x^2) / scale_const))
    res = colSums(gd) / k
  } else {
    res = colSums(gd) / k
  }
  
  # if (scale) {res = rowMeans(gd^2) / scale_const} else {res = rowMeans(gd^2)}
  return(res)
}

data_split = function(y, fold, k = max(y), seed = length(y))
{
  # k: the number of classes
  n_data = length(y)
  class_size = table(y)
  ran = rep(0, n_data) 
  if ( (min(class_size) < fold) & (fold != n_data) )
  {
    warning(' The given fold is bigger than the smallest class size. \n Only a fold size smaller than the minimum class size \n or the same as the sample size (LOOCV) is supported.\n')
    return(NULL)
  }
  
  if ( min(class_size) >= fold )
  {
    set.seed(seed)
    for (j in 1:k)
    {  
      ran[y == j] = ceiling(sample(class_size[j]) / (class_size[j] + 1) * fold) 
    }
  }
  else if ( fold == n_data)
  {
    ran = 1:n_data
  }
  return(ran)
}

kernelMat = function(x, y, kernel = "radial", kparam = 1.0) {
  
  if (NCOL(x) == 0) {
    x = matrix(1, nrow = nrow(x), ncol = 1)
  }
  
  if (NCOL(y) == 0) {
    y = matrix(1, nrow = nrow(y), ncol = 1)
  }
  
  if( kernel == "poly" ) {
    obj = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "radial" | kernel == "radial2") {
    # normx = drop((x^2) %*% rep(1.0, ncol(x)))
    # normy = drop((y^2) %*% rep(1.0, ncol(y)))
    # temp = x %*% t(y)
    # temp = (-2.0 * temp + normx) + outer(rep(1.0, nrow(x)), normy, "*")
    # obj = exp(-temp * kparam)
    obj = kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "spline") {
    K = 0
    spline_kernel = function(x, u)
    {
      x = as.matrix(x)
      u = as.matrix(u)
      K1x = (x - 1/2)
      K1u = (u - 1/2)
      K2x = (K1x^2 - 1/12) / 2
      K2u = (K1u^2 - 1/12) / 2
      ax = x %x% matrix(1, 1, nrow(u)) 
      au = u %x% matrix(1, 1, nrow(x))
      b = abs(ax - t(au))
      K1 = K1x %x% t(K1u)
      K2 = K2x %x% t(K2u) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
      return(list(K1 = K1, K2 = K2))
    }
    for(d in 1:p)
    {
      K_temp = spline_kernel(as.matrix(x[, d]), as.matrix(y[, d]))
      K = K + K_temp$K1 + K_temp$K2
    }
    obj = K
  } else if (kernel == "linear") {
    obj = tcrossprod(x, y)
  } else if (kernel == "anova_radial") {
      K = 0
      for (d in 1:NCOL(x))
      {
        A = as.matrix(x[,d])
        B = as.matrix(y[,d])
        K_temp = main_kernel(A, B, kernel = list(type = "radial", par = kparam))
        K = K + K_temp
      }
	obj = K
  } else {
    obj = NULL
  }
  return(obj)
}

interaction_svmfs = function(main_effect, interaction) 
{
  if (sum(interaction) != 0) {
    p = length(main_effect)
    comb_mat = combn(1:p, 2)
    ind_mat = comb_mat[, interaction == 1, drop = FALSE]
    for (i in 1:nrow(ind_mat)) {
      ind = ind_mat[i, ]
      main_effect[ind] = 1
    }
  }
  res = c(main_effect, interaction)
  return(res)
}


interaction_graph = function(comb, p, min = 3)
{
  int_mat = matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p - 1) / 2)
  int_mat[t(comb)] = 1
  int_mat = t(int_mat) + int_mat
  g = graph_from_adjacency_matrix(int_mat, mode = "undirected")
  cliques_list = max_cliques(g, min = 3)
  return(cliques_list)
}

interaction_kernel = function(x, u, kernel, active_set, interaction_set, clique_list)
{
  if (!is.matrix(x)) {
    x = as.matrix(x)
  } else {
    x = as.matrix(x)
  }
  u = as.matrix(u)
  dimx = ncol(x)
  
  scaled_kernel = function(x, u, kernel, active_set, index)
  {
    X1 = matrix(rowSums(x[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U1 = matrix(rowSums(u[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    X2 = matrix(rowSums(x[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U2 = matrix(rowSums(u[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    K = exp(-kernel$par * ((X1 + U1) - (X2 + U2)))
    K1 = exp(-kernel$par * (X1 + U1))
    K2 = exp(-kernel$par * (X2 + U2))
    K_mat = main_kernel(x[, index, drop = FALSE], u[, index, drop = FALSE], kernel)
    res = K * K_mat - K1
    return(list(res = res, K = K, K1 = K1, K2 = K2))
  }
  
  # numK = dimx * (dimx - 1) / 2
  # anova_kernel_temp = vector(mode = "list", dimx)
  main_effects = vector(mode = "list", dimx)
  high_order_kernel = vector(mode = "list", length(clique_list))
  # const_term = vector(mode = "list", dimx)
  # kernelCoord = vector(mode = "list", numK)
  
  
  for (j in 1:length(active_set)) {
    temp_kernel = scaled_kernel(x, u, kernel, active_set, active_set[j])
    main_effects[[active_set[j]]] = temp_kernel$res
    # const_term[[active_set[j]]] = temp_kernel[-1]
  }
  
  
  if (length(clique_list) != 0) {
    for (d in 1:length(clique_list)) {
      ind = sort(as.vector(clique_list[[d]]))
      # scale_const = scaler(x, u, kernel, active_set, ind)
      # anova_kernel[[d]] = scale_const$K * (main_kernel(x[, ind, drop = FALSE], u[, ind, drop = FALSE], kernel))
      clique_kernel = scaled_kernel(x, u, kernel, active_set, ind)$res
      temp_comb = combn(ind, 2)
      interaction_effects = lapply(1:ncol(temp_comb), FUN = function(i) {
        ind = temp_comb[, i]
        return(((main_effects[[ind[1]]]) * (main_effects[[ind[2]]])) / temp_kernel$K1)
      })
      sum_interaction = Reduce("+", interaction_effects)
      sum_main = Reduce("+", main_effects[ind])
      high_order_kernel[[d]] = clique_kernel - sum_main - sum_interaction
    }
  }
  
  interaction_kernel = list(0)
  if (ncol(interaction_set) != 0) {
    interaction_kernel = lapply(1:ncol(interaction_set), FUN = function(i) {
      ind = interaction_set[, i]
      return(((main_effects[[ind[1]]]) * (main_effects[[ind[2]]])) / temp_kernel$K1)
    })
  }
  
  if (length(clique_list) != 0) {
    K = temp_kernel$K1 + Reduce("+", main_effects[active_set]) + Reduce("+", interaction_kernel) + Reduce("+", high_order_kernel)
  } else {
    K = temp_kernel$K1 + Reduce("+", main_effects[active_set]) + Reduce("+", interaction_kernel)
  }
  return(K)
}
