kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0) {
  
  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)
  
  if (NCOL(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }
  
  if (NCOL(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }
  
  if (kernel == "poly") {
    K = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "gaussian" | kernel == "gaussian2") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # obj = kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "spline") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  } else if (kernel == "linear") {
    K = tcrossprod(x, y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else {
    K = NULL
  }
  return(K)
}

spline_kernel = function(x, y)
{
  x = as.matrix(x)
  y = as.matrix(y)
  K1x = (x - 1 / 2)
  K1y = (y - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2y = (K1y^2 - 1 / 12) / 2
  ax = x %x% matrix(1, 1, nrow(y)) 
  ay = y %x% matrix(1, 1, nrow(x))
  b = abs(ax - t(ay))
  K1 = K1x %x% t(K1y)
  K2 = K2x %x% t(K2y) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  list(K1 = K1, K2 = K2)
}


make_anovaKernel = function(x, y, kernel, kparam)
{
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)
  
  # calculate anova kernels for main effects
  if (kernel == "spline") {
    # assign the number of anova kernels
    numK = 2 * dimx
    # list of kernel matrices
    anova_kernel = vector(mode = "list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    } 
    
  } else if (kernel == 'spline2') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }  
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])            
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (kernel == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (kernel == 'spline-t2') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else if (kernel == "gaussian2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[index]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
    
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[d]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  return(list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel, kparam = kparam))
}

combine_kernel = function(anova_kernel, theta = rep(1, anova_kernel$numK))
{
  K = 0
  for (d in 1:anova_kernel$numK) {
    K = (K + theta[d] * anova_kernel$K[[d]])
  }
  return(K)
}



find_intercept = function(y, K, gamma = 0.5, cmat, lambda) 
{
  # standard LP form :
  # min a^T x , subject to A1x <= a1
  n_data = length(y)
  n_class = length(unique(y))
  # bvec = as.matrix(bvec)
  
  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, nobs = n_data, y = y)
  
  # calculate the 'a' matrix 
  # a_tmp = matrix((1 - gamma) / n_data, nrow = n_data, ncol = n_class)
  # a_tmp[cbind(1:n_data, y)] = gamma / n_data
  a_tmp = matrix((1 - gamma), nrow = n_data, ncol = n_class)
  a_tmp[cbind(1:n_data, y)] = gamma
  a = matrix(a_tmp, ncol = 1)
  a = rbind(a, matrix(rep(0, n_class - 1)))
  
  Y_code = Y_matrix_gen(n_class, nobs = n_class, y = 1:n_class)
  
  sign_mat = matrix(1, n_data, n_class)
  sign_mat[cbind(1:n_data, y)] = -1
  sign_vec = as.vector(sign_mat)
  # constraints
  
  for (j in 1:n_class) {
    if (j == 1) {
      temp_N = matrix(rowSums(K %*% cmat * matrix(Y_code[j, ], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE)))
    } else {
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



XI_gen = function(k) {
  
  tempA = - (1.0 + sqrt(k)) / ((k - 1.0)^(1.5))
  tempB = tempA + sqrt(k / (k - 1.0))
  
  XI = matrix(data = tempA, nrow = k - 1L, ncol = k)
  
  XI[, 1L] = 1.0 / sqrt(k - 1.0)
  
  for (ii in 2:k) XI[ii - 1, ii] = tempB
  
  return(XI)
}

Y_matrix_gen = function(k, nobs, y) {
  
  XI = XI_gen(k = k)
  
  Y_matrix = matrix(data = 0.0, nrow = nobs, ncol = k - 1)
  
  for (ii in 1:nobs) {Y_matrix[ii, ] = XI[, y[ii]]}
  
  return(Y_matrix)
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


dgaussian = function(X, xp, kparam = 1)
{
  gamma = kparam
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

ddgaussian = function(X, xp, kparam = 1, comb_set)
{
  gamma = kparam
  np = dim(X)
  diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
  # diff_mat2 = diff_mat[, comb_set]
  
  diff_mat2 = sapply(1:ncol(comb_set), function(i) diff_mat[, comb_set[, i][1]] * diff_mat[, comb_set[, i][2]])
  
  # tmp = exp(-rowSums(diff_mat^2) * gamma) * 4 * gamma^2 * (diff_mat[, ind1]) * (diff_mat[, ind2])
  tmp = exp(-rowSums(diff_mat^2) * gamma) * 4 * gamma^2 * diff_mat2
  # dgd = drop(crossprod(alpha, tmp))
  return(tmp)
}



dlinear = function(X, xp, kparam = 1)
{
  return(X)
}


dpoly = function(X, xp, kparam = 1)
{
  # degree = kernel_par$degree
  # scale = kernel_par$scale
  # offset = kernel_par$offset
  degree = kparam
  scale = 1
  offset = 0
  tmp = degree * (scale * drop(X %*% xp) + offset)^{degree - 1} * scale * X
  return(tmp)
}

ddspline = function(x, y)
{
  m1 = (x - (1 / 2)) + (1 / 2) * ((x - (1 / 2))^2 - (1 / 12)) * (y - (1 / 2))
  m2 = ifelse(x > y, -(1 / 6) * (x - y - (1 / 2))^3 + (1 / 24) * (x - y - (1 / 2)), 
              (4 * (y - x - (1 / 2))^3 - (y - x - (1 / 2))) / 24)
  return(m1 - m2)
}


pderiv = function(alpha, x, y, kernel = c("linear", "poly", "gaussian", "spline", "anova_gaussian"), kparam = 1)
{
  n = length(y)
  k = length(unique(y))
  
  dkernel = switch(kernel,
                   linear = dlinear,
                   poly = dpoly,
                   gaussian = dgaussian,
                   spline = dspline,
				           anova_gaussian = dgaussian)
  
  grad_mat = 0
  for (i in 1:n) {
    dK_sq = crossprod(alpha, dkernel(x, x[i, ], kparam))^2
    grad_mat = grad_mat + dK_sq / n
  }
  
  res = colSums(grad_mat)
  return(res)
}



pderiv_so = function(alpha, x, y, kernel = c("linear", "poly", "gaussian"), kparam = 1, active_set = NULL)
{
  n = length(y)
  k = length(unique(y))
  
  ddkernel = switch(kernel,
                    linear = ddlinear,
                    poly = ddpoly,
                    gaussian = ddgaussian)
  
  comb_set = combn(active_set, 2)
  
  grad_mat = 0
  for (i in 1:n) {
    dK_sq = crossprod(alpha, ddkernel(x, x[i, ], kparam, comb_set))^2
    grad_mat = grad_mat + dK_sq / n
  }
  
  res = colSums(grad_mat)
  return(res)
}

data_split = function(y, nfolds, seed = length(y))
{
  # k: the number of classes
  y = as.factor(y)
  n_data = length(y)
  n_class = length(levels(y))
  class_size = table(y)
  classname = names(class_size)
  ran = rep(0, n_data) 
  if ((min(class_size) < nfolds) & (nfolds != n_data))
  {
    warning('The given fold is bigger than the smallest class size. \n Only a fold size smaller than the minimum class size \n or the same as the sample size (LOOCV) is supported.\n')
    return(NULL)
  }
  
  if (min(class_size) >= nfolds) {
    set.seed(seed)
    for (j in 1:n_class) {  
      ran[y == classname[j]] = ceiling(sample(class_size[j]) / (class_size[j] + 1) * nfolds) 
    }
  }
  else if (nfolds == n_data) {
    ran = 1:n_data
  }
  return(ran)
}



strong_heredity = function(main_effect, interaction)
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


# interaction_graph = function(comb, p, min = 3)
# {
#   int_mat = Matrix::Matrix(0, nrow = p, ncol = p)
#   int_mat[t(comb)] = 1
#   int_mat = Matrix::t(int_mat) + int_mat
#   g = graph_from_adjacency_matrix(int_mat, mode = "undirected")
#   cliques_list = max_cliques(g, min = 3)
#   return(cliques_list)
# }

interaction_kernel = function(x, u, kernel, kparam, active_set, interaction_set)
{
  if (!is.matrix(x)) {
    x = as.matrix(x)
  } else {
    x = as.matrix(x)
  }
  u = as.matrix(u)
  dimx = ncol(x)
  
  scaled_kernel = function(x, u, kernel, kparam, active_set, index)
  {
    X1 = matrix(rowSums(x[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U1 = matrix(rowSums(u[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    X2 = matrix(rowSums(x[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U2 = matrix(rowSums(u[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    K = exp(-kparam * ((X1 + U1) - (X2 + U2)))
    K1 = exp(-kparam * (X1 + U1))
    K2 = exp(-kparam * (X2 + U2))
    K_mat = kernelMatrix(x[, index, drop = FALSE], u[, index, drop = FALSE], kernel = kernel, kparam = kparam)
    res = K * K_mat - K1
    return(list(res = res, K = K, K1 = K1, K2 = K2))
  }
  
  main_effects = vector(mode = "list", dimx)
  
  for (j in 1:length(active_set)) {
    temp_kernel = scaled_kernel(x, u, kernel = kernel, kparam = kparam, active_set = active_set, index = active_set[j])
    main_effects[[active_set[j]]] = temp_kernel$res
  }
  
  interaction_kernel = 0
  if (ncol(interaction_set) != 0) {
    
    for (i in 1:ncol(interaction_set)) {
      ind = interaction_set[, i]
      interaction_kernel = interaction_kernel + ((main_effects[[ind[1]]]) * (main_effects[[ind[2]]])) / temp_kernel$K1
    }
  }
  
  K = temp_kernel$K1 + Reduce("+", main_effects[active_set]) + interaction_kernel
  
  return(K)
}

code_ramsvm = function(y)
{
  n_class = length(unique(y))
  n = length(y)
  yyi = Y_matrix_gen(k = n_class, nobs = n, y = y)
  W = XI_gen(n_class)
  
  y_index = cbind(1:n, y)
  index_mat = matrix(-1, nrow = n, ncol = n_class)
  index_mat[y_index] = 1
  
  Hmatj = list()
  Lmatj = list()
  for (j in 1:(n_class - 1)) {
    Hmatj_temp = NULL
    Lmatj_temp = NULL
    for (i in 1:n_class) {
      temp = diag(n) * W[j, i]
      diag(temp) = diag(temp) * index_mat[, i]
      Hmatj_temp = rbind(Hmatj_temp, temp)
      Lmatj_temp = c(Lmatj_temp, diag(temp))
    }
    Hmatj[[j]] = Hmatj_temp
    Lmatj[[j]] = Lmatj_temp
  }
  return(list(yyi = yyi, W = W, Hmatj = Hmatj, Lmatj = Lmatj, y_index = y_index))
}

find_nonzero = function(Amat)
{
  nr = nrow(Amat)
  nc = ncol(Amat)
  Amat_compact = matrix(0, nr, nc)
  Aind = matrix(0, nr + 1, nc)
  for (j in 1:nc) {
    index = (1:nr)[Amat[, j] != 0]
    number = length(index)
    Amat_compact[1:number, j] = Amat[index, j]
    Aind[1, j] = number
    Aind[2:(number+1), j] = index
  }
  max_number = max(Aind[1, ])
  Amat_compact = Amat_compact[1:max_number, ]
  Aind = Aind[1:(max_number + 1), ]
  return(list(Amat_compact = Amat_compact, Aind = Aind))
}

fixit = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
{
  if (is_diag) {
    d = diag(A)
    tol = epsilon
    eps = max(tol * max(d), 0)
    d[d < eps] = eps
    Q = diag(d)
  } else {
    eig = eigen(A, symmetric = TRUE)
    tol = epsilon
    eps = max(tol * abs(eig$values[1]), 0)
    eig$values[eig$values < eps] = eps
    Q = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  return(Q)
}