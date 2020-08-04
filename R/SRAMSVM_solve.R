# dyn.load("./alpha_update2.dll")
SRAMSVM_solve = function(x = NULL, K = NULL, y, gamma, lambda, kernel, kparam,
                         weight = NULL, epsilon = 1e-4 * length(y) * length(unique(y)), warm = NULL) {

  if (is.null(K)) {
    K = kernelMat(x, x, kernel = kernel, kparam = kparam) + 1
    # K = kernelMatrix(rbfdot(sigma = kparam), x, x) + 1
  } else {
    K = K + 1
  }

  #------------------------------------------------------------------#
  # Convert labels to integers.                                      #
  #------------------------------------------------------------------#
  y_temp = as.factor(y)
  y_name = levels(y_temp)
  k = length(y_name)

  y_train = integer(length(y))
  for(ii in 1:k) y_train[which(y_temp %in% y_name[ii])] = ii
  if(is(y, "numeric")) y_name = as.numeric(y_name)

  #------------------------------------------------------------------#
  # Calculate kernel                                                 #
  #------------------------------------------------------------------#

  kdouble = as.double(k)
  nobs = length(y)
  nobsdouble = as.double(nobs)
  warm = matrix(data = 0.0, nrow = nobs, ncol = k)
  if (is.null(weight)) weight = numeric(nobs) + 1

  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  my = t(XI_gen(k = k, kd = kdouble))

  yyi = Y_matrix_gen(k = k,
                     kd = kdouble,
                     nobs = nobs,
                     y_train = y_train)

  fold = max(floor(min(summary(as.factor(y_train))) / 30), 1)

  # labidx = numeric(0)
  # len = numeric(k)
  # sam = numeric(0)
  # for(i in 1:k) {
  #   labidx = c(labidx, list(which(y_train == i)))
  #   len[k] = round( length(labidx[[i]]) / fold )
  #   sam = c(sam, list(sample(1L:length(labidx[[i]])) ))
  # }

  betaout = list()
  beta0out = list()

  for(count in 1:length(lambda)) {
    templambda = lambda[count]

    alpha_ij = warm
    alpha_yi = numeric(nobs)

    for(zz in 1:nobs) {alpha_yi[zz] = alpha_ij[zz, y_train[zz]]}

    erci = -diag(K) / 2 / nobsdouble / templambda

    aa = .C("alpha_update",
            as.vector(alpha_ij),
            as.vector(alpha_yi),
            as.vector(my),
            as.vector(yyi),
            as.vector(K),
            as.double(templambda),
            as.vector(weight),
            as.integer(nobs),
            as.double(nobsdouble),
            as.integer(k),
            as.double(kdouble),
            as.vector(erci),
            as.double(gamma),
            as.vector(y_train),
            as.double(epsilon),
            outalpha_ij = as.vector(numeric(nobs * k)))

    warm = matrix(data = aa$outalpha_ij, nrow = nobs, ncol = k)


    beta = beta_kernel(y = y_train,
                       k = k,
                       my = my,
                       warm = warm,
                       lambda = templambda)
    # drop(crossprod(beta[[1]][, 1], X))

    # tt = beta_linear(x = x, y = y_train, k = k, my = my, warm = warm, lambda = templambda)

    # beta0 = matrix(find_theta2(y, K, gamma = gamma, cmat = beta$beta, lambda = lambda), nrow = 1)


    betaout[[count]] = beta$beta
    beta0out[[count]] = beta$beta0
    # beta0out[[count]] = beta0

  }

  out = list()
  out$x = x
  out$y = y
  out$y_name = y_name
  out$k = k
  out$gamma = gamma
  out$weight = weight
  out$lambda = lambda
  out$kparam = kparam
  out$beta = betaout
  out$beta0 = beta0out
  out$epsilon = epsilon
  out$warm = warm
  out$kernel = kernel
  class(out) = "sramsvm"
  return(out)
}


predict.sramsvm = function(object, newx = NULL, newK = NULL, lambda = NULL, ...) {

  if (is.null(x = lambda)) {lambda = object$lambda}
  # if (is.null(newx)) {newx = object$x}
  if (is.null(newK)) {
    newK = kernelMat(newx, object$x, kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  if( !is.numeric(x = lambda) ) {
    stop("All lambda should be numeric.", call. = FALSE)
  }

  if( any(lambda < 0) ) {
    stop("All lambda should be non-negative.", call. = FALSE)
  }

  pred_y = numeric(0)

  nol = length(object$lambda)

  for( i in 1L:length(lambda) ) {
    temp = lambda[i]
    index = which(sapply(X = object$lambda,
                         FUN = function(x){isTRUE(all.equal(x, temp))}))
    if(length(index) == 1) {
      temp_beta = object$beta[[index]]
      temp_beta0 = object$beta0[[index]]
    } else if(length(index) == 0) {
      if( temp > object$lambda[1] ) {
        temp_beta = object$beta[[1]]
        temp_beta0 = object$beta0[[1]]
        cat(paste("Lambda",
                  temp,
                  "is bigger than the largest lambda in the",
                  "solution path.\nUsing the parameters",
                  "corresponding to the largest lambda in the",
                  "solution path.\n"))
      } else if( temp < object$lambda[nol] ) {
        temp_beta = object$beta[[nol]]
        temp_beta0 = object$beta0[[nol]]
        cat(paste("Lambda",
                  temp,
                  "is less than the smallest lambda in the",
                  "solution path.\nUsing the parameters",
                  "corresponding to the smallest lambda in the",
                  "solution path.\n"))
      } else if( temp > object$lambda[nol] &&
                 temp < object$lambda[1] ) {
        idx = max(which(object$lambda > temp))
        temp_beta = object$beta[[idx]] +
          (temp - object$lambda[idx]) /
          (object$lambda[(idx + 1)] - object$lambda[idx]) *
          (object$beta[[(idx + 1)]] - object$beta[[idx]])
        temp_beta0 = object$beta0[[idx]] +
          (temp - object$lambda[idx]) /
          (object$lambda[(idx + 1)] - object$lambda[idx]) *
          (object$beta0[[(idx + 1)]] - object$beta0[[idx]])
      }
    }

    temp_pred_y = predict_kernel(K_test = newK,
                                 beta = temp_beta,
                                 beta0 = temp_beta0,
                                 k = object$k)
    inner_prod = temp_pred_y$inner_prod
    temp_pred_y = temp_pred_y$class

    temp_pred_y2 = character(nrow(newK))
    for(i in 1:object$k) {
      temp_pred_y2[temp_pred_y == i] = object$y_name[i]
    }

    if( is.numeric(x = object$y) ) {
      temp_pred_y2 = as.numeric(temp_pred_y2)
    }

    pred_y = c(pred_y, list(temp_pred_y2))
  }

  names(pred_y) = lambda

  return(list(class = pred_y, inner_prod = inner_prod))
}


predict_kernel = function(K_test, beta, beta0, k)
{
  n = nrow(K_test)

  XI = XI_gen(k = k, kd = as.double(k))

  beta0 = matrix(beta0,
                 nrow = n,
                 ncol = ncol(beta),
                 byrow = TRUE)

  f_matrix = t(K_test %*% beta + beta0)

  inner_matrix = matrix(data = 0, nrow = n, ncol = k)

  for(ii in 1:k) inner_matrix[, ii] = colSums(f_matrix * XI[, ii])

  z = apply(X = inner_matrix, MARGIN = 1, FUN = pred)

  return(list(class = z, inner_prod = inner_matrix))

}

pred = function(f) {
  tst = sapply(f, function(i) {isTRUE(all.equal(i, max(f)))})
  y = min(which(tst))
  return(y)
}

beta_kernel = function(y, k, my, warm, lambda){

  nobs = length(y)
  dnobs = as.double(nobs)

  beta = matrix(data = 0, nrow = nobs, ncol = (k - 1))
  beta0 = matrix(data = 0, nrow = 1, ncol = (k - 1))

  for(q in 1:(k - 1)) {
    temp = numeric(nobs)
    temp0 = 0
    for( ii in 1L:nobs ) {
      for( jj in 1:k ) {
        t1 = warm[ii, jj] * my[jj, q]
        t2 = (2 * {y[ii] == jj} - 1)

        temp[ii] = temp[ii] + t2 * t1
        temp0 = temp0 + t2 * t1
      }
    }

    beta[, q] = temp / dnobs / lambda
    beta0[, q] = temp0 / dnobs / lambda
  }
  return(list(beta = beta, beta0 = beta0))
}

