\name{sramsvm}
\alias{sramsvm}
\title{
Structured angle-based multicategory support vector machine with reinforced multicategory loss.
}
\description{
Fitting a structured angle-based multicategory support vector machine with reinforced multicategory loss.
}
\usage{
sramsvm(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
	 lambda_seq = 2^{seq(-10, 10, length.out = 100)},
	 lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
	 kernel = c("linear", "gaussian", "gaussian2", "poly", "spline", "spline2", "spline-t"), kparam = 1,
	 scale = FALSE, criterion = c("0-1", "loss"),
	 isCombined = FALSE, cv_type = c("original", "osr"), type = c("type1", "type2"), 
	 optModel = FALSE, nCores = 1, verbose = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A \emph{n} x \emph{p} data matrix, where \emph{n} is the number of observations and \emph{p} is the number of variables.
  }
  \item{y}{
  A response variable with three and more labels.
  }
  \item{gamma}{
  The convex combination parameter of the loss function. 
  }
  \item{valid_x}{
  A validation data matrix for selecting \code{lambda} and threshold parameter \emph{v} (optional). If \code{valid_x=NULL}, \code{nfolds}-fold cross-validation is performed.
  }
  \item{valid_y}{
  A validation response vector (optional).
  }
  \item{nfolds}{
  The number of folds for cross-validation.
  }
  \item{lambda_seq}{
  A sequence of regularization parameter to control a level of \emph{l_2}-penalty.
  }
  \item{lambda_theta_seq}{
  A sequence of regularization parameter on theta's, which are weights for kernels.
  }
  \item{kernel}{
  A character string representing one of type of kernel.
  }
  \item{kparam}{
  A parameter needed for kernel.
  }
  \item{scale}{
  A logical indicating the variables to be scaled.
  }
  \item{criterion}{
  A type of criterion evaluating prediction performance of cross-validation.
  }
  \item{isCombined}{
  A logical. If \code{isCombined = TRUE}, the theta-step is combined with the following c-step. 
  }
  \item{cv_type}{
  A rule of selecting lambda_theta parameter. If \code{cv_type="osr"}, one standard error rule is adopted.
  }
  \item{type}{
  A type of optimization method for ramsvm. If \code{type = "type1"}, intercepts are included in a regularization term. If \code{type = "type2"}, intercepts are omitted in the regularization term and estimated by a linear programming 
  }
  \item{optModel}{
  A logical. Whether to obtain the optimal classification model.
  }
  \item{nCores}{
  The number of cores to use for parallel computing.
  }
  \item{verbose}{
  A integer value. If \code{verbose = 1}, the minimum values of cross-validation error for each step are printed. 
  }
  \item{...}{
  Other arguments that can be passed to ramsvm core function.
  }
}

\value{
  An S3 object of class "\code{sramsvm}" containing the following slots
  \item{opt_theta}{The theta vector at the optimal lambda_theta value.}
  \item{cstep_inform}{Information of c-step including cross-validation errors and optimal lambda value.}
  \item{theta_step_inform}{Information of theta-step including cross-validation errors and optimal lambda_theta value.}
  \item{opt_model}{If \code{optModel=TRUE}, classification model with the optimal lambda and lambda_theta is returned.}
  \item{opt_valid_err}{If \code{optModel=TRUE}, a minimum value of cross-validation errors based on the optimal lambda_theta is returned.}
  \item{valid_err}{If \code{optModel=TRUE}, cross-validation errors based on the optimal lambda_theta is returned.}
  \item{cv_type}{A rule of selecting lambda_theta.}
  \item{call}{The call of \code{sramsvm}.}
}

\references{}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(dbvsmsvm)
n = 100; p = 10;
data = dbvsmsvm:::sim_gen(n = n, p = p, type = "linear")
x = scale(data$x)
y = data$y
sigma = kernlab::sigest(y ~ x, scaled = FALSE)[3]

# Fitting sramsvm with the linear kernel
sram_linear = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                      lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                      lambda_theta_seq = 2^{seq(-20, 0, length.out = 100)},
                      kernel = "linear", scale = FALSE, criterion = "0-1",
                      isCombined = TRUE, cv_type = "osr", nCores = 1)

# Fitting the sramsvm with the Gaussian kernel
sram_radial = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                      lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                      lambda_theta_seq = 2^{seq(-20, 0, length.out = 100)},
                      kernel = "gaussian", kparam = sigma, scale = FALSE, criterion = "0-1",
                      isCombined = TRUE, cv_type = "osr", nCores = 1)

# Fitting the sramsvm with the Gaussian kernel with second-order interaction
sram_radial_interaction = sramsvm(x = x, y = y, gamma = 0.5, nfolds = 5,
                      lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                      lambda_theta_seq = 2^{seq(-20, 0, length.out = 100)},
                      kernel = "gaussian", kparam = sigma, scale = FALSE, criterion = "0-1",
                      isCombined = TRUE, cv_type = "osr", nCores = 1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.