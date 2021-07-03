\name{dbvsmsvm}
\alias{dbvsmsvm}
\title{
Derivative-based variable selection on the multicategory support vector machines
}
\description{
Derivative-based variable selection on the multicategory support vector machines.
}
\usage{
dbvsmsvm(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10, lambda_seq = c(2^{seq(-10, 15, length.out = 100)}, 1e+6),
		 thresh_Ngrid = 10, kernel = "linear", kparam = 1, scale = FALSE, criterion = "0-1", cv_type = "original", interaction = FALSE,
		 gd_scale = FALSE, optModel = FALSE, nCores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A \emph{n} x \emph{p} data matrix, where \emph{n} is the number of observations and \emph{p} is the number of variables.
  }
  \item{y}{
  A response vector with two and more label for each row of \emph{x}.
  }
  \item{gamma}{
  The convex combination parameter of the loss function.
  }
  \item{valid_x}{
  A validation data matrix for selecting \code{lambda} and threshold parameter \emph{v}. If \code{valid_x = NULL}, cross-validation is performed.
  }
  \item{valid_y}{
  A validation response vector.
  }
  \item{nfolds}{
  The number of folds for cross-validation.
  }
  \item{lambda_seq}{
  A sequence of regularization parameter to control a level of \emph{l_2}-penalty.
  }
  \item{thresh_Ngrid}{
  The number of grid for a sequence of threshold parameter.
  }
  \item{kernel}{
  A type of kernel.
  }
  \item{kparam}{
  A parameter needed for kernel.
  }
  \item{scale}{
  A logical indicating the variables to be scaled.
  }
  \item{criterion}{
  A type of criterion evaluating prediction performance.
  }
  \item{cv_type}{
  A type of cross-validation. If \code{cv_type = "osr"}, one standard error rule is adopted.
  }
  \item{interaction}{
  A logical. Whether to estimate interaction effect. 
  }
  \item{gd_scale}{
  A logical. Whether to scale the derivative-based measure.
  }
  \item{optModel}{
  A logical. Whether to obtain the optimal classification model.
  }
  \item{nCores}{
  The number of cores to use for parallel computing.
  }
  \item{...}{
  Other arguments that can be passed to ramsvm function.
  }
}

\value{
  An S3 object with the following slots
  \item{selected}{selected variables.}
  \item{gd}{Derivative-based measures}
  \item{lambda}{vector used for regularization path.}
  \item{v}{vector of network variability measured for each regularization level.}
  \item{opt_lambda}{The lambda that gives the optimal network.}
}

\references{Choi, H., J. Gim, S. Won, Y. J. Kim, S. Kwon, and C. Park, 2017: Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}, no. 1, 1-10. \cr
			Park, B., H. Choi, C. Park, 2021: Negative binomial graphical model with excess zeros.}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(ZILGM)
set.seed(1)
n = 100; p = 10; prob = 2 / p;
A = generate_network(p, prob, type = "random")
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 0.5, noise = 0.0)
					  
# Compute a sequence of regularization parameter
lambda_max = find_lammax(simul_dat$X)
lambda_min = 1e-4 * lambda_max
lambs = exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
nb2_fit = zilgm(X = simul_dat$X, lambda = lambs, family = "NBII", update_type = "IRLS", 
                do_boot = TRUE, boot_num = 30, sym = "OR", nCores = 10)

# To compute the regularization parameters automatically, use the argument nlambda
nb2_fit = zilgm(X = simul_dat$X, nlambda = 50, family = "NBII", update_type = "IRLS",
                do_boot = TRUE, boot_num = 30, sym = "OR", nCores = 10)
				
# Get estimated graph
est_graph = nb2_fit$network[[nb2_fit$opt_index]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.