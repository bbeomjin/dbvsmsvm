\name{dbvsmsvm}
\alias{dbvsmsvm}
\title{
Derivative-based variable selection for the angle-based multicategory support vector machines
}
\description{
Derivative-based variable selection for the angle-based multicategory support vector machines with reinforced multicategory loss.
}
\usage{
dbvsmsvm(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10, 
	 lambda_seq = 2^{seq(-10, 15, length.out = 100)},
	 v_seq = NULL, Nofv = 100, u_seq = NULL, Nofu = 100,
	 kernel = c("linear", "gaussian"), kparam = 1, scale = FALSE, 
	 criterion = c("0-1", "hinge"), cv_type = c("standard", "osr"), 
	 interaction = FALSE, optModel = FALSE, nCores = 1, ...)
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
  \item{v_seq}{
  A sequence of threshold parameter to select informative variables.
  }  
  \item{Nofv}{
  The number of threshold values. If \code{v_seq} is not given, a sequence of threshold parameter is automatically computed based on \code{Nofv}. 
  }
  \item{u_seq}{
  A sequence of threshold parameter for interaction selection.
  }
  \item{Nofu}{
  The number of threshold values for interaction selection. If \code{u_seq} is not given, a sequence of threshold parameter is automatically computed based of \code{Nofu}.
  }
  \item{kernel}{
  A character string representing one of type of kernel.
  }
  \item{kparam}{
  A parameter needed for kernel.
  }
  \item{scale}{
  A logical value indicating whether to scale the variables. If \code{scale=TRUE}, \code{x} is scaled to zero mean and unit variance.
  }
  \item{criterion}{
  A type of criterion evaluating prediction performance of cross-validation.
  }
  \item{cv_type}{
  A rule of selecting threshold parameter. If \code{cv_type="osr"}, one standard error rule is adopted.
  }
  \item{interaction}{
  A logical. Whether to select interaction terms. 
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
  An S3 object of class "\code{dbvsmsvm}" containing the following slots
  \item{selected}{A zero-one vector indicating the selected variables.}
  \item{lambda_cv_inform}{Information of cross-validation for lambda including cross-validation errors and optimal lambda value.}
  \item{selection_cv_inform}{Information of cross-validation for threshold parameter including cross-validation errors and optimal threshold value.}
  \item{interaction_selection_cv_inform}{Information of cross-validation for threshold parameter of interaction selection including cross-validation errors and optimal threshold value.}
  \item{opt_model}{If \code{optModel=TRUE}, classification model with the optimal lambda and threshold parameter is returned.}
  \item{opt_valid_err}{If \code{optModel=TRUE}, a minimum value of cross-validation errors based on the selected variables is returned.}
  \item{valid_err}{If \code{optModel=TRUE}, cross-validation errors based on the selected variables is returned.}
  \item{cv_type}{A rule of selecting threshold parameter.}
  \item{call}{The call of \code{dbvsmsvm}.}
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

# Fit the DBVS-MSVM with the linear kernel
dbvs_linear = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                       Nofv = 100, kernel = "linear", criterion = "0-1", scale = FALSE,
                       cv_type = "osr", interaction = FALSE, gamma = 0.5, optModel = FALSE, nCores = 1)

# Fit the DBVS-MSVM with the Gaussian kernel
dbvs_radial = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                       Nofv = 100, kernel = "gaussian", kparam = sigma, criterion = "0-1", scale = FALSE,
                       cv_type = "osr", interaction = FALSE, gamma = 0.5, optModel = TRUE, nCores = 1)

# Fit the DBVS-MSVM with the Gaussian kernel for selecting second-order interaction
dbvs_interaction = dbvsmsvm(x = x, y = y, nfolds = 5, lambda_seq = 2^{seq(-20, 0, length.out = 100)},
                            Nofv = 100, kernel = "gaussian", kparam = sigma, criterion = "0-1", scale = FALSE,
                            cv_type = "osr", interaction = TRUE, gamma = 0.5, optModel = FALSE, nCores = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.