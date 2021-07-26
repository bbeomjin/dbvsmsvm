\name{cv.ramsvm}
\alias{cv.ramsvm}
\title{
Cross-Validation for \code{ramsvm}.
}
\description{
Perform cross-validation for the optimal lambda of \code{ramsvm}.
}
\usage{
cv.ramsvm(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10,
	 lambda_seq = 2^{seq(-10, 15, length.out = 100)},
	 kernel = c("linear", "gaussian"), kparam = 1, 
	 scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A \emph{n} x \emph{p} data matrix, where \emph{n} is the number of observations and \emph{p} is the number of variables.
  }
  \item{y}{
  A response vector with three and more labels.
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
  An S3 object with "\code{ramsvm}" containing the following slots
  \item{opt_model}{If \code{optModel=TRUE}, classification model with the optimal lambda is returned.}
  \item{opt_valid_err}{If \code{optModel=TRUE}, a minimum value of cross-validation errors based on the selected variables is returned.}
  \item{valid_err}{If \code{optModel=TRUE}, cross-validation errors based on the selected variables is returned.}
  \item{call}{The call of \code{cv.ramsvm}.}
}

\references{}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(dbvsmsvm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.