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
  \item{selected}{Selected variables.}
  \item{gd}{Derivative-based measures.}
  \item{threshold_path}{A sequence of threshold parameter.}
  \item{opt_threshold}{Optimal threshold value.}
  \item{valid_err}{Validation errors.}
  \item{opt_valid_err}{A minimum value of validation errors.}
  \item{int_selected}{Selected variables of interaction selection.}
  \item{gd_interaction}{Derivative-based measures of interaction selection.}
  \item{threshold_path_int}{A sequence of threshold parameter of interaction selection.}
  \item{opt_threshold_int}{Optimal threshold value of interaction selection.}
  \item{int_valid_err}{Validation errors for interaction selection.}
  \item{int_opt_valid_err}{A minimum value of validation errors of interaction selection.}
  \item{cv_type}{A type of cross-validation.}
  \item{call}{The call of dbvsmsvm}
}

\references{}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(dbvsmsvm)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.