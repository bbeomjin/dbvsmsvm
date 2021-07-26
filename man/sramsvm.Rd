\name{sramsvm}
\alias{sramsvm}
\title{
Structured angle-based multicategory support vector machine with reinforced multicategory loss.
}
\description{
Fit a structured angle-based multicategory support vector machine with reinforced multicategory loss.
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
  A response vector with three and more labels.
  }
  \item{gamma}{
  The convex combination parameter of the loss function.
  }
  \item{lambda}{
  A regularization parameter to control a level of \emph{l_2}-penalty.
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
  \item{type}{
  A type of optimization method for ramsvm. If \code{type = "type1"}, intercepts are included in a regularization term. If \code{type = "type2"}, intercepts are omitted in the regularization term and estimated by a linear programming 
  }
  \item{...}{
  Other arguments that can be passed to ramsvm core function.
  }
}

\value{
  An S3 object with "\code{ramsvm}" containing the following slots
  \item{x}{Input \code{x}.}
  \item{K} {Input \code{K}.}
  \item{y}{Input \code{y}.}
  \item{y_name}{The class labels of \code{y}.}
  \item{gamma}{The convex combination parameter of the loss function.}
  \item{n_class}{The number of classes.}
  \item{lambda}{Given regularization parameter.}
  \item{kernel}{Given type of kernel.}
  \item{kparam}{Given parameter for kernel.}
  \item{cmat}{The corresponding coefficients.}
  \item{c0vec}{The intercepts.}
  \item{alpha}{The Lagrange multipliers.}
  \item{fit_class}{Fitted class.}
  \item{epsilon}{Convergence tolerance in the ramsvm core algorithm.}
}

\references{}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(dbvsmsvm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.