\name{ramsvm}
\alias{ramsvm}
\title{
Multicategory support vector machine with reinforced multicategory loss.
}
\description{
A function for multicategory support vector machine with reinforced multicategory loss.
}
\usage{
ramsvm(x, y, gamma = 0.5, lambda, kernel, kparam, scale = FALSE, type = c("type1", "type2"), ...)
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
  \item{lambda}{
  A regularization parameter to control a level of \emph{l_2}-penalty.
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
  \item{type}{
  A type of ramsvm. If \code{type = "type1"}, the intercepts included in regularization.
  }
  \item{...}{
  Other arguments that can be passed to ramsvm core function.
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