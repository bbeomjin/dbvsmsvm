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
  \item{x}{Input \emph{x}.}
  \item{y}{Input \emph{y}.}
  \item{y_name}{The class labels of \emph{y}.}
  \item{gamma}{The convex combination parameter of the loss function.}
  \item{n_class}{The number of classes.}
  \item{lambda}{The regularization parameter.}
  \item{kernel}{Given type of kernel.}
  \item{kparam}{Given parameter for kernel.}
  \item{cmat}{The estimated parameters of the classification functions.}
  \item{c0vec}{The estimated intercepts of the classification functions.}
  \item{alpha}{A solutions of quadratic programming.}
  \item{scale}{}
  \item{center}{}
  \item{scaled}{}
  \item{fit_class}{}
  \item{epsilon}{}
}

\references{}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(dbvsmsvm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.