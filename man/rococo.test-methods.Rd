\name{rococo.test-methods}
\docType{methods}
\alias{rococo.test}
\alias{rococo.test-methods}
\alias{rococo.test,numeric,numeric-method}
\alias{rococo.test,formula,data.frame-method}
\title{Robust Gamma Rank Correlation Test}
\description{
  Methods performing a robust gamma rank correlation test
}
\usage{
\S4method{rococo.test}{numeric,numeric}(x, y,
  similarity=c("linear", "exp", "gauss", "epstol", "classical"),
  tnorm="min", r=0, numtests=1000,
  alternative=c("two.sided", "less", "greater"))
\S4method{rococo.test}{formula,data.frame}(x, y, na.action, ...)
}
\arguments{
  \item{x}{a numeric vector or a formula; compulsory argument}
  \item{y}{compulsory argument; if \code{x} is a vector, \code{y}
    must be vector of the same length as \code{x}. If \code{x} is
    a formula, \code{y} must be a data frame.}
  \item{similarity}{a character string or a character vector identifying
    which type of similarity measure to use; see \code{\link{rococo}}
    for more details.}
  \item{tnorm}{t-norm used for aggregating results; see \code{\link{rococo}}
    for more details.}
  \item{r}{numeric vector defining the tolerances to be used;
    see \code{\link{rococo}} for more details.}
  \item{numtests}{number of random shuffles to perform; see details
    below.}
  \item{alternative}{indicates the alternative hypothesis and must be one of
    \code{"two.sided"}, \code{"greater"}, or
    \code{"less"}. Abbreviations are allowed as long as they are
    unique. \code{"greater"} corresponds to positive
    association, \code{"less"} to negative association.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}'s.  Defaults to \code{getOption("na.action")}.}
  \item{...}{all parameters specified are forwarded internally to the
    method \code{rococo.test} with signature \code{numeric,numeric}.}
}
\details{If called for numeric vectors, \code{rococo.test} computes the
  robust gamma rank correlation coefficient of \code{x} and \code{y}
  according to the specified parameters (see \code{\link{rococo}}) and
  then performs a permutation test to compute an estimated p-value.
  This is done by \code{numtests} random shuffles of \code{y} and
  computing the empirical standard deviation of the robust gamma
  correlation values for these shuffled data sets. Under the assumption
  that these values are normally distributed around mean zero, the
  p-value is then computed from this distribution in the usual way.

  Note that a too small choice of the number of shuffles (parameter
  \code{numtests}) leads to unreliable p-values.

  If \code{rococo.test} is called for a formula \code{x} and a data
  frame \code{y}, then the method checks whether the formula \code{x} correctly
  extracts two columns from \code{y} (see examples below). If so, the
  two columns are extracted and the robust gamma rank correlation test
  is applied to them according to the specified parameters.
}
\value{
  Upon successful completion, the function returns an object of class
  \code{\linkS4class{RococoTestResults}} containing the results.
}
\author{Martin Krone & Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}

U. Bodenhofer and F. Klawonn (2008).
Robust rank correlation coefficients on the basis of fuzzy orderings:
initial steps. \emph{Mathware Soft Comput.} \bold{15}(1):5-20.}
\seealso{\code{\link{rococo.test}}}
\examples{
## create data
f <- function(x)
{
  if (x > 0.9) x - 0.9
  else if (x < -0.9) x + 0.9
  else 0
}
x <- rnorm(25)
y <- sapply(x, f) + rnorm(25, sd=0.1)

## perform correlation tests
rococo.test(x, y, similarity="classical", alternative="greater")
rococo.test(x, y, similarity="linear", alternative="greater")
rococo.test(x, y, similarity=c("classical", "gauss"), r=c(0, 0.1),
            alternative="greater", numtests=10000)

## the formula variant
require(datasets)
data(iris)
rococo.test(~ Petal.Width + Petal.Length, iris, similarity="linear",
            alternative="two.sided")
}
\keyword{htest}
\keyword{methods}