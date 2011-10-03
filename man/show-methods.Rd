\name{show-methods}
\docType{methods}
\alias{show}
\alias{show-methods}
\alias{show,RococoTestResults-method}
\alias{show.RococoTestResults}
\title{Display RococoTestResults Objects}
\description{
  Display method for S4 class \code{\linkS4class{RococoTestResults}}
}
\usage{
\S4method{show}{RococoTestResults}(object)
}
\arguments{
  \item{object}{an object of class
    \code{\linkS4class{RococoTestResults}}
  }
}
\details{
  \code{show} displays the most important information stored in
  \code{object}.
}
\value{
  \code{show} returns an invisible \code{NULL}
}
\author{Martin Krone & Ulrich Bodenhofer
\email{rococo@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/rococo/}

U. Bodenhofer and F. Klawonn (2008).
Robust rank correlation coefficients on the basis of fuzzy orderings:
initial steps. \emph{Mathware Soft Comput.} \bold{15}(1):5-20.}
\seealso{\code{\link{rococo.test-methods}}, \code{\link{rococo}},
\code{\link{show-methods}}}
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
ret <- rococo.test(x, y, similarity="classical", alternative="greater")
show(ret)

ret <- rococo.test(x, y, similarity="linear", alternative="greater")
show(ret)

ret <- rococo.test(x, y, similarity=c("classical", "gauss"), r=c(0, 0.1),
                   alternative="greater", numtests=10000)
show(ret)
}
\keyword{htest}
\keyword{methods}
