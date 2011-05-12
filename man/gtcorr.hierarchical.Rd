\name{gtcorr.hierarchical}
\alias{gtcorr.hierarchical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the efficiency of hierarchical group testing procedures}
\description{
\sQuote{gtcorr.hierarchical} calculates the efficiencies of hierarchical group
testing procedures, allowing for correlation between units and test error
}
\usage{
gtcorr.hierarchical(n, m = 1, p, sigma = 0, se = 1, sp = 1, arrangement
= c("nested", "random"), model = c("beta-binomial", "Madsen",
"Morel-Neerchal"), ...)
}

\arguments{
  \item{n}{
    a numeric vector of pool sizes where \code{n[s]} is the size of a pool in stage
    \code{s}.  The size of a pool in the last stage is 1, which can be omitted.
  }
  \item{m}{
    cluster size. 
  }
  \item{p}{
    probability of a unit testing positive.
  }
  \item{sigma}{
    pairwise correlation of two units in a cluster.
  }
  \item{se}{
    sensitivity.  The probability that a pool of units tests positive given
    than at least one unit in that pool is positive
  }
  \item{sp}{
    specificity.  The probability that a pool of units tests negative
    given that at least one unit in that pool is negative 
  }
  \item{arrangement}{
    how clusters are arranged. Should be \sQuote{nested} or \sQuote{random}.
  }
  \item{model}{
    probability model for clusters. Should be \sQuote{beta-binomial},
  \sQuote{Madsen}, or \sQuote{Morel-Neerchal}.
  }
  \item{\dots}{
    \describe{
      \item{\code{runs}}{for a random arrangement, number of Monte Carlo simulations to
	perform to calculate the probability of a pool having no positive
	units. Default is 1000.
      }
    }
  }
}
\details{
  One of \code{m}, \code{p}, \code{sigma}, \code{se}, or \code{sp} can have
  more than one value.
  \code{m} should not be greater than \code{n[1]}.  For a
  \sQuote{nested} arrangement, \code{m} should be divisible by \code{n[s]} or
  \code{n[s]} should be divisible by \code{m} for all \code{s}.
  See Lendle et. al. 2011 for more information.
}

\value{
  \item{n}{number of units per pool at each stage.}
  \item{param.grid}{a data frame containing the values of \code{p},
    \code{sigma}, \code{se}, \code{sp}, and \code{m} for each value of
    \code{efficiency}.}
  \item{arrangement}{arrangement.}
  \item{model}{model.}
  \item{efficiency}{a vector of efficiencies, one for each row of
    \code{param.grid}.}
}
\references{
Samuel D. Lendle, Michael Hudgens, and Bahjat F. Qaqish, "Group Testing
for Case Identification with Correlated Responses" Submitted 2011.
Biometrics.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##Plot efficiencies of a Dorfman (2 stage hierarchical) algorithm 
##by cluster size and sigma
m <- 2^(0:8)
sig.0 <- gtcorr.hierarchical(n=256, p=.001, m=m, sigma=0)$efficiency
sig.05 <- gtcorr.hierarchical(n=256, p=.001, m=m, sigma=0.05)$efficiency
sig.5 <- gtcorr.hierarchical(n=256, p=.001, m=m, sigma=0.5)$efficiency
sig.99 <- gtcorr.hierarchical(n=256, p=.001, m=m, sigma=.99)$efficiency

plot(m, sig.99, log="x", type='b', ylab="Efficiency", axes=FALSE)
box()
axis(1, at=m)
axis(2)
lines(m, sig.5, type='b', pch=22)
lines(m, sig.05, type='b', pch=23)
lines(m, sig.0, type='b', pch=24)
legend('bottomleft', c("sigma=0", "sigma=0.05", "sigma=0.5",
"sigma=0.99"), pch=21:24)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
