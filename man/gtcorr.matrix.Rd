\name{gtcorr.matrix}
\alias{gtcorr.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the efficiency of matrix group testing procedures}
\description{
\sQuote{gtcorr.matrix} calculates the efficiencies of matrix group
testing procedures, allowing for correlation between units and test error
}
\usage{
gtcorr.matrix(r, c, m = 1, p, sigma = 0, se = 1, sp = 1, r.prime, c.prime, arrangement = c("rectangular", "diagonal", "random"), model = c("beta-binomial", "Madsen", "Morel-Neerchal"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
    number of rows in the pooling matrix.
  }
  \item{c}{
    number of columns in the pooling matrix.
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
  \item{r.prime}{
    for a \sQuote{rectangular} arrangement, the number of rows in a
    rectangular cluster.
  }
  \item{c.prime}{
    for a \sQuote{rectangular} arrangement, the number of columns in a
    rectangular cluster.
  }
  \item{arrangement}{
    how clusters are arranged. Should be \sQuote{rectangular},
    \sQuote{diagonal} or \sQuote{random}.
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
  For a diagonal arrangement, \code{r}, \code{c}, and \code{m} should be
  equal.
  For a rectangular arrangement, \code{m} should be
  \code{r.prime*c.prime}.
  See Lendle et. al. 2011 for more information.  
}
\value{
  \item{r}{}
  \item{c}{}
  \item{m}{cluster size.}
  \item{r.prime}{number of rows in the pooling matrix.}
  \item{c.prime}{number of columns in the pooling matrix.}
  \item{param.grid}{a data frame containing the values of \code{p},
    \code{sigma}, \code{se}, and \code{sp} for each value of \code{efficiency}.}
  \item{arrangement}{arrangement.}
  \item{model}{model.}
  \item{efficiency}{a vector of efficiencies, one for each row of \code{param.grid}.}
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
##Plot efficiencies of a 16 by 16 matrix procedure by arrangement
sigma <- seq(0,.99, length.out=100)
sig2 <- seq(0, .99, length.out=10)
diag <- gtcorr.matrix(r=16, c=16, m=16, r.prime=1, c.prime=16,
                      arr='diag', p=.05, sigma=sigma)$efficiency
rand <- gtcorr.matrix(r=16, c=16, m=16, r.prime=1, c.prime=16,
                      arr='rand', p=.05, sigma=sig2)$efficiency 
rect1 <- gtcorr.matrix(r=16, c=16, m=16, r.prime=1, c.prime=16, p=.05,
                       sigma=sigma)$efficiency 
rect2 <- gtcorr.matrix(r=16, c=16, m=16, r.prime=2, c.prime=8, p=.05,
                       sigma=sigma)$efficiency 
rect3 <- gtcorr.matrix(r=16, c=16, m=16, r.prime=4, c.prime=4, p=.05,
                       sigma=sigma)$efficiency 

plot(sigma, diag, ylim=c(0, max(diag)), type='l', ylab="Efficiency", xlab="sigma")
lines(sig2, rand, col=2)
lines(sigma, rect3, col=3)
lines(sigma, rect2, col=4)
lines(sigma, rect1, col=5)
legend("bottomleft", c("Diagonal", "Random", "4x4 rect.", "2x8 rect.",
                       "1x16 rect."), lty=1, col=1:5) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
