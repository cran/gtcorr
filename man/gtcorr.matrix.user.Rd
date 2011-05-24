\name{gtcorr.matrix.user}
\alias{gtcorr.matrix.user}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the efficiency of matrix group testing procedures for
  user specified arrangements
}
\description{
\sQuote{gtcorr.matrix.user} calculates the efficiencies of
matrix based group testing procedures, allowing for correlation between
units and test error. Cluster arrangements can are specified by the user
and probability of testing postive and correlation can vary by cluster.
}
\usage{
gtcorr.matrix.user(clusters, p, sigma=0, se=1, sp=1, model = c("beta-binomial", "Madsen", "Morel-Neerchal"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
    a matrix of integers from 1 up to the total number of clusters.
    \code{cluster[i,j]} is the ith, jth unit's cluster.
}
  \item{p}{
    probability of a unit testing positive.  If the length is one, then
    all clusters have the same probability of testing positive.  If the
    length is the total number of clusters, then \code{p[k]} is the
    probability that a unit in the kth cluster tests positive. 
}
  \item{sigma}{
    pairwise correlation of two units in a cluster. If the length is one, then
    all clusters have the same pairwise correlation.  If the length is
    the total number of clusters, then \code{sigma[k]} is the pairwise
    correlation for the kth cluster.
}
  \item{se}{
    sensitivity.  The probability that a pool of units tests positive given
    than at least one unit in that pool is positive
}
  \item{sp}{
    specificity.  The probability that a pool of units tests negative
    given that at least one unit in that pool is negative 
}
  \item{model}{
    probability model for clusters. Should be \sQuote{beta-binomial},
    \sQuote{Madsen}, or \sQuote{Morel-Neerchal}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 A length one vector containing the efficiency is returned
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

\seealso{
\code{\link{gtcorr.matrix}}
}
\examples{
##Calculate the efficiency of a 4 by 4 matrix procedure where the first
##two rows are in the first cluster, the bottom left 2 by 2 submatrix
##makes up the second cluster, and the bottom right 2 by 2 submatrix
##makes up the third cluster.
clusters <- matrix(NA, 4, 4)
clusters[1:2, ] <- 1
clusters[3:4, 1:2] <- 2
clusters[3:4, 3:4] <- 3
p <- .01
sigma <- .4
gtcorr.matrix.user(clusters, p, sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
