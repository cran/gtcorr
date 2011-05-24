\name{gtcorr.hierarchical.user}
\alias{gtcorr.hierarchical.user}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the efficiency of hierarchical group testing procedures
for user specified arrangements}
\description{
\sQuote{gtcorr.hierarchical.user} calculates the efficiencies of
hierarchical group testing procedures, allowing for correlation between
units and test error. Cluster arrangements can are specified by the user
and probability of testing postive and correlation can vary by cluster.
}
\usage{
gtcorr.hierarchical.user(n, clusters, p, sigma=0, se=1, sp=1, model =
c("beta-binomial", "Morel-Neerchal", "Madsen"))
}
\arguments{
  \item{n}{
    a numeric vector of pool sizes where \code{n[s]} is the size of a pool in stage
    \code{s}.  The size of a pool in the last stage is 1, which can be omitted.
}
  \item{clusters}{
    a vector of length \code{n[1]} of integers from 1 up to the total
    number of clusters.  \code{cluster[i]} is the ith unit's cluster.
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
\details{
  Units are ordered such that the first pool in stage s contains units
  \code{1, ..., n[s]}, and the second pool in stage s has units
  \code{(n[s]+1),...,(2*n[s])} and so forth.
}
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
\code{\link{gtcorr.hierarchical}}
}
\examples{
##Calculate the efficiency of a two stage hierarchical procedure with
##n[1]=20, where the first cluster has 10 units and the second and third
##clusters have 5 units each.
n <- 20
clusters <- c(rep(1,10), rep(2, 5), rep(3, 5))
p <- .1
sigma <- .3
gtcorr.hierarchical.user(n, clusters, p, sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
