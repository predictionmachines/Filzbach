\name{filzbach-package}
\alias{filzbach-package}
\alias{filzbach}
\docType{package}
\title{
Filzbach MCMC sampler
}
\description{
Filzbach is a flexible, fast, robust, parameter estimation engine that allows you to parameterize arbitrary, 
	non-linear models, of the kind that are necessary in biological sciences, against multiple, heterogeneous data sets. 
	Filzbach allows for Bayesian parameter estimation, maximum likelihood analysis, priors, latents, 
	hierarchies, error propagation and model selection, from just a few lines of code.
}
\details{
\tabular{ll}{
Package: \tab filzbach\cr
Type: \tab Package\cr
Version: \tab 0.7-1\cr
Date: \tab 2015-07-04\cr
License: \tab MIT + file LICENSE\cr
}

Index:
\preformatted{
runMCMC                 Generates a sample from target distribution
                        function
filzbach-package        Filzbach MCMC
}

}
\author{
Vassily Lyutsarev [aut, cre],
  Drew Purves [aut] (creator of the original filzbach C code)

Maintainer: Vassily Lyutsarev<vassilyl@microsoft.com>
}
\references{
\url{http://research.microsoft.com/filzbach}
}
\keyword{ package mcmc }
\seealso{
\code{\link{runMCMC}}
~~ Optional links to other man pages, e.g. ~~
}
\examples{
## NORMAL example
#################

y<-rnorm(999, 37, 9)
iter<-200#00 # use bigger value to get better results
r<-runMCMC(iter, iter,
    function(mean, sigma)
        sum(log(dnorm(y,mean,sigma))),
    999,
    list(
        mean=c(0.0010,2000.0,100.0,0,0,1), 
        sigma=c(0.0010,2000.0,20.0,1,0,1)
    ))

print(r)
print(summary(r))
}
