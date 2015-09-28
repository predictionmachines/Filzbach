\name{runMCMC}
\alias{runMCMC}
\title{
Generates a sample from target distribution function
}
\description{
Run Metropolis-Hastings algorithm to sample parameter space.
}
\usage{
runMCMC(burnin, eststeps, loglikelihood, samplesize, paramdefs, thinning = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{burnin}{
number of iterations in burn-in phase
%%     ~~Describe \code{burnin} here~~
}
  \item{eststeps}{
number of iterations in a phase that collects samples (~ number of samples times 100, where 100 is the default value for thining).
}
  \item{loglikelihood}{
a function that computes a logarithm of un-normalized likelihoodgiven the model parameter values (see examples below).
}
  \item{samplesize}{
an integer number that is used in BIC computation. Usually, it corresponds to the number of terms in the likelihood.
}
  \item{paramdefs}{
a list of parameter definitions. 
      The number and position of parameter definitions must correspond to parameters in the loglikelihood functions.
	  Each list item contains data for parameter_create function in the original Filzbach package.
}
  \item{thinning}{
minimum number of iteration that produce one sample in the resulting table.
}
}
\value{
returns a matrix of  with columns corresponding to parameters.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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


## SEED_RAIN example
###################

phref <- 7.0
depthref <- 10.0
mossref <- 0.0

lambda <- function (m0, beta, gph, gdepth, gmoss, depth, ph, moss)
    beta + m0*exp(gph*(ph-phref) + gdepth*(depth-depthref) + gmoss*(moss-mossref))

depth <- runif(999,0,20)
ph <- 5.0 + 2.0*(depth/20.0) + runif(999,0,2)
moss <-  as.integer(runif(999)<0.5)
count <- sapply(lambda(2.0, 0.010, 0.40, 0.00, 0.20, depth, ph, moss), function(lam) rpois(1,lam))
iter<-200#00 # use bigger value to get better results

r<-runMCMC(iter, iter,
    function(m0, beta, gph, gdepth, gmoss)
    {
        s <- 0
        for (i in 1:999) {
            p <- dpois(count[i],lambda(m0, beta, gph, gdepth, gmoss, depth[i], ph[i], moss[i]))
            if (is.na(p) || p<1e-8) p <- 1e-8
            s <- s+log(p)
        }
        s
    },
    999,
    list(
        m0 = c(0.000010,50.0,5.0,1,0,1),
        beta = c(0.000010,50.0,5.0,1,0,1),
        gph = c(-3.0,3.0,0.10,0,0,1),
        gdepth = c(-1.0,1.0,0.10,0,0,1),
        gmoss = c(-3.0,3.0,0.10,0,0,1)
    ))

print(summary(r))

## POISSON_MULTISPP sample
##########################

numsp <- 10
true_lambda <- runif(numsp, 1.0, 1.9)
print(true_lambda)

spid <- as.integer(runif(999,0, numsp)) + 1
count <- as.numeric(sapply(spid, function(i){rpois(1, true_lambda[i])}))
iter<-200#00 # use bigger value to get better results

r<-runMCMC(iter, iter,
function(lambda, lambda_mean, lambda_sdev)
{
    s <- 0
    for (i in 1:999) {
        s <- s + dpois(count[i], lambda[spid[i]], log=TRUE)
    }
    for (i in 1:numsp) {
        s <- s + dnorm(lambda[i], lambda_mean, lambda_sdev, log=TRUE)
    }
    s
},
999,
list(
    "lambda" = c(0.0010, 200.0, 100.0, 1, 0, 1, numsp),
    "lambda_mean" = c(0.0010, 200.0, 100.0, 1, 0, 1),
    "lambda_sdev" = c(0.0010, 2.0, 0.20, 1, 0, 1)
))

print(summary(r))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
