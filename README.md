R demo and functions implementing the approach described in the paper
"Continuous time-interaction processes for population size estimation, 
with an application to drug dealing in Italy"

The files are:

-- timEM.R: it contains function tipEM, implementing the Expectation-Maximization algorithm for estimating a Time Interaction Process by maximizing the conditional likelihood. Inputs are:

obsData: list of length n, each element is capture times for i-th individual (see demo.R) 
XX: matrix of covariates (see demo.R)
TT: scalar with observation time span (must be larger than or equal to the last capture time)
CC: number of latent classes

tol: tolerance for convergence of EM, defaults to 0.001

maxit: maximum number of EM iterations, defaults to Inf

verbose: whether print progress information, defaults to FALSE,

inits: "det" (default) for deterministic initial solution,
list with user specified inits for parameters (see output for names), or
"alt" for alternative deterministic initial solution

mx: maximum number of iterations for numerical optimization algorithm at the M step, defaults to 200 

o, te, b: whether to include observed covariates, time-effects, behavioural effects (all default to TRUE)

stdErr: whether to compute standard errors (defaults to FALSE)

Output is a list with elements:

lk: the log-likelihood at convergence
alpha.v, beta.v, eta, theta.v, gamm, mu.v: parameter estimates (in the notation of the paper)
piC: parameter estimates for the prior latent mass probabilities
Nhat: population size estimate
seN: standard error for population size estimate (NA if stdErr=FALSE) 
J: observed information matrix 
aic: AIC 

-- demo.R: generate a data set and estimate parameters and population size 

-- utils.R: internal functions 

Suitable comments are directly included in the demo script in order to illustrate use of the functions in R.

The library(statmod) must be installed in R for the scripts to work. 