
source("utils.R")
source("tip_EM.R")

#data format: for n subjects observed at least once, a list of length n
#             each element collects the capture times of a subject in increasing order
#             between 0 and Tmax

set.seed(123)

Tmax=1
n=100
fake.data=vector("list", length=n)
for(i in 1:n)
  fake.data[[i]]=runif(min=0, max=Tmax, n=sample(1:3, size=1))

fake.cov=matrix(sample(0:1, size=n, replace=T), ncol=1)
#for p covariates, build a n x p matrix

#fitting (time taken: about one minute)
fake.fit = tipEM(obsData = fake.data, #list of captured subjects
                 XX = fake.cov,       #covariates, see below
                 TT = Tmax,           #time limit
                 CC = 2,              #number of latent classes for model coefficients
                 tol = 0.001,         #tolerance for convergence of EM algorithm
                 maxit = Inf,         #maximum number of EM iterations
                 verbose = TRUE,      #prints info during fitting process
                 inits = "det",       #initial parameter values, see below
                 infmet = "EM",       #inference method (currently, only EM is possible)
                 mx = 200,            #maximum number of iterations for parameter optimization at each EM step
                 met = "NM",          #currently fixed to "NM" for EM algorithm
                 o = TRUE,            #TRUE: model with covariates (Mo models)
                 te = TRUE,           #TRUE: model with time effect (Mt models)
                 b = TRUE,            #TRUE: model with behavioural effect (Mb models)
                 stdErr = FALSE)      #currently not implemented, leave F

#Note: XX cannot be empty. If you want to fit a model without covariates, add a random
#      covariate as XX and then set o =F

#Note: if you have no initial values to set, leave the default option
#      otherwise inits can either be the output of a previous tipEM(), or
#      a list with all coefficients given as alpha.v, beta.v, theta.v, mu.v, eta, piC, gamm
#      where piC are the initial latent class probabilities

