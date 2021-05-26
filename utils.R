library(statmod)
dd=gauss.quad(39,"chebyshev1")
fint=function(x,eta,theta,lower,upper) {
  cv=(upper-lower)*x/2+(lower+upper)/2
  cv^(eta+1)*exp(theta*eta*cv^(eta-1))*sqrt(1-x^2)*(upper-lower)/2}

####################################################################
#1) lambda(t)

lambda=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,tt){ 
  lam=alpha*sum(exp(-beta*(tt-individual.hist[individual.hist<tt])))+
    eta*(tt^(eta-1))*exp((XX%*%gamm)+mu)*exp(-theta*(sum(individual.hist<=tt)-eta*(tt^(eta-1))))
  return(lam)
}

####################################################################
#2) integral of lambda in [0,T]

lambda.integral=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT){
  #self-exciting
  somma.sep=-alpha/beta*sum(exp(-beta*(TT-individual.hist))-1)
  #self-correcting
  fuori=exp((XX%*%gamm)+mu)
  tt2=c(0,individual.hist,TT)
  vett.int=(sapply(1:(length(tt2)-1), function(j) sum(fint(dd$nodes,eta,theta,tt2[j],tt2[j+1])*dd$weights))) ##as.numeric?
  somma.scp=fuori*sum(exp(-theta*0:length(individual.hist))*vett.int)
  return(as.numeric(somma.sep+somma.scp))
}

####################################################################
#4) integral of lambda for unobserved subjects

lambda.integral.0capture=function(eta,theta.v,gamm,mu.v,XX,TT,WW){
  XX=as.matrix(XX)
  fuori=sapply(1:length(mu.v), function(j) exp((XX%*%gamm)+mu.v[j]))
  integr=sapply(1:length(theta.v), 
                function(j) sum(dd$weights*fint(dd$nodes,eta,theta.v[j],0,TT)))
  if(!all(dim(fuori)==dim(WW))) warning("Problem with function *lambda.integral*")
  out=as.numeric((fuori*WW)%*%integr) 
  return(out)
}

####################################################################
#4b) integral of lambda for unobserved subjects - matrix for all latent classes

lambda.integral.0capture.mat=function(eta,theta.v,gamm,mu.v,XX,TT,CC){
  XX=as.matrix(XX)
  fuori=sapply(1:CC, function(j) exp((XX%*%gamm)+mu.v[j]))
  integr=sapply(1:CC, function(j) sum(fint(dd$nodes,eta,theta.v[j],0,TT)*dd$weights))
  t(integr*t(matrix(fuori,ncol=CC)))
}


####################################################################
#5a) likelihood conditional to class c - one subject

lik.unit.condC=function(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT){
  lik=0
  for(i in 1:length(individual.hist))
    lik=lik+log(lambda(alpha,beta,eta,theta,gamm,mu,XX,individual.hist[1:i],individual.hist[i]))
  lik=lik-lambda.integral(alpha,beta,eta,theta,gamm,mu,XX,individual.hist,TT)
  return(lik) 
}

####################################################################
#5b) unconditional likelihood - one subject

lik.unit.compl=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,individual.hist,TT,WW,pi.C){
  lik.v=(sapply(1:length(alpha.v), function(j) lik.unit.condC(alpha.v[j],beta.v[j],eta,theta.v[j],
                                                              gamm,mu.v[j],XX,individual.hist,TT)))
  lik.v*(WW)+(WW)*log(pi.C)
}

####################################################################
#6) opposite of likelihood - all subjects

lik.all=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C){
  sum(unlist(lapply(1:length(obs.data), function(j) -lik.unit.compl(alpha.v,beta.v,eta,
                                                                    theta.v,gamm,mu.v,XX[j,],unlist(obs.data[[j]]),TT,WW[j,],pi.C))))}

####################################################################
#7) opposite of likelihood conditional on class c - all subjects

cond.lik.all=function(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C) {
  lik.all(alpha.v,beta.v,eta,theta.v,gamm,mu.v,XX,obs.data,TT,WW,pi.C)+
    log(1-exp(-lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,WW)))
}

####################################################################
#8) function for likelihood maximization

wrap.lik=function(par.vec, XX, pp, obs.data, TT, CC, WW, pi.C, o, te, b) {
  if(b) alpha.v=exp(par.vec[1:CC]) else alpha.v=rep(0,CC)
  if(b) beta.v=exp(par.vec[CC+1:CC])+alpha.v else beta.v=rep(1,CC)
  if(b & te) eta=exp(par.vec[2*CC+1]) else if(!te) eta=1 else if(!b & te) eta=exp(par.vec[1])
  if(b & te) theta.v=exp(par.vec[(2*CC+1)+1:CC]) else if(b & !te) theta.v=exp(par.vec[(2*CC)+1:CC]) else if(!b) theta.v=rep(0,CC)
  if(!o) gamm=rep(0, pp) else gamm=par.vec[(length(par.vec)-(CC-1))-pp:1]
  mu.v=par.vec[(length(par.vec)-(CC-1)):length(par.vec)]
  
  sum(unlist(lapply(1:length(obs.data), function(j) 
    -lik.unit.compl(alpha.v,beta.v,eta, theta.v,gamm,mu.v,matrix(XX[j,], nrow = 1),obs.data[[j]],
                    TT,matrix(WW[j,], nrow = 1),pi.C))))+
    sum((log(1-exp(-lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,
                    WW)[!is.na(lambda.integral.0capture(eta,theta.v,gamm,mu.v,XX,TT,WW))])+1e-30)))
}



#####################################################################
#OTHER FUNCTIONS SUPPORTING TIP

safe_sumlog = function(x) {

s <- (try(snipEM::sumlog(x), silent=T))
if(inherits(s,"try-error")) {
s <- max(x)}

return(s)}

unconstrain = function(alpha.v, beta.v, eta, theta.v, gamm, mu.v, o, te, b) {
  par.vec = c()

  if (b) par.vec = c(par.vec, log(alpha.v), log(beta.v - alpha.v))
  if (te) par.vec = c(par.vec, log(eta))
  if (b) par.vec = c(par.vec, log(theta.v))
  if (o) par.vec = c(par.vec, gamm)

  par.vec = c(par.vec, mu.v)

  return(par.vec)
}

unconstrainPiC = function(alpha.v, beta.v, eta, theta.v, gamm, mu.v, piC, o, te, b) {
  par.vec = c()

  if (b) par.vec = c(par.vec, log(alpha.v), log(beta.v - alpha.v))
  if (te) par.vec = c(par.vec, log(eta))
  if (b) par.vec = c(par.vec, log(theta.v))
  if (o) par.vec = c(par.vec, gamm)

  par.vec = c(par.vec, mu.v, log(piC[-1] / piC[1]))

  return(par.vec)
}

stDirectVec=function(par.vec,CC,XX,TT,obsData,b,te,o,inits) {

   pp = ncol(XX)
   n  = nrow(XX)
    i = 0

    if(!b) {alpha=inits$alpha
            beta=inits$beta
	    theta=inits$theta}

    if(!te) {eta = inits$eta}

    if(!o) {gamm = inits$gamm}

    if (b) {
      alpha = exp(par.vec[i + 1:CC])
      i = i + CC
      beta = exp(par.vec[i + 1:CC]) + alpha
      i = i + CC
    }

    if (te) {
      eta = exp(par.vec[i+1])
      i = i + 1
    }

    if(b) {
      theta = exp(par.vec[i + 1:CC])
      i = i + CC
    }

    if (o) {
      gamm = par.vec[i + 1:pp]
      i = i + pp
    }

    mu = par.vec[i + 1:CC]
    i = i + CC

    piC = c(1, exp(par.vec[i + 1:(CC - 1)]))
    piC = piC / sum(piC)

    or = order(mu)
    piC = matrix(piC[or], n, CC, byrow = TRUE)
    mu = mu[or]
    alpha = alpha[or]
    beta = beta[or]
    theta = theta[or]

obj=list()
obj$alpha=alpha[or]
obj$beta=beta[or]
obj$eta=eta
obj$theta=theta[or]
obj$gamm=gamm
obj$mu=mu[or]
obj$piC=piC

stDirect(obj, XX, TT, CC, obsData)[1]}

stDirect = function(object, XX, TT, CC, obsData) {
  n = nrow(XX)
  XX = as.matrix(XX)
  pp = ncol(XX)
  ex = expand.grid(1:CC, 1:n)
  WW = matrix(1 / CC, n, CC)
  eta = object$eta
  piC = object$piC
  theta.v = object$theta
  gamm = object$gamm
  mu.v = object$mu
  alpha.v = object$alpha
  beta.v = object$beta
  ma = 1 - exp(-lambda_integral_0capture_mat(eta, theta.v, gamm, mu.v, XX, TT, CC))

  ap = matrix(apply(ex, 1, function(jj) lik_unit_condC(alpha.v[jj[1]], beta.v[jj[1]], eta, theta.v[jj[1]], gamm, mu.v[jj[1]], XX[jj[2], ],
                                                       obsData[[jj[2]]], TT)), nrow = n, byrow = TRUE)

  sa = -log(1 - exp(t(t(-lambda_integral_0capture_mat(eta, theta.v, gamm, mu.v, XX, TT, CC)))))

  if (CC > 1)
    sa = sa + matrix(log(piC), n, CC, byrow = TRUE)

  if (CC == 1)
    sa = sa + log(piC)

  condlik.mat = ap + sa

  obj = apply(condlik.mat, 1, safe_sumlog)
  log.WW = condlik.mat - obj
  WW = exp(log.WW)
  wm = apply(WW, 1, which.max)
  Nhat = sum(1 / diag(ma[, wm]))
  vett.condiz = sapply(1:nrow(XX), function(j) 1 - exp(-lambda_integral_0capture(eta, theta.v, gamm, mu.v, matrix(XX[j, ], nrow = 1), TT, matrix(WW[j, ], nrow = 1))))

  Nhat2 = sum(1 / vett.condiz)

  return(c(Nhat, Nhat2))
}
