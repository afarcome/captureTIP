#UPDATED MAY 25, 2021

tipEM=
  function (obsData, XX, TT, CC, tol = 0.001, maxit = Inf, verbose = FALSE, 
            inits = "det", mx = 200, 
            o = TRUE, te = TRUE, b = TRUE, stdErr = FALSE) { 
    
    ###INITIAL CHECKS 
    
    if (all(sapply(obsData, length)) != TRUE) 
      stop("Only subjects with at least one observation must be included")
    if (max(sapply(obsData, max)) > TT) 
      stop("Some subjects observed after time horizon TT")
    if ((infmet == "EM" | infmet == "EMTMB") & stdErr) 
      stop("Standard Errors with EM not yet implemented. Fit EM with stdErr=FALSE and then use results as inits with a direct method")
    if ((infmet=="direct"|infmet=="directTMB"|infmet=="Laplace") & CC==1)
      stop("For 1 latent class, only EM and EMTMB are implemented")
    if (is.list(inits)) {
      if (!o && !all(inits$gamm == 0)) {
        warning("Requesting model without covariates but given non-zero initial regression coefficients")
      }
      if (!b && !all(c(inits$theta.v, inits$alpha.v) == 0)) {
        warning("Requesting no behavioural effect but non-zero theta and alpha given")
      }
      if (!te && inits$eta != 1) {
        warning("Requesting no time effects but initial eta given is not equal to unity")
      }
      if (any(alpha.v < 0)) 
        stop("Negative values for alpha are not admitted")
      if (any(beta.v < 0)) 
        stop("Negative values for beta are not admitted")
      if (any(theta.v < 0)) 
        stop("Negative values for theta are not admitted")
      if (eta < 0) 
        stop("Negative values for eta are not admitted")
      if (any(alpha.v > beta.v)) 
        stop("The value of alpha must be less than or equal to the value of beta")
      if (all(round(apply(matrix(piC,ncol=CC),1,sum),5)!=1)) 
        stop("The latent class probabilities must sum to 1")
      if (length(piC) != length(alpha.v)) 
        stop("The number of latent classes is inconsistent across initial values; check the length of piC, alpha.v, beta.v, theta.v, mu.v")
     }
    
    n = length(obsData)
    if (!is.matrix(XX)) {
        XX = as.matrix(XX)
    }
    if (nrow(XX) != n) 
      stop("nrow(XX) != number of subjects")

    ##SET INITIAL VALUES FOR PARAMETERS
    
    if (verbose) cat("Initializing the model...", fill=T)
    
    pp = ncol(XX)
    is.alt = ifelse(!is.list(inits) && (inits == "alt" || 
                                          inits == "detAlt"), TRUE, FALSE)
    use.clara = ifelse(!is.list(inits) && (inits == "detBig" || 
                                             inits == "detAlt"), TRUE, FALSE)
    sePar = seN = seN1 = seN2 = J = gr = NA
    

    #when initial values not given
    if (!is.list(inits)) {
      inits = list() 
      
      if (CC == 1) {
        inits$eta = 1
        inits$gamm = rep(0, pp)
        inits$alpha.v = ifelse(b, 0.5, 0)
        inits$beta.v = 1
        inits$mu.v = 0
        inits$theta.v = ifelse(b, 1, 0)
        inits$piC = 1
      }
      
      if (CC > 1) {
        d.inits = obsData
        len = max(sapply(d.inits, length))
        d.inits = cbind(0, t(sapply(unname(d.inits), 
                                    function(x) `length<-`(unname(unlist(x)), len))))
        if (use.clara) 
        {initial.class = cluster::clara(d.inits, CC, correct.d = FALSE)$clustering} else {
          d.inits1 = kml::clusterLongData(d.inits)
          try(invisible(capture.output(kml::kml(d.inits1, CC, parAlgo = kml::parALGO(saveFreq = Inf)))),
              silent=T)
          initial.class = try(unclass(kml::getClusters(d.inits1, CC)),silent=T)
          if(inherits(initial.class,"try-error")) 
            initial.class=sample(1:CC, size=length(obsData),replace=T)
          }
        
        inits$piC = as.numeric(prop.table(table(initial.class))) 

        for (cl in 1:CC) {
          if(length(which(initial.class==cl))==1)
            XX.cl=matrix(XX[initial.class == cl, ], nrow=1, ncol=pp) else
              XX.cl=as.matrix(XX[initial.class == cl, ])
          jnk = tipEM(obsData[initial.class == cl], XX.cl, TT, CC=1, tol = 0.001, 
                      maxit = Inf, verbose = F, inits = ifelse(use.clara, "detBig", "det"), 
                      o = o, te = te, b = b)
          if(verbose) cat(paste("Done for class", cl), fill=T)
          
          inits$alpha.v = c(inits$alpha.v, jnk$alpha)
          inits$beta.v = c(inits$beta.v, jnk$beta)
          inits$eta = c(inits$eta, jnk$eta)
          inits$gamm = rbind(inits$gamm, jnk$gamm)
          inits$mu.v = c(inits$mu.v, jnk$mu)
          inits$theta.v = c(inits$theta.v, jnk$theta)
        }
        inits$eta = max(mean(inits$eta), 0.1)
        inits$gamm = apply(inits$gamm, 2, mean)
        
        #to ensure identifiability
        or = order(inits$mu.v) 
        inits$mu.v = inits$mu.v[or]
        inits$alpha.v = inits$alpha.v[or]
        inits$beta.v = inits$beta.v[or]
        inits$theta.v = inits$theta.v[or]
        inits$piC = inits$piC[or]
        
        if (is.alt) {
          jnk2 = tipEM(obsData, XX, TT, 1, inits = ifelse(use.clara, "detBig", "det"), 
                       o = o, te = te, b = b)
          inits$gamm = jnk2$gamm
          inits$eta = jnk2$eta
        }
        if (inits$eta < 0.1) 
          inits$eta = 1
      }
      
      if (!o) {
        inits$gamm = rep(0, pp)
      }
      if (!te) {
        inits$eta = 1
      }
      if (!b) {
        inits$theta.v = rep(0, CC)
        inits$alpha.v = rep(0, CC)
      }
    }
    
    ##when starting values are given, or after the lines above
    if (is.list(inits)) {
      alpha.v = inits$alpha.v
      beta.v = inits$beta.v
      theta.v = inits$theta.v
      mu.v = inits$mu.v
      eta = inits$eta
      piC = inits$piC
      gamm = inits$gamm
      par.vec = unconstrain(alpha.v, beta.v, eta, theta.v, gamm, mu.v, o, te, b)
      ##NB the length of par.vec depends on the model
        #b=F no behavioural -> no alpha, no beta, no theta (each as many as CC)
        #o=F no covariates -> no gamma (as many as ncol(XX))
        #te=F no time effects -> no eta (1 parameter)
        #mu always present (as many as CC)
    }
    
    if (verbose) cat("Initialization completed", fill=T)
    if (verbose) cat("Starting parameter optimization...", fill=T)
    
   ###INITIAL COMPUTATION OF LOK LIKELIHOOD

    ex = expand.grid(1:CC, 1:n)
    ap = matrix(apply(ex, 1, 
                      function(jj) lik.unit.condC(alpha.v[jj[1]], beta.v[jj[1]], eta, theta.v[jj[1]], 
                                                  gamm, mu.v[jj[1]], XX[jj[2], ], obsData[[jj[2]]], TT)), 
                nrow = n, byrow = TRUE)

    sa = -log((1 - exp(t(t(-lambda.integral.0capture.mat(eta, theta.v, gamm, mu.v, XX, TT, CC)))))+1e-30)

    if (CC > 1) 
      sa = sa + matrix(log(piC), n, CC, byrow = TRUE)
    condlik.mat = ap + sa
    obj = apply(condlik.mat, 1, safe_sumlog) 
    log.WW = condlik.mat - obj
    obj = sum(obj) 
    WW = exp(log.WW)
    for(cc in 1:CC)
    {
    if(all(WW[,cc]==0)) WW[,cc]=1e-05
    }
    WW=WW/apply(WW,1,sum)
    
    ###PARAMETER OPTIMIZATION
    
    obj.old = obj - 10
    it = 0
    mx = ifelse(CC > 1, mx, 1000)
    
    
    while (abs(obj - obj.old) > tol && it < maxit){  
      
      it = it + 1
      obj.old = obj
      piC = apply(WW, 2, sum)
      piC = piC/sum(piC)
      piC[piC==0] =1e-05
      piC = piC/sum(piC)

      #EM INFERENCE METHOD
      
      par.vec.OLD=par.vec
      par.vec = try(optim(par.vec.OLD, function(x) wrap.lik(x,XX, pp, obsData, TT, CC, WW, piC, o, te, b), 
                          method = "BFGS", control = list(maxit = mx))$par, silent=T)
     if(inherits(par.vec,"try-error")) 
        {
        return(list(par.vec=par.vec.OLD, XX=XX, dati=obsData, piC=piC, WW=WW))
        stop("Issue with parameter optimization, see returned values to check the last parameter values
             and whether it is the global optimization or the initialization of a single latent class")
        }
      
      i = 0
      if (b) {
        alpha.v = exp(par.vec[i + 1:CC])
        i = i + CC
        beta.v = exp(par.vec[i + 1:CC]) + alpha.v
        i = i + CC
      }
      if (te) {
        eta = exp(par.vec[i + 1])
        i = i + 1
      }
      if (b) {
        theta.v = exp(par.vec[i + 1:CC])
        i = i + CC
      }
      if (o) {
        gamm = par.vec[i + 1:pp]
        i = i + pp
      }
      mu.v = par.vec[i + 1:CC]
      
      if(CC>1)
      {
      or = order(mu.v)
      alpha.v = alpha.v[or]
      beta.v = beta.v[or]
      theta.v = theta.v[or]
      mu.v = mu.v[or]
      piC = piC[or]
      par.vec=unconstrain(alpha.v, beta.v, eta, theta.v, gamm, mu.v, o, te, b)
      }
     
      ap = matrix(apply(ex, 1, 
                        function(jj) lik.unit.condC(alpha.v[jj[1]], beta.v[jj[1]], eta, theta.v[jj[1]], 
                                                    gamm, mu.v[jj[1]], XX[jj[2], ], obsData[[jj[2]]], TT)), 
                  nrow = n, byrow = TRUE)
      sa = -log((1 - exp(t(t(-lambda.integral.0capture.mat(eta, theta.v, gamm, mu.v, XX, TT, CC)))))+1e-30)
      
      if (CC > 1) 
        sa = sa + matrix(log(piC), n, CC, byrow = TRUE)
      condlik.mat.old=condlik.mat
      condlik.mat = ap + sa
      condlik.mat[is.na(condlik.mat)]=condlik.mat.old[is.na(condlik.mat)]
      condlik.mat[is.infinite(condlik.mat)]=condlik.mat.old[is.infinite(condlik.mat)]
      obj = apply(condlik.mat, 1, safe_sumlog)
      obj[is.na(obj)]=mean(obj,na.rm=T) 
      log.WW = condlik.mat - obj
      WW = exp(log.WW)
      for(cc in 1:CC)
      {
        if(all(WW[,cc]==0)) WW[,cc]=1e-05
      }
      WW=WW/apply(WW,1,sum)
      
      obj = sum(obj)
      if (verbose) {
        print(paste0("Iteration= ", it))
        print(paste0("obj= ", round(obj, 2), " diff= ", round(abs(obj - obj.old), 3)))
      }
    }
    
    if (verbose) cat("Done", fill=T)
    
    #ESTIMATION OF N
    
    ma = 1 - exp(-lambda.integral.0capture.mat(eta, theta.v, gamm, mu.v, XX, TT, CC))
    wm = apply(WW, 1, which.max)
    Nhat = sum(1/diag(ma[, wm]))
    
    #vett.condiz = sapply(1:nrow(XX), 
    #                     function(j) 1 - exp(-lambda.integral.0capture(eta, theta.v, gamm, mu.v, 
    #                                        matrix(XX[j, ], nrow = 1), TT, matrix(WW[j,], nrow = 1))))
    #Nhat2 = sum(1/vett.condiz)
    
    numpar=CC; #mu
      if(b) numpar=numpar+3*CC; #alpha, beta, theta
      if(te) numpar=numpar+1; #eta
      if(o) numpar=numpar+ncol(XX) #gamma
    
    if(verbose)
      print(paste("Nhat =", round(Nhat,0), "; AIC = ", round(-2*obj + 2*numpar,2)))
    # STORAGE
    
    
    return(list(lk = obj, alpha.v = alpha.v, beta.v = beta.v, 
                eta = eta, theta.v = theta.v, gamm = gamm, mu.v = mu.v, 
                piC = piC, Nhat = Nhat, #Nhat2 = Nhat2, 
                WW = WW, ma = ma,
                #sePar = sePar, seN = seN, seN1 = seN1, seN2 = seN2, J = J, gr = gr
                aic = round(-2*obj + 2*numpar,2)))
  }


