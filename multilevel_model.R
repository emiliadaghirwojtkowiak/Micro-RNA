## Codes for multilevel Bayesian modeling ##

# In the description of methodology we use standard deviation while the codes applied in JAGS use #precision model
{
  for(i in 1:I){
    sETA [i,1:K] ~ dmnorm(vzeros[1:K], Omega.inv[1:K, 1:K])
  }
  
  for(i in 1:I){
    for(k in 1:K){
      ETA[i,k] <- xi[k]*sETA[i,k]
    }}
  
  for(n in 1:N){
    logC[i] <- miu[MIRNAID[n]] + ETA[ID[n], MIRNAID[n]] + betaDIS[MIRNAID[n]]*DIS[ID[n]] + betaAGE[MIRNAID[n]]*AGE[ID[n]]  + betaBW[MIRNAID[n]]*BW[ID[n]]
    Y[n] ~ dnorm(logC[n], tausigma[MIRNAID[n]]) 	
    Ycond[n] ~ dnorm(logC[n], tausigma[MIRNAID[n]])	
  }
  
  # missing values for DIS, AGE, BW
  for (i in 1:I) {
    DIS[i]~dbern(0.3333)
    AGE[i] ~ dnorm(0,1)
    BW[i] ~ dnorm(0,1)
  } 
  
  # priors 
  Omega.inv[1:K, 1:K] ~ dwish(Omegainvprior[1:K, 1:K], df)
  somega[1:K, 1:K] <- inverse(Omega.inv[1:K, 1:K])
  for(k in 1:K){
    for (k.prime in 1:K){
      rho[k,k.prime] <- somega [k,k.prime]/sqrt(somega [k,k]* somega[k.prime,k.prime])
    }
    omega[k] <- abs(xi[k])*sqrt(somega[k,k])	# scaled SD of parameter k
  }
  
  df <- K + 1
  
  for (k in 1:K ) {
    miu[k]  ~ dnorm(0, 0.04) # 1/5^2
    sigma[k] ~ dnorm(0, 1)T(0,)
    tausigma[k] <- 1/(sigma[k]*sigma[k])
    
    xi[k] ~ dnorm(0, 1)T(0,)
    
    betaDIS[k] ~ dnorm(0, 1)    
    betaAGE[k]  ~ dnorm(0, 1)
    betaBW[k]    ~ dnorm(0, 1)
    
    d[k] <- betaDIS[k]/pow(sigma[k]^2+omega[k]^2,0.5)
    
  }
}
