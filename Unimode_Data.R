#############################Unimodal###########################
generate_data <- function(obsin, obsout, p, s, dist) {
  if (dist=="Normal") {
    gennorm = function(obsin, p, s, obsout, deltap = 0.023){
      require(mvnfast)
      require(mvtnorm)
      mu0 = (rep(0, p))
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvnorm(1-deltap, sigma  = sigma)$quantile
      mu1a = rep(mu1a, p) 
      rsigns = sample(c(-1,1),p, replace = T)
      mu1 = mu1a * rsigns  ##  Shifting Random Dimensions
      #mu1 = mu1a 
      din = rmvn(obsin, mu0, sigma)
      if (obsout > 0){
        dout = rmvn(obsout, mu1, sigma)
      }
      else
      {
        dout = c()
      }
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p))
    }
    data <- gennorm(obsin,p,s,obsout)
    
  } else if (dist=="Lognormal") { 
    genlognorm = function(obsin, p, s, obsout,  deltap = 0.023){
      require(mvnfast)
      require(mvtnorm)
      mu0 = (rep(0, p))
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvnorm(1-deltap, sigma  = sigma)$quantile
      mu1a = rep(mu1a, p) 
      rsigns = sample(c(-1,1),p, replace = T)
      mu1 = mu1a * rsigns  ## Shifting Random Dimensions
      # mu1 = mu1a 
      din = exp(rmvn(obsin, mu0, sigma))
      if (obsout > 0){
        dout = exp(rmvn(obsout, mu1, sigma))
      }
      else
      {
        dout = c()
      }
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p))
    }
    data <-  genlognorm(obsin,p,s,obsout)
    
  } else if (dist=="T"){
    gent = function(obsin, p, s, obsout, deltap = 0.023){
      require(mvnfast)
      require(mvtnorm)
      mu0 = (rep(0, p))
      sigma <-matrix(s, p, p) + diag((1-s), p)
      mu1a = qmvt(1-deltap, sigma  = sigma, df = 10)$quantile
      mu1a = rep(mu1a, p) 
      rsigns = sample(c(-1,1),p, replace = T)
      mu1 = mu1a * rsigns  ## Shifting Random Dimensions
      # mu1 = mu1a 
      din = mvnfast::rmvt(obsin, mu0,sigma, df = 10)
      if (obsout > 0){
        dout = mvnfast::rmvt(obsout, mu0,sigma, df = 10) + matrix(rep(mu1,obsout), ncol = p, byrow = T)
      }
      else
      {
        dout = c()
      }
      return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p))  
    }
    data <- gent(obsin,p,s,obsout)
  }
  return(data)
}


N_mu_unimode_s45_5_cor0<-NULL

for(i in 1:1000){
  print(i)
  N_mu_unimode_s45_5_cor0[[i]] <- generate_data(obsin=45, obsout=5, p=100, s=0, dist="Normal")
}

saveRDS(N_mu_unimode_s45_5_cor0,file="E:\\Miami\\Research\\Professors\\Data\\Unimodal\\Outlier\\N_mu_unimode_s45_5_cor0.rds")

###

N_mu_unimode_s45_5_cor75<-NULL

for(i in 1:1000){
  print(i)
  N_mu_unimode_s45_5_cor75[[i]] <- generate_data(obsin=45, obsout=5, p=100, s=0.75, dist="Normal")
}

saveRDS(N_mu_unimode_s45_5_cor75,file="E:\\Miami\\Research\\Professors\\Data\\Unimodal\\Outlier\\N_mu_unimode_s45_5_cor75.rds")

###
LN_mu_unimode_s45_5_cor0<-NULL

for(i in 1:1000){
  print(i)
  LN_mu_unimode_s45_5_cor0[[i]] <- generate_data(obsin=45, obsout=5, p=100, s=0, dist="Lognormal")
}


saveRDS(LN_mu_unimode_s45_5_cor0,file="E:\\Miami\\Research\\Professors\\Data\\Unimodal\\Outlier\\LN_mu_unimode_s45_5_cor0.rds")

####
LN_mu_unimode_s45_5_cor75<-NULL

for(i in 1:1000){
  print(i)
  LN_mu_unimode_s45_5_cor75[[i]] <- generate_data(obsin=45, obsout=5, p=100, s=0.75, dist="Lognormal")
}


saveRDS(LN_mu_unimode_s45_5_cor75,file="E:\\Miami\\Research\\Professors\\Data\\Unimodal\\Outlier\\LN_mu_unimode_s45_5_cor75.rds")

