generate_bimodal <- function(mu, IC, outlier, p, corr, dist_mix) {
  
  generate_data_mm <- function(o, obsin, obsout, p, s, dist) {
    
    if (dist=="Normal") {
      gennorm = function(obsin, p, s, obsout, deltap = 0.023){
        require(mvnfast)
        require(mvtnorm)
        mu0 = (rep(o, p))
        sigma <-matrix(s, p, p) + diag((1-s), p)
        mu1a = qmvnorm(1-deltap, sigma = sigma)$quantile
        mu1a = rep(mu1a, p) 
        rsigns = sample(c(-1,1),p, replace = T)
        mu1 = mu1a * rsigns  ##  Shifting Random Dimensions
        #mu1 = mu1a 
        din = rmvn(obsin, mu0, sigma)
        if (obsout > 0){
          dout = rmvn(obsout, mu0, sigma)+matrix(rep(mu1,obsout), ncol = p, byrow = T)
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
        mu0 = (rep(mu[1], p))
        sigma <-matrix(s, p, p) + diag((1-s), p)
        mu1a = qmvnorm(1-deltap, sigma  = sigma)$quantile 
        mu1a = rep(mu1a, p) 
        rsigns = sample(c(1,1),p, replace = T)
        mu1 = mu1a * rsigns  ## Shifting Random Dimensions
        
        if (o == mu[1]){        
          din = exp(rmvn(obsin, mu0, sigma)) 
          if (obsout > 0){
            dout = exp(rmvn(obsout, mu1, sigma))
          } 
          else
          {
            dout = c()}}
        else if (o == mu[2]){        
          din = exp(rmvn(obsin, mu0, sigma))+mu[2]
        if (obsout > 0){
          dout = exp(rmvn(obsout, mu1, sigma))+mu[2]
          }

        else
          {
  dout = c()
}}
        return(list(datain = din, dataout = dout, data = rbind(din,dout), dim_shifted = p))
      }
      data <-  genlognorm(obsin,p,s,obsout)
      
    } else if (dist=="T"){
      gent = function(obsin, p, s, obsout, deltap = 0.023){
        require(mvnfast)
        require(mvtnorm)
        mu0 = (rep(o, p))
        sigma <-matrix(s, p, p) + diag((1-s), p)
        mu1a = qmvt(1-deltap, sigma  = sigma, df = 10)$quantile
        mu1a = rep(mu1a, p) 
        rsigns = sample(c(-1,1),p, replace = T)
        mu1 = mu1a * rsigns  ## Shifting Random Dimensions
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
  
  y1 <- generate_data_mm(o=mu[1], s=corr[1], obsin=IC[1], obsout=outlier[1], p, dist=dist_mix[1])  
  y2 <- generate_data_mm(o=mu[2], s=corr[2], obsin=IC[2], obsout=outlier[2], p, dist=dist_mix[2])  
  
  data <- rbind(y1[["data"]], y2[["data"]])   
  datain <- rbind(y1[["datain"]], y2[["datain"]])   
  dataout <- rbind(y1[["dataout"]], y2[["dataout"]])  
  
  return(list(data=data, datain=datain, dataout=dataout))
}

#test run
td <- generate_bimodal( mu= c(0, 10), IC = c(8, 8), outlier= c(2, 2), 
                        p=2, corr= c(0,0), dist_mix= c("T", "T"))
b <- OCPtextunt(td$data)
d<-data.frame(b$RKD_OCP, c(rep(1,10), rep(2,10)))
boxplot(d$b.RKD_OCP~d$c.rep.1..10...rep.2..10..)
