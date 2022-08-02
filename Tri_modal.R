
generate_trimodal <- function(mu, IC, outlier, p, corr, dist_mix) {
  
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
        rsigns = sample(c(-1,1),p, replace = T)
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
  y3 <- generate_data_mm(o=mu[3], s=corr[3], obsin=IC[3], obsout=outlier[3], p, dist=dist_mix[3]) 
  
  data <- rbind(y1[["data"]], y2[["data"]],y3[["data"]])   
  datain <- rbind(y1[["datain"]], y2[["datain"]],y3[["datain"]])   
  dataout <- rbind(y1[["dataout"]], y2[["dataout"]],y3[["dataout"]])  
  
  return(list(data=data, datain=datain, dataout=dataout))
}




Norm_45_5_cor0 <- generate_trimodal( mu= c(0, 5,10), IC = c(45, 45,45), outlier= c(5, 5,5), 
                                     p=100, corr= c(0,0,0), dist_mix= c("Normal", "Normal","Normal"))



dip.test(Norm_45_5_cor0$data)

library(diptest)

library(LaplacesDemon)

is.unimodal(Norm_45_5_cor0$data)

is.multimodal(Norm_45_5_cor0$data)

is.bimodal(Norm_45_5_cor0$data)

library(mousetrap)
bimodality_coefficient(Norm_45_5_cor0$data)


is.trimodal(Norm_45_5_cor0$data)

?is.trimodal



a11<-generate_data(obsin=45, obsout=5, p=100, s=0, dist="Normal")

is.unimodal(a11$data)