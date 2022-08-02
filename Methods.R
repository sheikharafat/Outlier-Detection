##DBSCAN
datafnn_dbscan<- function (outliers, datag,n){
  
  require(dbscan)
  res <- NULL
  flags <- NULL
  
  for(i in 1:n){
    res[[i]] <- dbscan(datag[[i]]$data,eps = 0.5, minPts = 5)
    flags[[i]]<-as.logical(!res[[i]][["cluster"]])
  }
  
  simdat <-data.frame( Flag=unlist(flags), Outlier=outliers)
  # Please check the markers
  simdat$Status <- ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="In", "FP",
                          ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="Out", "TP",
                                 ifelse(simdat$Flag=="FALSE" & simdat$Outlier=="In", "TN","FN")))
  return(simdat)
}


DB_N_mu_00_s4505_4505_cor0_final_p100<- datafnn_dbscan(outliers=ICOC_4505out,datag=N_mu_00_s4505_4505_cor0_final_p100,n=1000)

table(DB_N_mu_00_s4505_4505_cor0_final_p100$Status)

Detection_Rate<- (length(which(DB_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))/(length(which(DB_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(DB_N_mu_00_s4505_4505_cor0_final_p100$Status=='FN'))))*100

Detection_Rate

Correctly_Classified<-((length(which(DB_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(DB_N_mu_00_s4505_4505_cor0_final_p100$Status=='TN')))/length(DB_N_mu_00_s4505_4505_cor0_final_p100$Status))*100

Correctly_Classified


## ISOFORREST

datafnn_isotree <- function (outliers, datag,n,h){
  
  require(isotree)
  iso <- NULL
  pred<- NULL
  flags <- NULL
  
  for(i in 1:n){
    iso[[i]] <- isolation.forest(datag[[i]]$data,ndim=1, ntrees=100,missing_action="fail",nthreads=1)
    pred[[i]] <-predict(iso[[i]],datag[[i]]$data)
    flags[[i]]<-pred[[i]]>h
  }
  
  simdat <-data.frame( Flag=unlist(flags), Outlier=outliers)
  # Please check the markers
  simdat$Status <- ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="In", "FP",
                          ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="Out", "TP",
                                 ifelse(simdat$Flag=="FALSE" & simdat$Outlier=="In", "TN","FN")))
  return(simdat)
}

ISOTREE_N_mu_00_s4505_4505_cor0_final_p100<- datafnn_isotree(outliers=ICOC_4505out,datag=N_mu_00_s4505_4505_cor0_final_p100, h=0.5,n=1000)

table(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status)


Detection_Rate<- (length(which(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))/(length(which(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status=='FN'))))*100

Detection_Rate

Correctly_Classified<-((length(which(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status=='TN')))/length(ISOTREE_N_mu_00_s4505_4505_cor0_final_p100$Status))*100

Correctly_Classified


### New

OCPtextunt = function(dataframe, degree = 2, peeled_to_obs = 2, h = NULL, standardize = TRUE, na.rm=TRUE, plot = FALSE, nu=0.0001){
  p = dim(dataframe)[2]
  require(kernlab)
  if(standardize == TRUE){
    dataframe = data.frame(scale(dataframe, center = TRUE, scale = TRUE))
  }
  df_oc = data.frame(type=1,dataframe)
  N = dim(df_oc)[1]
  m = dim(df_oc)[1]
  invsigma = 1/p
  d = 0
  o1 = 0
  bagdat = df_oc;
  KD_OCP = c()
  while (m > peeled_to_obs)
  {
    o1 = o1+1;
    d = d+1 ;
    if (class(bagdat[,-1]) =="numeric"){mu_OCP=mean(bagdat[,-1])} else {mu_OCP = colMeans(bagdat[,-1])}
    w = ksvm(type ~., data=bagdat, type="one-svc", kernel="rbfdot", nu=nu, scaled = F, kpar=list(sigma=invsigma)) #Find boundary
    
    dec_bound_dist = predict(w, dataframe, type = 'decision')*-1 #negative if inside (modified)
    # dec_bound_dist[dec_bound_dist < 0] = 0
    KD_OCP = cbind(KD_OCP, dec_bound_dist) #grow by the peel, more negative~closer to the center
    
    sv = SVindex(w);
    bagdat = bagdat[-sv,]
    m = dim(bagdat)[1]
  }
  
  KD_OCP = rowMeans(KD_OCP) #na.rm=TRUE?
  MAD = median(abs(KD_OCP-median(KD_OCP)))
  RKD_OCP = (KD_OCP  - median(KD_OCP))/MAD
  
  if(is.null(h)) {h <- quantile(RKD_OCP, 0.75) + 1.5 * IQR(RKD_OCP)}
  OCP_Flag = RKD_OCP > h
  if (plot == TRUE){
    plot (RKD_OCP, type = 'l', col = 'grey', lwd = 2,cex.main=1,  xlab="Observation Index", ylab="sRKDs")
    par(new = TRUE)
    plot (RKD_OCP, col = ifelse(OCP_Flag > 0,'red', NA),  xaxt='n', ann=FALSE)
    par(new = TRUE)
    abline (h = h, lty=3, col = 'red', lwd = 1.8)
  }
  return(list(RKD_OCP= RKD_OCP, OCP_Flag = OCP_Flag))
}

##############################################################################################

#                             Performance Evaluation Function                                #

##############################################################################################
datafnn_new<- function (outliers, datag, fn,n){
  
  distances <- NULL
  flags <- NULL
  
  for(i in 1:n){
    print(i)
    distances[[i]] <- fn(datag[[i]]$data)$RKD_OCP
    flags[[i]] <-fn(datag[[i]]$data)$OCP_Flag
  }
  
  simdat <-data.frame(Distance=unlist(distances), Flag=unlist(flags), Outlier=outliers)
  # Please check the markers
  simdat$Status <- ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="In", "FP",
                          ifelse(simdat$Flag=="TRUE" & simdat$Outlier=="Out", "TP",
                                 ifelse(simdat$Flag=="FALSE" & simdat$Outlier=="In", "TN","FN")))
  return(simdat)
}

# FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100<- datafnn_new(outliers=ICOC_4505out,datag=N_mu_00_s4505_4505_cor0_final_p100, fn=OCPtextunt,n=1000)
# 
# saveRDS(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100,file="D:\\Outlier\\Performance\\FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100.rds")


FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100<-readRDS("D:\\Outlier\\Performance\\FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100.rds")

table(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status)

Detection_Rate<- (length(which(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))/(length(which(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status=='FN'))))*100

Detection_Rate

Correctly_Classified<-((length(which(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status=='TP'))+length(which(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status=='TN')))/length(FUN_OCPN_N_mu_00_s4505_4505_cor0_final_p100$Status))*100

Correctly_Classified
