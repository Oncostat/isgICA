# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

#####################################################################
########################### library #################################

library(mvtnorm)

#####################################################################
########################### generate data############################


simu_z_fixe <- function(P=50, N=100, Kmax=50,
                        w_var_prior=1,
                        phi_var_prior=1, error_var_prior=0.1,
                        rho=NULL, seed=1,Zth=NULL, K=10,
                        mean_phi=NULL,mean_W=NULL, Z_Phi=NULL){
  
  set.seed(seed)
  nbSimFactors = K
  if(is.null(Zth)){ 
    Zth = 1*(matrix(rnorm(K*N), nrow=K)>0)
  }
  if( K!= dim(Zth)[1]){
    stop("Dimension error between K and the row of Zth")
  }
  
  if(is.null(mean_phi)){ mean_phi <-  rep(0, nbSimFactors) }
  if(is.null(mean_W)){ mean_W <-  rep(0, nbSimFactors) }
  w_var <- rep(w_var_prior,nbSimFactors) # weight variance
  eta_var <- error_var_prior #error variance
  phi_var <- rep(phi_var_prior,nbSimFactors)
  
  Wth <- W <- t(rmvnorm(n=N, mean=mean_W, sigma = diag(w_var)))
  Eth <- E <- t(matrix(rnorm(N*P,0,eta_var),nrow=N,ncol=P))
  Phith <- phi<-matrix(rnorm(K*P,mean_phi,phi_var),nrow=P,ncol=K,byrow = T)
  
  if(is.null(Z_Phi)){
    B<-matrix(0,nrow = P,ncol = K)
    offset=0.05*P
    length=round(P/(K-1)+offset,1)
    end<-1
    for(i in (1:K)){
      start<-end-(i!=1)*offset
      end<-start+length-1
      B[start:min(end,nrow(B)),i]<-phi[start:min(end,nrow(B)),i]
    }
    end<-1
    for(i in (K:1)){
      start<-end-(i!=K)*offset
      end<-start+length-1
      B[start:min(end,nrow(B)),i]<-phi[start:min(end,nrow(B)),i]
    }
  }else{
    B = Z_Phi
    B = B*phi
  }
  Phith <- phi<- B
  
  X <- Xth <- Phith%*%(Wth*Zth)+Eth
  
  return(list("X"=X, 
              "Z"=Zth,
              "W"=Wth,
              "ZW"=Zth*W,
              "phi"=Phith,
              "w_var"= w_var_prior,
              "eta_var"= error_var_prior,
              "phi_var"= phi_var_prior,
              "K_sim"=K))
}

