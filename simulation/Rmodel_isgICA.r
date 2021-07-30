# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

#####################################################################
########################### library #################################

library(fBasics) # tr()
library(MASS) 
library(R.matlab)
library(invgamma) # digamma()

#####################################################################
########################### Code ####################################

# The following arguments are required to customize your model:
#   
# - the data in the dimension (P,N): `X` with P the number of genes and N the number of observations
# - the allowed number of latents factors: `Kmax`
# - the number of iterations of the variational inference algorithm: `maxit`
# - the seed of the algorithm: `seed`
# - the variance hyperparameter of W matrix: `inivarW`
# - the variance hyperparameter of $\Phi$ matrix: `iniVarPhi`
# - the variance hyperparameter of $\E$ matrix: `tau_error`
# - the variance hyperparameter of the gamma distribution: for $\Phi$ (`cc` and `dd`) ; for W (`e0` and `f0`)
# - the variance hyperparameter of the beta distribution: `a` and `b`

isgICA_model <- function(X,Kmax,maxit, seed, val_change=5e-6, iniVarW, iniVarPhi,a=1,b=Kmax,
                         cc=10^-1,dd=10^-1,e0=10^-1,f0=10^-1,tau_error=1){ 
  EWTW <- c()
  J <- c()
  W_sigma <- list()
  Phi_sigma <- list()
  set.seed(seed)
  eps=2e-8
  P = dim(X)[1]
  N = dim(X)[2]
  Kmax=min(Kmax,N,P)
  
  ######## hyperparameter for A matrix (covariance)  ########
  c0 =cc*matrix(1,nrow = P, ncol = Kmax) 
  d0 = dd*matrix(1,nrow = P, ncol = Kmax) 
  c = c0; d = d0;
  lgc0 = lgamma(c0) 
  
  e=e0;f=f0;
  
  # standardize the data
  X = t(scale(t(X)))
  tX = t(X)
  
  ######## Initialize A and S   ######## 
  #  X = U D V'
  
  SVD <- svd(X)
  Phi = SVD$u %*% diag(SVD$d)
  Phi = Phi[,1:Kmax]
  W = t(SVD$v)
  W = W[1:Kmax,]
  
  ##### Init random
  
  #Phi = matrix(rnorm(Kmax*P, sd=2), nrow = P, ncol = Kmax)
  #W = matrix(rnorm(N*Kmax, sd=2), nrow = Kmax, ncol = N)
  
  Z=matrix(T,nrow = Kmax, ncol = N) 
  ZW=Z*W;  
  sigmav = matrix(1, nrow = P, ncol = Kmax)
  
  tZW = t(ZW)
  
  V_error = matrix(tau_error, nrow = P, ncol = 1)
  vec_tr=c()
  
  ######## initialize pi  ######## 
  
  KmaxZ = Kmax-1 # mise en place du profil moyen
  
  pia0=c(1,a/KmaxZ* matrix(1, nrow = 1,ncol = KmaxZ) )
  pib0=c(1, b*(KmaxZ-1)/KmaxZ*matrix(1, nrow = 1,ncol = KmaxZ) )
  pia=pia0;pib=pib0;
  pia0_profil = pia0[-1]
  pib0_profil = pib0[-1]
  lgpia0=lgamma(pia0_profil)
  lgpib0=lgamma(pib0_profil)
  lgpiab0=lgamma(pia0_profil+pib0_profil)
  pai= 0.1*matrix(1, nrow = 1,ncol = Kmax) 
  
  Phi_sigma=matrix(1e-6, nrow = P,ncol = Kmax) #convariance of Phi_jk
  EPhiPhi = Phi*Phi+Phi_sigma;
  PhiV_error = Phi*matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax)
  PhiTPhi=t(PhiV_error)%*%Phi
  EPhiTPhi=t(Phi)%*%Phi+diag(colSums(Phi_sigma))
  sigmaw=1
  DWvar=matrix(0,Kmax,N)
  
  result <- list(
    W=list(),
    Z=list(),
    ZW=list(),
    Phi=list(),
    pai=list(),
    V_error=list(),
    mse=numeric(maxit), #c(),
    K=numeric(maxit), # c(),
    elbo=numeric(maxit), # c(),
    sigmav = list(),
    sigmaw = list(),
    hyperparameter = list(),
    lastElbo = c()
  )
  iter = 0; change = 1; 
  
  ######## While ########
  while (iter<maxit){
    start_time <- Sys.time()
    iter = iter + 1;
    
    ######## update Z  ########
    
    Zcnt = rowSums(Z) 
    exp_pi = cbind(digamma(pia0 + Zcnt) - digamma(pia0 + pib0 + N),
                   digamma(pib0 + N - Zcnt) - digamma(pia0 + pib0 + N) )
    
    PhiV_error = Phi*matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax)
    PhiTPhi=t(PhiV_error)%*%Phi
    EPhiPhiV_error = colSums(EPhiPhi * matrix(rep(V_error,Kmax), nrow=P, ncol = Kmax))
    nonZeroPos=which(Zcnt> eps)
    Z[which(Zcnt<=eps),]=0
    nonZeroPos=nonZeroPos[nonZeroPos!=1] 
    if(length(nonZeroPos)!=0){
      if(length(nonZeroPos)>1){
        obj=tX%*%PhiV_error[,nonZeroPos] - tZW%*%PhiTPhi[,nonZeroPos]  + t(diag(PhiTPhi)[nonZeroPos]*ZW[nonZeroPos,])
      }else{
        #insure the matrix object format
        obj=tX%*%matrix(PhiV_error[,nonZeroPos], ncol=1) - tZW%*%matrix(PhiTPhi[,nonZeroPos], ncol=1)  + matrix(t(diag(PhiTPhi)[nonZeroPos]*ZW[nonZeroPos,]), ncol=1)
      }
    }
    t1 = exp_pi[nonZeroPos,2]-(exp_pi[nonZeroPos,1] - 0.5*(( (W[nonZeroPos,]^2) + DWvar[nonZeroPos,])*EPhiPhiV_error[nonZeroPos]) + t(obj)*W[nonZeroPos,])
    Z[nonZeroPos,]=1/(1+exp(t1))  
    
    ######## Update for W    ######## 
    logdetWsig=0;
    t2=PhiTPhi
    diag(t2)=EPhiPhiV_error+(sigmaw+1);
    
    midval=tryCatch(solve(t2), 
                    error = function(e){
                      return(ginv(t2))
                    }) 
    
    t1 = t(PhiV_error)%*%X 
    diagMidval=diag(midval)
    DWvar=diagMidval*matrix(1,nrow=Kmax,ncol=N)
    W = midval%*%t1
    WW = W*W
    EWTW=colSums(WW) + sum(diagMidval) 
    logdetWsig=0.5*N*log(2*pi*(det(midval) +eps)) 
    
    ZW=Z*W
    tZW = t(ZW)
    ZWZW = ZW%*%tZW
    
    ######## Update  Phi_jk    ######## 
    
    magWZ = rowSums(Z*(WW + DWvar))
    Phi_sigma=1/(V_error%*%magWZ+sigmav)
    for(k in 1:Kmax){
      Phi[,k]=0
      Xmid = X%*%ZW[k,] - Phi%*%ZWZW[,k]
      Phi[,k]=Phi_sigma[,k]*V_error*Xmid 
    }
    EPhiTPhi = t(Phi)%*%Phi + diag(colSums(Phi_sigma))
    EPhiPhi = Phi*Phi + Phi_sigma
    
    ######## Update  pi_k    ######## 
    
    pia=rowSums(Z)+pia0
    pib=N-rowSums(Z)+pib0
    pai=pia/(pia+pib)
    
    ######## Update update sigma Phi    ######## 
    
    c=c0+0.5;
    d=d0+0.5*EPhiPhi;
    sigmav=c/d;  
    
    ######## Update sigma W    ######## 
    
    e=e0 + Kmax*N/2; #sum(apply(Z>=0.5,1, sum)>=1)*N/2; #
    f=sum(EWTW)/2 + f0;
    sigmaw = e/f;
    Zb = Z > 0.5
    
    ######## Save samples.   ########
    
    res=X-Phi%*%ZW
    result$mse[iter] = mean(sqrt(colSums(res*res)))
    result$K[iter] = sum(rowSums(Zb) >0) #sum(Zcnt > 0)
    
    print(paste("Elapsed time: ", round(Sys.time()-start_time,4),"; Iteration", iter))
    print(paste("K", result$K[iter],"; mse", round(result$mse[iter],4))) 
    
    if((iter%%100) ==0){ 
      result$W = W
      result$Phi = Phi
      result$Z = Z
      result$ZW = ZW
      result$pai = pai
      result$V_error = V_error
      result$sigmav = sigmav
      result$sigmaw = sigmaw
      result$hyperparameter = list(pia,pib,c,d,e,f)
      save(result,X,a,b,cc,dd,e0,f0, file = paste0("", "N_",N,"_P_",P, "_Kmax_",Kmax,"_a_prior_",a, "_b_prior_",b, "_ongoing", ".rdata"))
    }
  }
  result$W = W
  result$Phi = Phi
  result$Z = Z
  result$ZW = ZW
  result$pai = pai
  result$V_error = V_error
  result$sigmav = sigmav
  result$sigmaw = sigmaw
  result$hyperparameter = list(pia,pib,c,d,e,f) 
  ######## Lower Bound    ########
  J9 = 0
  tmp2 = digamma(c) - log(d)
  tmp3 = digamma(e) - log(f)      
  
  pia_profil = pia[-1]
  pib_profil = pib[-1]
  pai_profil = pai[-1]
  
  # log likelihood 
  predX = Phi%*%ZW
  midval = apply(X*(X-2*predX),1,sum)
  tmp5 = ZWZW-diag(diag(ZWZW))+diag(magWZ)
  vec_tr = c()
  for(j in 1:P){
    tmp4 = Phi[j,]%*%t(Phi[j,]) # attention au probleme de R sur les dimensions
    tmp4 = tmp4-diag(diag(tmp4))+diag(EPhiPhi[j,])
    trtmp = tr(tmp5%*%tmp4)
    vec_tr[j] = 0.5*(trtmp +midval[j])
  }
  j1bis = -P*N/2*log(2*pi)+sum(N*P/2*log(V_error)-V_error*vec_tr) 
  
  # KLD of gamma distribution
  j3 = - sum( c*log(d) - c0*log(d0) - lgamma(c) +lgc0 + (c-c0)*(tmp2) - c*(1-d0/d) )
  j11 = -(e*log(f)-e0*log(f0)-lgamma(e)+lgamma(e0)+(e-e0)*tmp3-e*(1-f0/f))
  
  # KL de Pi_k
  piab=pia_profil+pib_profil
  dgpia=digamma(pia_profil)
  dgpib=digamma(pib_profil)
  dgpiab=digamma(piab)
  j5 = sum(lgpiab0-lgamma(piab)+lgamma(pia_profil)-lgpia0+lgamma(pib_profil)-lgpib0-(pia_profil-pia0_profil)*(dgpia-dgpiab)-(pib_profil-pib0_profil)*(dgpib-dgpiab))
  j4 = (dgpia-dgpiab)%*%rowSums(Z[-1,])+(dgpib-dgpiab)%*%(N-rowSums(Z[-1,]));
  j2 =-0.5*P*Kmax*log(2*pi) + 0.5*sum(tmp2) - 0.5*sum((c/d)*(d-d0));
  j6 = (N*Kmax/2)*tmp3-0.5*sigmaw*sum(EWTW)
  j8 = 0.5*sum(log(2*pi*abs(Phi_sigma)))
  tmpZ = Z
  tmpZ[which(Z<eps)]=eps
  j9 = sum(tmpZ*log(tmpZ))
  j10 = logdetWsig
  result$lastElbo = j2+j3+j4+j5+j6+j8+j9+j10+j11 + j1bis # j7
  return(result)
}
