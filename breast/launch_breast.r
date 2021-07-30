# write by Sarah-Laure Rincourt
# Corresponding author : damien.drubay@gustaveroussy.fr 
# licence GNU3 

#####################################################################
########################### library #################################

library(zeallot) # to use :  %<-%
library(glue)
library(GPfit)
library(magrittr) # for %>%
library(ggplot2)
library(gridExtra)
library(reshape2) #melt 
library(ParBayesianOptimization)
library(doParallel)
library(compiler)

#####################################################################
############################ define setwd ###########################

setwd("C:/Users/Utilisateur/Downloads/SarahLaure/These/axe1_adap/2021_article_CMMM_code/99_code_article_2021/Z_BBP/")

#####################################################################
###################### Initialization ###############################

# load the model
source("Rmodel_isgICA.r")
cmp_isgICA_model = cmpfun(isgICA_model)

# load the data
Breast = read.table("Breast_filtration.txt")
bb = Breast 
dim(bb)
seed = 1

X = t(bb)
P = dim(X)[1] # covariables
N = dim(X)[2] # individuals
Kmax = 100L # maximum of latent components
vec_seed = rep(0,2)
maxit=10000 # number of iterations 



#####################################################################
####################  Variational inference #########################

# define the list of required elements for the inference
initPar = list(X=X, N=N, P=P,Kmax=Kmax,
               maxit = maxit ,
               seed=0, val_change = 1e-6,
               a_prior=1,
               cc=1,
               dd=1,
               e=1e-6,
               f=1e-6,
               vec_seed = vec_seed)

getTestFun <- function(initPar) {
  fBO=function(U){
    mu = qbeta(U,0.25,2.5)
    b_prior=1/mu-1
    initPar$b_prior = b_prior
    fit = cmp_BPFA_VB_profil_svd(initPar$X,Kmax=initPar$Kmax,
                                 maxit = initPar$maxit,seed=initPar$seed,
                                 val_change = initPar$val_change, 
                                 iniVarW = 1, iniVarPhi = 1,a=1,b=b_prior,
                                 initPar$cc,initPar$dd,initPar$e,initPar$f)
    save(fit,initPar,
         file = paste0("breast_",
           "_N_",initPar$N,"_P_",initPar$P,"_a_prior_",initPar$a_prior, "_b_prior_",
           initPar$b_prior,
           "_c_",initPar$cc, #"_g_",gg,"_h_",hh,
           "_ite",initPar$maxit, ".rdata"))
    return(list(Score=fit$lastElbo))
  }
  fBO
}


Do_Opt=function(initPar) {
  Test_Fun <- getTestFun(initPar)
  beta = c(1,10,100,1000,10000,100000) 
  mu = 1 / (beta +1)
  U = pbeta(mu, 0.25,2.5)
  bornes = pbeta(1 / (c(1e6,1e-6)+1), 0.25,2.5)
  
  tryCatch( 
    fits <- bayesOpt(FUN = Test_Fun, bounds = list(U=bornes), #0.999999)), 
                     # initPoints = 4,
                     initGrid = data.frame(U),
                     iters.n = 24,iters.k=1, 
                     gsPoints = 1000,parallel=F,errorHandling="continue"),
    error = function(e){
      write(as.character(paste0("vec_seed=", initPar$vec_seed ,
                                " value_var_error=",
                                initPar$value_var_error,
                                " typeZ=", initPar$typeZ,
                                " Erreur: stop bayesOpt")),
            "test2.txt", append=TRUE)
      write(as.character(e), "test2.txt", append=TRUE)
    })
  if(exists("fits")){
    save(fits,initPar, file=paste0(
      "fit_breast","_N_",initPar$N,"_P_",initPar$P,"_ite",initPar$maxit, ".rdata"))
    return(fits)
  }else{
    return(NA)
  }
}

#  Complete process of the gaussian optimization
Do_Opt(initPar = initPar)


# ## Exemple based on one epoch with beta = 1e4
# 
# b_prior=1e4
# initPar$b_prior = b_prior
# fit = cmp_BPFA_VB_profil_svd(initPar$X,Kmax=initPar$Kmax,
#                              maxit = initPar$maxit,seed=initPar$seed,
#                              val_change = initPar$val_change,
#                              iniVarW = 1, iniVarPhi = 1,a=1,b=b_prior,
#                              initPar$cc,initPar$dd,initPar$e,initPar$f)
# save(fit,initPar,
#      file = paste0("breast_",
#                    "_N_",initPar$N,"_P_",initPar$P,"_a_prior_",initPar$a_prior, "_b_prior_",
#                    initPar$b_prior,
#                    "_c_",initPar$cc, #"_g_",gg,"_h_",hh,
#                    "_ite",initPar$maxit, ".rdata"))



##################################################################

