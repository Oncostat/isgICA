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
library(class)
library(RcppHungarian)
library(psych)

#####################################################################
################# Path to the code directory ########################

setwd("C:/Users/Utilisateur/Downloads/SarahLaure/These/axe1_adap/2021_article_CMMM_code/99_code_article_2021/Z_BBP/")

#####################################################################
###################### Initialization ###############################

source("simulation_isgICA.r")
source("Rmodel_isgICA.r")

cmp_isgICA_model = cmpfun(isgICA_model)

N=600 # individuales
P=1200 # covariables
K = 10 # latent components
value_var_error = 2 # error variance 
value_sd_error = sqrt(value_var_error)
nb_seed = 1
vec_seed = rep(0,2)

Zth<-rbind(matrix(c(rep(1,N), # baseline component
                    rep(c(rep(1,round(N*0.3,0)),rep(0,round(N*0.6,0)),rep(1,round(N*0.1,0))  ),2),
                    rep(c(rep(0,round(N*0.2,0)),rep(1,round(N*0.3,0)),rep(0,round(N*0.5,0))),3),
                    rep(c(rep(0,round(N*0.4,0)),rep(1,round(N*0.2,0)),rep(0,round(N*0.4,0))),1),
                    rep(c(rep(0,round(N*0.6,0)),rep(1,round(N*0.1,0)),rep(0,round(N*0.3,0))),1),
                    rep(c(rep(0,round(N*0.7,0)),rep(1,round(N*0.3,0))),2)  ),nrow = K,ncol = N,byrow = T))

print(dim(Zth))
BZ = Zth 
table(c(BZ))

val = 3
mean_phi = seq(-K*val/K,K*val/K, length=K)
mean_W = rep(0,K) 
mean_W = mean_W[order(abs(mean_W))]

B<-matrix(0,nrow = P,ncol = K)
end<-1
set.seed(vec_seed[2]); (group = rbinom(n = K,size = 10,prob = 0.5)+1)
pinit = sapply(group, function(i) sample(1:P,size = i,replace = T)) 
offset = sapply(group, function(i) sample(1:150,size = i,replace = T)) 

for(i in (1:K)){
  group[1]
  pinit[[1]]
  offset[[1]]
  for(j in 1:group[i]){
    B[pinit[[i]][j]:min(nrow(B),pinit[[i]][j]+offset[[i]][j] ) ,i] = 1
  }
}

sim = simu_z_fixe(N = N,P = P,
                  Kmax = Kmax,
                  Zth=Zth,K=10, mean_phi = mean_phi,
                  mean_W = mean_W,Z_Phi = B,
                  error_var_prior =value_sd_error)

X=sim$X
Zth=(sim$Z>0.5)*1
Wth=sim$W
ZWth=sim$ZW
phith=sim$phi
w_var_Sim= sim$w_var
eta_var_Sim= sim$eta_var
phi_var_Sim= sim$phi_var
K_sim=sim$K_sim
Kmax = 100
maxit = 10000
type = 1# 1 == "diag" ; 2=="random"
typeZ = "diag"


#####################################################################
############### Representation of Zth and B ######################### 

meltedZth <- melt(Zth > 0.5)
meltedB <- melt(B > 0.5)
grid.arrange(
  ggplot(meltedZth, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +
    scale_fill_manual(values = c("white", "black")) +
    theme_bw() +
    theme(legend.position = "none")+ ggtitle("Real Z")+
    ylab("Latent Components") + xlab("Individuals")
  ,
  ggplot(meltedB, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +
    scale_fill_manual(values = c("white", "black")) +
    theme_bw() +
    theme(legend.position = "none")+ ggtitle("Real Phi")+
    ylab("Covariables") + xlab("Latent components")
)

save(sim,B,BZ,vec_seed,type, typeZ,value_var_error,value_sd_error, file = paste0( "simu_typeZ_",typeZ,"_var_error_",value_var_error,"_N_",N,"_P_",P,"_K_",K, "_Z_",table(c(BZ))[2] / sum(table(c(BZ))), "_Zphi_",table(c(B))[2] / sum(table(c(B))),"_seedZ_",vec_seed[1], "_seedPhi_",vec_seed[2] ,".rdata"))

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
               vec_seed = vec_seed,
               value_var_error=value_var_error,
               type=type,  
               typeZ=typeZ,
               sim = sim)

getTestFun <- function(initPar) {
  fBO=function(U){
    mu = qbeta(U,0.25,2.5)
    b_prior=1/mu-1
    initPar$b_prior = b_prior
    fit = cmp_isgICA_model(initPar$X,Kmax=initPar$Kmax,
                           maxit = initPar$maxit,seed=initPar$seed,
                           val_change = initPar$val_change, 
                           iniVarW = 1, iniVarPhi = 1,a=1,b=b_prior,
                           initPar$cc,initPar$dd,initPar$e,initPar$f)
    
    save(fit,initPar,
         file = paste0( 
           "typeZ_",initPar$typeZ,"_simu_sd_error_",initPar$value_var_error,
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
      "fit/typeZ_",initPar$typeZ,"_simu_sd_error_",initPar$value_var_error,"_seed_", initPar$vec_seed[1],
      "_N_",initPar$N,"_P_",initPar$P,"_ite",initPar$maxit, ".rdata"))
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
# fit = cmp_isgICA_model(initPar$X,Kmax=initPar$Kmax,
#                              maxit = initPar$maxit,seed=initPar$seed,
#                              val_change = initPar$val_change,
#                              iniVarW = 1, iniVarPhi = 1,a=1,b=b_prior,
#                              initPar$cc,initPar$dd,initPar$e,initPar$f)
# save(fit,initPar,
#      file = paste0( 
#        "typeZ_",initPar$typeZ,"_simu_sd_error_",initPar$value_var_error,
#        "_N_",initPar$N,"_P_",initPar$P,"_a_prior_",initPar$a_prior, "_b_prior_",
#        initPar$b_prior,
#        "_c_",initPar$cc, #"_g_",gg,"_h_",hh,
#        "_ite",initPar$maxit, ".rdata"))

#####################################################################
###################### Results for one epoch#########################

sorting <- function(TH,new){
  # TH = phith
  # new = pp
  colnames(TH) = paste0("K",1:ncol(TH))
  colnames(new) = paste0("K",1:ncol(new))
  cor <- psych::corr.test( TH, new,method="pearson", adjust="none")
  cost = abs(cor$r) 
  soln <- HungarianSolver(1-cost)
  return(list(new[soln$pairs[,2],],
              soln$pairs[,2]))
}

load("typeZ_diag_simu_sd_error_2_N_600_P_1200_a_prior_1_b_prior_10000_c_1_ite10000.rdata")

zz = fit$Z
ligne = round(apply(zz>0.5,1,sum),1)
zz1 = (zz[ligne>=1,] > 0.5)*1
zw = fit$W * (fit$Z >0.5)
zw = zw[ligne>=1,]
pp = fit$Phi
pp = pp[,which(ligne>=1)]

tryCatch({
  res_sort = sorting(phith,pp)
  arr = res_sort[[2]]
  zz2 = zz1[arr,]
}, error=function(e){print(ll[ii]) ; cat("ERROR :",conditionMessage(e), "\n") })

pp = pp[,arr]
colnames(phith) = paste0("K",1:ncol(phith))
colnames(pp) = paste0("K",1:ncol(pp))
cor <- psych::corr.test(phith, pp,method="pearson", adjust="none")

ttable = table(zz2,Zth)
Se = ttable[2,2] / sum(ttable[,2])
Sp = ttable[1,1] / sum(ttable[,1])
acc = (ttable[1,1]+ttable[2,2]) / sum(ttable)

recap = c(value_var_error, 
          initPar$a_prior,
          initPar$b_prior,
          initPar$cc,
          initPar$dd,
          initPar$e,
          initPar$f,
          fit$lastElbo,
          fit$K[length(fit$K)],
          length(fit$K),
          sum(round(apply(fit$Z>0.5,1,sum),1)[round(apply(fit$Z>0.5,1,sum),1)>0] == ncol(fit$Z)) != fit$K[length(fit$K)],
          sum(round(apply(fit$Z>0.5,1,sum),1)[round(apply(fit$Z>0.5,1,sum),1)>0] == ncol(fit$Z)),
          fit$mse[length(fit$K)],
          Se,Sp,acc,
          mean(diag(abs(cor$r))))

names(recap) = c("var_error","a_prior","b_prior","c","d","e","f", "elbo", "K", "ite","sparse", "sparse2", "mse", "Se", "Sp","acc", "corPhi")
recap





