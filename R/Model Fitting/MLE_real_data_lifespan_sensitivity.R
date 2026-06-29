library(optimx)
library(rootSolve)
library(tidyverse)
library(numDeriv)
library(dplyr)

# # vary the fixed values for gamma and omega
gamma=1/(0.5*365)
# gamma=1/(2*365)
omega=1/(5*365)
# omega=1/(15*365)
delta=1/21

# Load the village level data
data_village <- read.csv("village_level_data_w_FOI.csv", header=TRUE)
data <- data_village[,c(1:16)]

# Load the 200 randomly drawn lifespans
lifespan <- read.csv(file="Random_Lifespan.csv", header=TRUE)

#### Estimate Parameter for each set of lifespans
MLE_list <- list()
MLE_start <- list()

for(i in 1:200){
mu_goat = 1/(lifespan$Goat[i]*365)
mu_sheep = 1/(lifespan$Sheep[i]*365)
mu_bovine = 1/(lifespan$Cow[i]*365)
mu_camel = 1/(lifespan$Camel[i]*365)

data$b_goat <- data$size_goat * mu_goat
data$b_sheep <- data$size_sheep * mu_sheep
data$b_bovine <- data$size_bovine * mu_bovine
data$b_camel <- data$size_camel * mu_camel

# Ordinary differential equation model 3 with composite parameter A in rootsolve for each county

####Marsabit
initial_MAR_reduced = c(Sw=100, Sx=100, Sy=100, Sz=100,
                        Iw=100, Ix=100, Iy=100, Iz=100,
                        Rw=100, Rx=100, Ry=100, Rz=100,
                        Qw=100, Qx=100, Qy=100, Qz=100)
bruce_SP_model_MAR_reduced <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    muw=mu_goat
    mux=mu_sheep
    muy=mu_bovine
    muz=mu_camel
    delta=delta
    gamma=gamma
    omega=omega
    
    Nw = Sw + Iw + Rw + Qw
    Nx = Sx + Ix + Rx + Qx
    Ny = Sy + Iy + Ry + Qy
    Nz = Sz + Iz + Rz + Qz
    N = Nw + Nx + Ny + Nz
    POSw = (Iw + Rw)
    POSx = (Ix + Rx)
    POSy = (Iy + Ry)
    POSz = (Iz + Rz)
    
    dS = function(S,b,mu){ b - beta/N*S*(Iw+rho*(Rw+Qw)+
                                           Ix+rho*(Rx+Qx)+
                                           Iy+rho*(Ry+Qy)+
                                           Iz+rho*(Rz+Qz) ) - 
        A*(rho*(Qw+Qx+Qy+Qz+Rw+Rx+Ry+Rz)+Iw+Ix+Iy+Iz)*S - mu*S }
    dI = function(S,I,mu){ beta/N*S*( Iw+rho*(Rw+Qw)+
                                        Ix+rho*(Rx+Qx)+
                                        Iy+rho*(Ry+Qy)+
                                        Iz+rho*(Rz+Qz) ) + 
        A*(rho*(Qw+Qx+Qy+Qz+Rw+Rx+Ry+Rz)+Iw+Ix+Iy+Iz)*S - gamma*I - mu*I }
    dR = function(I,R,mu){ gamma*I - omega*R - mu*R }
    dQ = function(R,Q,mu){ omega*R - mu*Q }
    
    dSw = dS(Sw,bw,muw)
    dSx = dS(Sx,bx,mux)
    dSy = dS(Sy,by,muy)
    dSz = dS(Sz,bz,muz)
    
    dIw = dI(Sw,Iw,muw)
    dIx = dI(Sx,Ix,mux)
    dIy = dI(Sy,Iy,muy)
    dIz = dI(Sz,Iz,muz)
    
    dRw = dR(Iw,Rw,muw)
    dRx = dR(Ix,Rx,mux)
    dRy = dR(Iy,Ry,muy)
    dRz = dR(Iz,Rz,muz)
    
    dQw = dQ(Rw,Qw,muw)
    dQx = dQ(Rx,Qx,mux)
    dQy = dQ(Ry,Qy,muy)
    dQz = dQ(Rz,Qz,muz)
    
    list(c(dSw,dSx,dSy,dSz, 
           dIw,dIx,dIy,dIz,
           dRw,dRx,dRy,dRz,
           dQw,dQx,dQy,dQz), 
         Size=c(Nw,Nx,Ny,Nz), 
         Pos=c(POSw,POSx,POSy,POSz), 
         SeroPrev = c(POSw/Nw,POSx/Nx,POSy/Ny,POSz/Nz) )
  })
}
####Kajaido
initial_KAJ_reduced = c(Sw=100, Sx=100, Sy=100,
                        Iw=100, Ix=100, Iy=100,
                        Rw=100, Rx=100, Ry=100,
                        Qw=100, Qx=100, Qy=100)
bruce_SP_model_KAJ_reduced <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    muw=mu_goat
    mux=mu_sheep
    muy=mu_bovine
    delta=delta
    gamma=gamma
    omega=omega
    
    Nw = Sw + Iw + Rw + Qw
    Nx = Sx + Ix + Rx + Qx
    Ny = Sy + Iy + Ry + Qy
    N = Nw + Nx + Ny
    POSw = (Iw + Rw)
    POSx = (Ix + Rx)
    POSy = (Iy + Ry)
    
    dS = function(S,b,mu){ b - beta/N*S*(Iw+rho*(Rw+Qw)+
                                           Ix+rho*(Rx+Qx)+
                                           Iy+rho*(Ry+Qy)) - 
        A*(rho*(Qw+Qx+Qy+Rw+Rx+Ry)+Iw+Ix+Iy)*S - mu*S }
    dI = function(S,I,mu){ beta/N*S*( Iw+rho*(Rw+Qw)+
                                        Ix+rho*(Rx+Qx)+
                                        Iy+rho*(Ry+Qy)) + 
        A*(rho*(Qw+Qx+Qy+Rw+Rx+Ry)+Iw+Ix+Iy)*S - gamma*I - mu*I }
    dR = function(I,R,mu){ gamma*I - omega*R - mu*R }
    dQ = function(R,Q,mu){ omega*R - mu*Q }
    
    dSw = dS(Sw,bw,muw)
    dSx = dS(Sx,bx,mux)
    dSy = dS(Sy,by,muy)
    
    dIw = dI(Sw,Iw,muw)
    dIx = dI(Sx,Ix,mux)
    dIy = dI(Sy,Iy,muy)
    
    dRw = dR(Iw,Rw,muw)
    dRx = dR(Ix,Rx,mux)
    dRy = dR(Iy,Ry,muy)
    
    dQw = dQ(Rw,Qw,muw)
    dQx = dQ(Rx,Qx,mux)
    dQy = dQ(Ry,Qy,muy)
    
    list(c(dSw,dSx,dSy, 
           dIw,dIx,dIy,
           dRw,dRx,dRy,
           dQw,dQx,dQy), 
         Size=c(Nw,Nx,Ny), 
         Pos=c(POSw,POSx,POSy), 
         SeroPrev = c(POSw/Nw,POSx/Nx,POSy/Ny) )
  })
}

##############Function to calculate R0
R0_4spp <- function(size_goat,size_sheep,size_bovine,size_camel,N_t,A,beta,gamma,rho){
  (beta/N_t+A)*
    ( size_goat/mu_goat*(1+gamma*rho/mu_goat) / (1+gamma/mu_goat) +
        size_sheep/mu_sheep*(1+gamma*rho/mu_sheep) / (1+gamma/mu_sheep) +
        size_bovine/mu_bovine*(1+gamma*rho/mu_bovine) / (1+gamma/mu_bovine) +
        size_camel/mu_camel*(1+gamma*rho/mu_camel) / (1+gamma/mu_camel) )
}

#####################Function to calculate the negative log likelihood for real data
nloglik_4spp <- function(theta){
  A=theta[1]
  beta=theta[2]
  rho=theta[3]
  if(A<0|beta<0|rho<0|rho>10e-04){
    loglik=log(0)}else{
      
      try_R0 <- R0_4spp(size_goat=data$size_goat, 
                        size_sheep=data$size_sheep, 
                        size_bovine=data$size_bovine, 
                        size_camel=data$size_camel,
                        N_t=data$size_livestock,
                        A = theta[1]/1000, 
                        beta = theta[2], 
                        gamma = gamma, 
                        rho = theta[3]*1000)
      
      if (all(try_R0 > 1 & try_R0 < 2)){
        
        data_MAR <- subset(data, size_camel!=0)
        goatSeropos_MAR=vector()
        sheepSeropos_MAR=vector()
        bovineSeropos_MAR=vector()
        camelSeropos_MAR=vector()
        for (i in 1:nrow(data_MAR)){
          Pars <- c(theta*c(1/1000,1,1000),
                    bw=data_MAR$b_goat[i], bx=data_MAR$b_sheep[i],
                    by=data_MAR$b_bovine[i], bz=data_MAR$b_camel[i])
          SS <- stode(y=initial_MAR_reduced, time = 1000,
                      func = bruce_SP_model_MAR_reduced, parms = Pars, pos = TRUE, verbose = FALSE)
          goatSeropos_MAR[i] <- SS$SeroPrev[1]
          sheepSeropos_MAR[i] <- SS$SeroPrev[2]
          bovineSeropos_MAR[i] <- SS$SeroPrev[3]
          camelSeropos_MAR[i] <- SS$SeroPrev[4]
        }
        
        data_KAJ <- subset(data, size_camel==0)
        goatSeropos_KAJ=vector()
        sheepSeropos_KAJ=vector()
        bovineSeropos_KAJ=vector()
        for (i in 1:nrow(data_KAJ)){
          Pars <- c(theta*c(1/1000,1,1000),
                    bw=data_KAJ$b_goat[i], bx=data_KAJ$b_sheep[i],
                    by=data_KAJ$b_bovine[i])
          SS <- stode(y=initial_KAJ_reduced, time = 1000,
                      func = bruce_SP_model_KAJ_reduced  , parms = Pars, pos = TRUE, verbose = FALSE)
          goatSeropos_KAJ[i] <- SS$SeroPrev[1]
          sheepSeropos_KAJ[i] <- SS$SeroPrev[2]
          bovineSeropos_KAJ[i] <- SS$SeroPrev[3]
        }
        
        loglik_goat_MAR = dbinom(data_MAR$n_positive_goat, size=data_MAR$n_tested_goat, prob=goatSeropos_MAR, log=TRUE)
        loglik_sheep_MAR  = dbinom(data_MAR$n_positive_sheep, size=data_MAR$n_tested_sheep, prob=sheepSeropos_MAR, log=TRUE)
        loglik_bovine_MAR  = dbinom(data_MAR$n_positive_bovine, size=data_MAR$n_tested_bovine, prob=bovineSeropos_MAR, log=TRUE)
        loglik_camel_MAR  = dbinom(data_MAR$n_positive_camel, size=data_MAR$n_tested_camel, prob=camelSeropos_MAR, log=TRUE)
        
        loglik_MAR = sum(c(loglik_goat_MAR,loglik_sheep_MAR,loglik_bovine_MAR,loglik_camel_MAR))
        
        loglik_goat_KAJ = dbinom(data_KAJ$n_positive_goat, size=data_KAJ$n_tested_goat, prob=goatSeropos_KAJ, log=TRUE)
        loglik_sheep_KAJ = dbinom(data_KAJ$n_positive_sheep, size=data_KAJ$n_tested_sheep, prob=sheepSeropos_KAJ, log=TRUE)
        loglik_bovine_KAJ = dbinom(data_KAJ$n_positive_bovine, size=data_KAJ$n_tested_bovine, prob=bovineSeropos_KAJ, log=TRUE)
        
        loglik_KAJ = sum(c(loglik_goat_KAJ,loglik_sheep_KAJ,loglik_bovine_KAJ))
        
        loglik = loglik_MAR + loglik_KAJ 
      } else{
        loglik=log(0)
      }
    }
  return(-loglik)
}


########Multi-start MLE


library(foreach)
library(doParallel)
n_cores <- parallel::detectCores() - 1  # use all but one core
cl <- makeCluster(n_cores)
registerDoParallel(cl)

n_jobs <- 10  # number of MLE fits

mle <- foreach(i = 1:n_jobs, .packages =  c("optimx", "rootSolve"), .combine = rbind) %dopar% {
  omega=omega
  gamma=gamma
  repeat {
    start <- c(
      A     = runif(1, 0.00001 * 0.00001 * delta, 0.005 * 0.005 * delta) * 1000,
      beta  = runif(1, 0.0001, 0.005),
      rho   = runif(1, 0.3, 0.6) / 1000
    )
    if (!is.infinite(nloglik_4spp(start))) break
  }
  
  m <- try(optimx(start, nloglik_4spp, method = "Nelder-Mead",
                  hessian = FALSE,
                  control = list(maxit = 1000, trace = 0, kkt = FALSE)), silent = TRUE)
  
  names(start) <- c("A_start", "beta_start","rho_start")
  
  if (inherits(m, "try-error")) {
    return(c(start,
             A_hat = NA, 
             beta_hat = NA, 
             rho_hat = NA,
             nloglik = NA, 
             code = NA))
  } else {
    return(c(start,
             A_hat = m$A / 1000,
             beta_hat = m$beta,
             rho_hat = m$rho * 1000,
             nloglik = m$value,
             code = m$convcode))
  }
}
stopCluster(cl)
mle <- data.frame(mle)
rownames(mle) <- NULL
MLE_start[[i]] <- mle

mle_min <- mle[which.min(mle$nloglik), ]

MLE_list[[i]] <- mle_min

  print(i)
}


library(purrr)
real_data_results <- do.call(rbind, Map(cbind, index = seq_along(MLE_list), MLE_list))

# write.csv(real_data_results, file="Lifespan_sensitivity_estimation_gammahalf_omega5.csv")
# write.csv(real_data_results, file="Lifespan_sensitivity_estimation_gamma2_omega5.csv")
# write.csv(real_data_results, file="Lifespan_sensitivity_estimation_gammahalf_omega15.csv")
# write.csv(real_data_results, file="Lifespan_sensitivity_estimation_gamma2_omega15.csv")


