library(optimx)
library(rootSolve)
library(tidyverse)
library(numDeriv)
library(dplyr)

# Fixed parameter values
mu_goat=1/(2*365)
mu_sheep=1/(3*365)
mu_bovine=1/(6*365)
mu_camel=1/(16*365)
delta=1/21

# Ordinary differential equation model 3 in rootsolve for each county

####Marsabit
initial_MAR = c(Sw=100, Sx=100, Sy=100, Sz=100,
                Iw=10, Ix=10, Iy=10, Iz=10,
                Rw=10, Rx=10, Ry=10, Rz=10,
                Qw=10, Qx=10, Qy=10, Qz=10,
                B=10)
bruce_SP_model_MAR <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    muw=mu_goat
    mux=mu_sheep
    muy=mu_bovine
    muz=mu_camel
    delta=delta
    
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
                                           Iz+rho*(Rz+Qz) ) - xi*B*S - mu*S }
    dI = function(S,I,mu){ beta/N*S*( Iw+rho*(Rw+Qw)+
                                        Ix+rho*(Rx+Qx)+
                                        Iy+rho*(Ry+Qy)+
                                        Iz+rho*(Rz+Qz) ) + xi*B*S - gamma*I - mu*I }
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
    
    dB = alpha*(Iw+rho*(Rw+Qw)+
                  Ix+rho*(Rx+Qx)+
                  Iy+rho*(Ry+Qy)+
                  Iz+rho*(Rz+Qz) ) - delta*B
    
    list(c(dSw,dSx,dSy,dSz, 
           dIw,dIx,dIy,dIz,
           dRw,dRx,dRy,dRz,
           dQw,dQx,dQy,dQz,
           dB), 
         Size=c(Nw,Nx,Ny,Nz), 
         Pos=c(POSw,POSx,POSy,POSz), 
         SeroPrev = c(POSw/Nw,POSx/Nx,POSy/Ny,POSz/Nz) )
  })
}
####Kajaido
initial_KAJ = c(Sw=100, Sx=100, Sy=100,
                Iw=10, Ix=10, Iy=10,
                Rw=10, Rx=10, Ry=10,
                Qw=10, Qx=10, Qy=10,
                B=10)
bruce_SP_model_KAJ <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    muw=mu_goat
    mux=mu_sheep
    muy=mu_bovine
    delta=delta
    
    Nw = Sw + Iw + Rw + Qw
    Nx = Sx + Ix + Rx + Qx
    Ny = Sy + Iy + Ry + Qy
    N = Nw + Nx + Ny
    POSw = (Iw + Rw)
    POSx = (Ix + Rx)
    POSy = (Iy + Ry)
    
    dS = function(S,b,mu){ b - beta/N*S*(Iw+rho*(Rw+Qw)+
                                           Ix+rho*(Rx+Qx)+
                                           Iy+rho*(Ry+Qy)) - xi*B*S - mu*S }
    dI = function(S,I,mu){ beta/N*S*( Iw+rho*(Rw+Qw)+
                                        Ix+rho*(Rx+Qx)+
                                        Iy+rho*(Ry+Qy)) + xi*B*S - gamma*I - mu*I }
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
    
    dB = alpha*(Iw+rho*(Rw+Qw)+
                  Ix+rho*(Rx+Qx)+
                  Iy+rho*(Ry+Qy) ) - delta*B
    
    list(c(dSw,dSx,dSy, 
           dIw,dIx,dIy,
           dRw,dRx,dRy,
           dQw,dQx,dQy,
           dB), 
         Size=c(Nw,Nx,Ny), 
         Pos=c(POSw,POSx,POSy), 
         SeroPrev = c(POSw/Nw,POSx/Nx,POSy/Ny) )
  })
}

############function to calc R0

R0_4spp <- function(size_goat,size_sheep,size_bovine,size_camel,N_t,A,beta,gamma,rho){
  (beta/N_t+A)*
    ( size_goat/mu_goat*(1+gamma*rho/mu_goat) / (1+gamma/mu_goat) +
        size_sheep/mu_sheep*(1+gamma*rho/mu_sheep) / (1+gamma/mu_sheep) +
        size_bovine/mu_bovine*(1+gamma*rho/mu_bovine) / (1+gamma/mu_bovine) +
        size_camel/mu_camel*(1+gamma*rho/mu_camel) / (1+gamma/mu_camel) )
}

#### Load real data to use herd level information
data_village <- read.csv("village_level_data_w_FOI.csv", header=TRUE)
data <- data_village[,c(1:12)]
data$b_goat <- data$size_goat * mu_goat
data$b_sheep <- data$size_sheep * mu_sheep
data$b_bovine <- data$size_bovine * mu_bovine
data$b_camel <- data$size_camel * mu_camel

#####################Make simulated data

# number of datasets
set=50
#herd size
n=15
# vary the fixed values for gamma and omega
gamma=1/(0.5*365)
# gamma=1/(2*365)
# omega=1/(5*365)
omega=1/(15*365)

# generate simulate data according to criteria
data_list <- list()
sim_pars_list <- list()
for(j in 1:set){
  # select parameter combinations with 1<R0<2
  repeat {
    sim_pars <- c(
      alpha = runif(1, 0.00001, 0.005),
      beta  = runif(1, 0.0001, 0.005),
      gamma = gamma,
      omega = omega,
      rho   = runif(1, 0.2,0.8),
      xi    = runif(1, 0.00001, 0.005)
    )
    R0_herd_true <- R0_4spp(size_goat=data$size_goat, 
                    size_sheep=data$size_sheep, 
                    size_bovine=data$size_bovine, 
                    size_camel=data$size_camel,
                    N_t=data$size_livestock,
                    A=(sim_pars[1]*sim_pars[6])/delta,
                    beta=sim_pars[2],
                    gamma=sim_pars[3],
                    rho=sim_pars[5])
    
    if (all(R0_herd_true > 1 & R0_herd_true < 2)) break
  }                 
  # make number of seropositive in each herd with the selected parameters
  df_seropos <- data.frame(seropositive_goat = NA_real_,
                           seropositive_sheep = NA_real_,
                           seropositive_bovine  = NA_real_,
                           seropositive_camel  = NA_real_)
  for (i in 1:n) {
    if(data$size_camel[i] !=0){
      Pars <- c(sim_pars,
                bw = data$b_goat[i], bx = data$b_sheep[i],
                by = data$b_bovine[i], bz = data$b_camel[i])
      
      SS <- stode(y = initial_MAR, time = 1000, func = bruce_SP_model_MAR,
                  parms = Pars, pos = TRUE, verbose = FALSE)
      df_seropos[i,] <- SS$SeroPrev
    } else {
      Pars <- c(sim_pars,
                bw = data$b_goat[i], bx = data$b_sheep[i],
                by = data$b_bovine[i])
      
      SS <- stode(y = initial_KAJ, time = 1000, func = bruce_SP_model_KAJ,
                  parms = Pars, pos = TRUE, verbose = FALSE)
      
      df_seropos[i, ] <- c(SS$SeroPrev,0)
    }
  } 
  
  df_numPos <- data.frame(
    n_positive_goat = rbinom(n, data$n_tested_goat, df_seropos$seropositive_goat),
    n_positive_sheep = rbinom(n, data$n_tested_sheep, df_seropos$seropositive_sheep),
    n_positive_bovine = rbinom(n, data$n_tested_bovine, df_seropos$seropositive_bovine),
    n_positive_camel = rbinom(n, data$n_tested_camel, df_seropos$seropositive_camel)
  )
  df_numPos$n_positive_livestock = df_numPos$n_positive_goat+df_numPos$n_positive_sheep+
    df_numPos$n_positive_bovine+df_numPos$n_positive_camel
  
  data_list[[j]] <- cbind(data, R0_herd_true, df_seropos, df_numPos)
  sim_pars["A"] <- (sim_pars["alpha"]*sim_pars["xi"])/delta
  sim_pars_list[[j]] <- sim_pars
  print(j)
}




