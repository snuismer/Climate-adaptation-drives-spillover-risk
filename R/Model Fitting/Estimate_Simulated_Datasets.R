# Must run Make_Simulated_Datasets.R first

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
    omega=omega
    gamma=gamma

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
    omega=omega
    gamma=gamma
    
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


#####################Function to calculate the negative log likelihood
nloglik_4spp <- function(theta){
  A=theta[1]
  beta=theta[2]
  rho=theta[3]
  if(A<0|beta<0|rho<0|rho>1e-03){
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

# Estimate the parameters in each dataset
MLE_list <- list()
# Calculate R0 from the MLE
R0_hat_list <- list()
# Keep random starting values
MLE_start <- list()

############### Run Optimization for all datasets
library(foreach)
library(doParallel)

for(j in 1:set){
  data <- data_list[[j]]
  true_params <- c(sim_pars_list[[j]]["A"],
                   sim_pars_list[[j]]["beta"],
                   sim_pars_list[[j]]["gamma"],
                   sim_pars_list[[j]]["omega"],
                   sim_pars_list[[j]]["rho"])
  names(true_params) <- c("A_true", "beta_true", "gamma_true", "omega_true", "rho_true")
  n_cores <- parallel::detectCores() - 1  # use all but one core
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  n_jobs <- 10  # number of starting conditions
  
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
      return(c(true_params, start,
               dataset=j,
               A_hat = NA, beta_hat = NA, 
               rho_hat = NA,
               nloglik = NA, code = NA))
    } else {
      return(c(true_params, start,
               dataset=j,
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
  MLE_start[[j]] <- mle
  # keep starting value with best likelihood
  mle_min <- mle[which.min(mle$nloglik), ]
  
  MLE_list[[j]] <- mle_min
  R0_hat_list[[j]] <- R0_4spp(size_goat=data_list[[j]]$size_goat,
                              size_sheep=data_list[[j]]$size_sheep,
                              size_bovine=data_list[[j]]$size_bovine,
                              size_camel=data_list[[j]]$size_camel,
                              N_t=data_list[[j]]$size_livestock,
                              A = mle_min[["A_hat"]], 
                              beta = mle_min[["beta_hat"]], 
                              gamma = gamma, 
                              rho = mle_min[["rho_hat"]])
  
  print(j)
}



results <- as.data.frame(do.call(rbind, MLE_list))

data_list <- Map(function(df, val) {
  df$R0_herd_hat <- val
  return(df)
}, data_list, R0_hat_list)


mean_df <- data.frame(
  mean_R0_true = sapply(data_list, function(df) mean(df$R0_herd_true, na.rm = TRUE)),
  mean_R0_hat = sapply(data_list, function(df) mean(df$R0_herd_hat, na.rm = TRUE))
)

MLE_results <- cbind(results, mean_df)

View(MLE_results)

# write.csv(MLE_results, file="Sim_study_results_gammahalf_omega5.csv")
# write.csv(MLE_results, file="Sim_study_results_gamma2_omega5.csv")
# write.csv(MLE_results, file="Sim_study_results_gammahalf_omega15.csv")
# write.csv(MLE_results, file="Sim_study_results_gamma2_omega15.csv")
