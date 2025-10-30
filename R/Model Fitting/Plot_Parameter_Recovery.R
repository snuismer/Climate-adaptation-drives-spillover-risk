library(ggplot2)
library(latex2exp)
library(ggpubr)

# Load parameter recover results and remove rows that did not converge during optimization
MLE_set_1 <- read.csv("Sim_study_results_gammahalf_omega5.csv", header=TRUE)
MLE_set_1 <- subset(MLE_set_1, code==0)
MLE_set_2 <- read.csv("Sim_study_results_gamma2_omega5.csv", header=TRUE)
MLE_set_2 <- subset(MLE_set_2, code==0)
MLE_set_3 <- read.csv("Sim_study_results_gammahalf_omega15.csv", header=TRUE)
MLE_set_3 <- subset(MLE_set_3, code==0)
MLE_set_4 <- read.csv("Sim_study_results_gamma2_omega15.csv", header=TRUE)
MLE_set_4 <- subset(MLE_set_4, code==0)


###########Plot results

#####SET 1
# Make linear prediction for plots
Amod <- summary(lm(A_hat ~ A_true , data = MLE_set_1))
betamod <- summary(lm(beta_hat ~ beta_true , data = MLE_set_1))
rhomod <- summary(lm(rho_hat ~ rho_true , data = MLE_set_1))

## Parameter recovery plots
p1 <- ggplot(MLE_set_1, aes(x=A_true, y=A_hat)) +
  geom_point( color="steelblue4",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$A$")) + ylab(TeX("$\\hat{A}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = Amod$coefficients[1,1], slope = Amod$coefficients[2,1], 
              linetype=1, color="steelblue4", size=1)+
  xlim(c(0,4e-07))+ylim(c(0,4e-07))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p2 <- ggplot(MLE_set_1, aes(x=beta_true, y=beta_hat)) +
  geom_point(color="steelblue3", shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\beta$")) + ylab(TeX("$\\hat{\\beta}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = betamod$coefficients[1,1], slope = betamod$coefficients[2,1], 
              linetype=1, color="steelblue3", size=1)+
  xlim(c(0,0.005))+ylim(c(0,0.005))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p3 <- ggplot(MLE_set_1, aes(x=rho_true, y=rho_hat)) +
  geom_point( color="steelblue2",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = rhomod$coefficients[1,1], slope = rhomod$coefficients[2,1], 
              linetype=1, color="steelblue2", size=1)+
  xlim(c(0,1))+ylim(c(0,1))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

set_1_plot <- ggarrange(p1,p2,p3, nrow=1, ncol = 3, heights = c(1,1,1))



#####SET 2
# Make linear prediction for plots
Amod <- summary(lm(A_hat ~ A_true , data = MLE_set_2))
betamod <- summary(lm(beta_hat ~ beta_true , data = MLE_set_2))
rhomod <- summary(lm(rho_hat ~ rho_true , data = MLE_set_2))

## Parameter recovery plots
p4 <- ggplot(MLE_set_2, aes(x=A_true, y=A_hat)) +
  geom_point( color="cadetblue4",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$A$")) + ylab(TeX("$\\hat{A}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = Amod$coefficients[1,1], slope = Amod$coefficients[2,1], 
              linetype=1, color="cadetblue4", size=1)+
  xlim(c(0,4e-07))+ylim(c(0,4e-07))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p5 <- ggplot(MLE_set_2, aes(x=beta_true, y=beta_hat)) +
  geom_point(color="cadetblue3", shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\beta$")) + ylab(TeX("$\\hat{\\beta}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = betamod$coefficients[1,1], slope = betamod$coefficients[2,1], 
              linetype=1, color="cadetblue3", size=1)+
  xlim(c(0,0.005))+ylim(c(0,0.005))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p6 <- ggplot(MLE_set_2, aes(x=rho_true, y=rho_hat)) +
  geom_point( color="cadetblue2",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = rhomod$coefficients[1,1], slope = rhomod$coefficients[2,1], 
              linetype=1, color="cadetblue2", size=1)+
  xlim(c(0,1))+ylim(c(0,1))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

set_2_plot <- ggarrange(p4,p5,p6, nrow=1, ncol = 3, heights = c(1,1,1))



#####SET 3
# Make linear prediction for plots
Amod <- summary(lm(A_hat ~ A_true , data = MLE_set_3))
betamod <- summary(lm(beta_hat ~ beta_true , data = MLE_set_3))
rhomod <- summary(lm(rho_hat ~ rho_true , data = MLE_set_3))

## Parameter recovery plots
p7 <- ggplot(MLE_set_3, aes(x=A_true, y=A_hat)) +
  geom_point( color="khaki4",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$A$")) + ylab(TeX("$\\hat{A}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = Amod$coefficients[1,1], slope = Amod$coefficients[2,1], 
              linetype=1, color="khaki4", size=1)+
  xlim(c(0,4e-07))+ylim(c(0,4e-07))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p8 <- ggplot(MLE_set_3, aes(x=beta_true, y=beta_hat)) +
  geom_point(color="khaki3", shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\beta$")) + ylab(TeX("$\\hat{\\beta}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = betamod$coefficients[1,1], slope = betamod$coefficients[2,1], 
              linetype=1, color="khaki3", size=1)+
  xlim(c(0,0.005))+ylim(c(0,0.005))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p9 <- ggplot(MLE_set_3, aes(x=rho_true, y=rho_hat)) +
  geom_point( color="khaki2",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = rhomod$coefficients[1,1], slope = rhomod$coefficients[2,1], 
              linetype=1, color="khaki2", size=1)+
  xlim(c(0,1))+ylim(c(0,1))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

set_3_plot <- ggarrange(p7,p8,p9, nrow=1, ncol = 3, heights = c(1,1,1))


#####SET 4
# Make linear prediction for plots
Amod <- summary(lm(A_hat ~ A_true , data = MLE_set_4))
betamod <- summary(lm(beta_hat ~ beta_true , data = MLE_set_4))
rhomod <- summary(lm(rho_hat ~ rho_true , data = MLE_set_4))

## Parameter recovery plots
p10 <- ggplot(MLE_set_4, aes(x=A_true, y=A_hat)) +
  geom_point( color="palevioletred4",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$A$")) + ylab(TeX("$\\hat{A}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = Amod$coefficients[1,1], slope = Amod$coefficients[2,1], 
              linetype=1, color="palevioletred4", size=1)+
  xlim(c(0,4e-07))+ylim(c(0,4e-07))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

p11 <- ggplot(MLE_set_4, aes(x=beta_true, y=beta_hat)) +
  geom_point(color="palevioletred3", shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\beta$")) + ylab(TeX("$\\hat{\\beta}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = betamod$coefficients[1,1], slope = betamod$coefficients[2,1], 
              linetype=1, color="palevioletred3", size=1)+
  xlim(c(0,0.005))+ylim(c(0,0.005))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )
p12 <- ggplot(MLE_set_4, aes(x=rho_true, y=rho_hat)) +
  geom_point( color="palevioletred2",shape=20, size=3, alpha = 0.7) +
  xlab(TeX("$\\rho$")) + ylab(TeX("$\\hat{\\rho}$"))+
  geom_abline(intercept = 0, slope = 1, linetype=2, color="black")+
  geom_abline(intercept = rhomod$coefficients[1,1], slope = rhomod$coefficients[2,1], 
              linetype=1, color="palevioletred2", size=1)+
  xlim(c(0,1))+ylim(c(0,1))+
  theme_bw()+theme(text = element_text(size = 12),
                   axis.title =element_text(size = 16) )

set_4_plot <- ggarrange(p10,p11,p12, nrow=1, ncol = 3, heights = c(1,1,1))



#########Make final combined plot
combined_plot <- ggarrange(set_1_plot, set_2_plot, set_3_plot, set_4_plot,
                           nrow = 4, ncol = 1,
                           labels = c("a", "b", "c", "d") ,
                           font.label = list(size = 25, face = "bold"))
combined_plot
