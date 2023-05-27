#Learning about EDM
#Ruby Krasnow
#5/25/23

library(rEDM)
library(Kendall)



# # time series of red noise and logistic map --------------------------------------------------------------------
# dat <- read.csv('ESM2_Data_noise.csv',header=T)
# 
# plot(dat$R)
# Red <- ((dat[,"R"]-mean(dat[,"R"]))/sd(dat[,"R"]))
# plot(Red)
# 
# 
# plot(dat$L)
# Logi <- ((dat[,"L"]-mean(dat[,"L"]))/sd(dat[,"L"]))
# plot(Logi)
# 
# ## Simplex projection for red noise and logistic map
# sim_r <- simplex(Red,lib=c(1,500),pred=c(501,1000),E=c(2:8))
# sim_l <- simplex(Logi,lib=c(1,500),pred=c(501,1000),E=c(2:8))
# 
# # i=c(2:8)
# # ## Plot predictive skill (rho) vs embedding dimension (E) - gives error that rho cannot be of type list
# # par(mfrow=c(2,1),mar=c(4,4,1,1))
# # plot(rho[i]~i,data=sim_r,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.3,0.4),col=2,main="Red noise")
# # plot(rho~E,data=sim_l,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.95,1.02),col=4,main="Logistic map")
# 
# ## The optimal embedding dimension determined by maximizing rho
# (E_r <-sim_r[which.max(sim_r$rho),"E"][1])# The optimal E of red noise
# (E_l <-sim_l[which.max(sim_l$rho),"E"][1])# The optimal E of logistic map
# 
# # S map for Red Noise & logistic map
# smap_r <- s_map(Red,E=E_r,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))
# smap_l <- s_map(Logi,E=E_l,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))
# 
# ## Plot predictive skill (rho) vs state-dependency parameter (theta) 
# plot(rho~theta,data=smap_r,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.4,0.5),col=2,main="Red noise")
# plot(rho~theta,data=smap_l,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.6,1),col=4,main="Logistic map")
# 
# ## The optimal theta determined by maximizing rho
# (the_r <- smap_r[which.max(smap_r$rho),"theta"][1])
# (the_l <- smap_l[which.max(smap_l$rho),"theta"][1])
# 
# 
# # time series of Moran effect model --------------------------------------------------------------------
# # Loading the time series for the Moran effect and mirage correlation models 
# dam <- read.csv('ESM3_Data_moran.csv',header=T) # Moran effect
# dac <- read.csv('ESM4_Data_competition.csv',header=T) # Mirage correlation
# 
# # Data normalization
# dac.n <- scale(dac[,-1], center = TRUE, scale = TRUE)
# dam.n <- scale(dam[,-1], center = TRUE, scale = TRUE)
# 
# ## CCM analysis of the Moran effect model, N1 and N2
# # Design a sequence of library size
# libs <- c(seq(20,80,20),seq(100,1000,100))
# 
# # Moran effect model: N1 cross-mapping N2 (i.e. testing N2 as a cause of N1)
# # Determine the embedding dimension - gives error
# E.test.n1=NULL
# for(E.t in 2:8){
#   cmxy.t <- ccm(dam.n, E = E.t, lib_column = "N1", target_column = "N2",
#                 lib_sizes = 1000, num_samples = 1, tp=-1,random_libs = F)
#   E.test.n1=rbind(E.test.n1,cmxy.t)}
# (E_n1 <- E.test.n1$E[which.max(E.test.n1$rho)[1]]) # the optimal E
# 
# n1_xmap_n2 <- ccm(dam.n, E=E_n1,lib_column="N1", target_column="N2",
#                   lib_sizes=libs, num_samples=200, replace=T, RNGseed=2301)
# 
# # Calculate the median, maximum, and 1st & 3rd quantile of rho for each L
# n12q=as.matrix(aggregate(n1_xmap_n2[,c('rho')],by = list(as.factor(n1_xmap_n2$lib_size)), quantile)[,'x'])
# apply(n12q[,2:5],2,MannKendall)
# 
# ###########################################################
# # Moran effect model: N2 cross-mapping N1 (i.e. testing N1 as a cause of N2)
# # Determine the embedding dimension
# E.test.n2=NULL
# for(E.t in 2:8){
#   cmxy.t <- ccm(dam.n, E = E.t, lib_column = "N2", target_column = "N1",
#                 lib_sizes = 1000, num_samples = 1, tp=-1, random_libs = F)
#   E.test.n2=rbind(E.test.n2,cmxy.t)}
# (E_n2 <- E.test.n2$E[which.max(E.test.n2$rho)[1]])
# 
# # CCM analysis
# n2_xmap_n1 <- ccm(dam.n, E=E_n2,lib_column="N2", target_column="N1",
#                   lib_sizes=libs, num_samples=200, replace=T, RNGseed=2301)
# 
# # Calculate the (25%,50%,75%,100%) quantile for predictive skills
# n21q=as.matrix(aggregate(n2_xmap_n1[,c('rho')],by = list(as.factor(n2_xmap_n1$lib_size)), quantile)[,'x'])
# apply(n21q[,2:5],2,MannKendall)
# 
# # Plot forecast skill vs library size
# # Plot N1 cross-mapping N2
# plot(n12q[,3]~libs,type="l",col="red",ylim=c(0,1),lwd=2,
#      main="Convergent cross mapping CCM",xlab="Library size",ylab=expression(rho)) # median predictive skill vs library size (or we can use mean predictive skill)
# lines(n12q[,2]~libs,col="red",lwd=1,lty=2) # 1st quantile 
# lines(n12q[,4]~libs,col="red",lwd=1,lty=2) # 3rd quantile



## Load package and data
d <- read.csv("ESM5_Data_5spModel.csv")
# Please reduce the number of data points if the calculation needs long time
data_used <- 1:1000

# Specify the length of time series to be used to reconstruct state space (Library length)
lib_point <- c(1,floor(max(data_used)/2))

# Specify which points will be predicted based on the reconstructed state space
pred_point <- c(floor(max(data_used)/2)+1, max(data_used))

# Time series of C1 is normalized
C1 <- as.numeric(scale(d[data_used,'C1']))

# Estimate the best embedding dimension
simp_C1_tmp <- simplex(C1, E=1:10, silent = T)
plot(simp_C1_tmp$E, simp_C1_tmp$mae, type="l", xlab="E", ylab="MAE")

# Best E = 3
bestE_C1 <- simp_C1_tmp[which.min(simp_C1_tmp$mae),"E"]

# Perform univariate simplex projection
# We need to specify time series (C1), embedding dimension (E), library length (lib), predictee (pred) and which output we need (stats_only). If you do not want to see warning message, "silent" option should be set as "T".
simp_C1 <- simplex(C1, E=bestE_C1, lib=lib_point, pred=pred_point, stats_only = F, silent = T)
C1_pred_uni <- na.omit(simp_C1$model_output$E3$Predictions)
C1_obs_uni <- na.omit(simp_C1$model_output$E3$Observations)
plot(C1_obs_uni, C1_pred_uni, xlab="Observed", ylab="Predicted")
abline(0,1) # add 1:1 line

# Make multivariate embedding
Embedding <- c("C1", "R", "P1")
block <- d[,Embedding]


block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))

# Do multivariate simplex projection using block_lnlp() function
# We need to specify time series, method (simplex or s-map), library length (lib), predictee (pred) and which output we need (stats_only).
mult_simp_C1 <-  block_lnlp(block[data_used,], method = "simplex", lib = lib_point, pred = pred_point,stats_only = F, silent = T)

C1_pred_mult <- na.omit(mult_simp_C1$model_output$Predictions)
C1_obs_mult <- na.omit(mult_simp_C1$model_output$Observations)
plot(C1_obs_mult, C1_pred_mult, xlab="Observed", ylab="Predicted")
abline(0,1) # add 1:1 line
