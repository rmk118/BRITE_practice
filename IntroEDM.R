#Learning about EDM
#Ruby Krasnow
#5/27/23

library(rEDM)
library(Kendall)
library(dplyr)
library(ggplot2)
library(ggplotify)
library(gridExtra)


# TentMap example from vignette (univariate) -------------------------------------------
str(TentMap) #view data

simplex_out <- Simplex(dataFrame = TentMap, lib = "1 100", pred = "201 500", columns = "TentMap",target = "TentMap", E = 3) #simplex projection using known E

simplex_out
simplex_out[c(1:2, 300:301), ] #it would use Time 200 to predict the value at Time 201, but that's not included in the set so Pred at 201 is NaN

ComputeError(simplex_out$Observations, simplex_out$Predictions) #rho is 0.94

rho_E <- EmbedDimension(dataFrame = TentMap, lib = "1 100", pred = "201 500", columns = "TentMap",target = "TentMap") #determines optimal embedding dimension by calling simplex fn with E values from 1 to 10 and determining max. rho value

rho_E[which.max(rho_E$rho),"E"][1] #2

rho_Tp <- PredictInterval(dataFrame = TentMap, lib = "1 100", pred = "201 500", target = "TentMap",columns = "TentMap", E = 2) #forecasting skill declines as time interval increases

rho_theta <- PredictNonlinear(dataFrame = TentMapNoise, lib = "1 100", pred = "201 500",
                              target = "TentMap", columns = "TentMap", E = 2) #evaluates S-map forecasting skill as theta increases. Since rho isn't maximized at theta=0 and increases as theta increases, we have evidence of nonlinear deterministic behavior. We also see rho decline at high values of theta as the local linear map overfits to insufficient nearest neighbors



# 3-species example from vignette (multivariate) -------------------------------------------------

head(block_3sp,3) #preview data
smplx_3species = Simplex(dataFrame = block_3sp, lib = "1 100", pred = "101 190",
                         E = 3, columns = "x_t x_t-1 z_t", target = "x_t", embedded = TRUE) #by default, simplex does time-lagging on univariate or multivariate data. We have time-lagged data already (E=3, tau = -1), so we need to set embedded=TRUE so the input data are assumed to constitute a valid multidimensional embedding and no time-delay embedding will be performed

err = ComputeError(smplx_3species$Observations, smplx_3species$Predictions)

plot(smplx_3species$Observations, smplx_3species$Predictions, pch = 19, cex = 0.5,
     xlab = "Observations", ylab = "Predictions", main = "3 Species x_t")
abline(a = 0, b = 1, lty = 2, col = "blue")
text(-1, 1, paste(capture.output(cbind(err)), collapse = "\n"))

#Construct all possible embeddings of dimension E with lag up to E-1, ranked by rho over the library portion of the data; individual forecasts for the top m embeddings are then averaged together, where m=sqrt(C) and C is the # of E-dimensional combinations created from all data vectors.
Mview = Multiview(dataFrame = block_3sp, lib = "1 100", pred = "101 190", E = 3,
                  columns = "x_t y_t z_t", target = "x_t")

Mview$View
Mview$Predictions #final averaged multiview projections
Mview$View[which(Mview$View$rho > 0.91), ]


# CCM anchovy example from vignette (SST affects anchovy pop but not vice verse)
cmap <- CCM(dataFrame = sardine_anchovy_sst, E = 3, Tp = 0, columns = "anchovy",
            target = "np_sst", libSizes = "10 70 5", sample = 100, showPlot = TRUE)


# Apple pest example from vignette (multivariate) ----------------------------------------

head(Thrips, 2) #view data
rho_E <- EmbedDimension(dataFrame = Thrips, columns = "Thrips_imaginis", target = "Thrips_imaginis",lib = "1 72", pred = "1 72", showPlot = TRUE)

E = 8 #although could also use 3

#test for nonlinearity
rho_theta_e3 = PredictNonlinear(dataFrame = Thrips, columns = "Thrips_imaginis",
                                target = "Thrips_imaginis", lib = "1 73", pred = "1 73", E = E) #clear non-linearity, indicates doesn't passively track seasonality

#create matrix showing cross-mapping of all 2-variable pairs
vars = colnames(Thrips[3:6])
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time
libSize = paste(NROW(Thrips) - E, NROW(Thrips) - E, 10, collapse = " ")
ccm_matrix = array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,
                                                                            vars))
for (i in 1:ncol(var_pairs)) {
  ccm_out = CCM(dataFrame = Thrips, columns = var_pairs[1, i], target = var_pairs[2,
                                                                                  i], libSizes = libSize, Tp = 0, E = E, sample = 100)
  outVars = names(ccm_out)
  var_out = unlist(strsplit(outVars[2], ":"))
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2]
  var_out = unlist(strsplit(outVars[3], ":"))
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3]
}

#create correlation matrix for comparison
corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,
                                                                              vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(Thrips[, ccm_from], Thrips[, ccm_to], type = "correlation",
                   lag.max = 6, plot = FALSE)$acf
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out)) }}

ccm_matrix #non-symmetric
corr_matrix #symmetric
    
    
thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0, columns = "Thrips_imaginis",
                        target = "maxT_degC", libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix["Thrips_imaginis", "maxT_degC"], col = "black", lty = 2)

thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0, columns = "Thrips_imaginis",
                        target = "Rain_mm", libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix["Thrips_imaginis", "Rain_mm"], col = "black", lty = 2)

thrips_xmap_maxT <- CCM(dataFrame = Thrips, E = E, Tp = 0, columns = "Thrips_imaginis",
                        target = "Season", libSizes = "13 73 3", sample = 300, showPlot = TRUE)
abline(h = corr_matrix["Thrips_imaginis", "Season"], col = "black", lty = 2)

# Create matrix with temperature and rain surrogates (1000 time series vectors)
surr_maxT = SurrogateData(Thrips$maxT_degC, method = "seasonal", T_period = 12, num_surr = 1000,
                          alpha = 3)
surr_rain = SurrogateData(Thrips$Rain_mm, method = "seasonal", T_period = 12, num_surr = 1000,
                          alpha = 3)

# Rain cannot be negative
surr_rain = apply(surr_rain, 2, function(x) {
  i = which(x < 0)
  x[i] = 0
  x
})

# data.frame to hold CCM rho values between Thrips abundance and variable
rho_surr <- data.frame(maxT = numeric(1000), Rain = numeric(1000))

# data.frames with time, Thrips, and 1000 surrogate climate variables for CCM()
maxT_data = as.data.frame(cbind(seq(1:nrow(Thrips)), Thrips$Thrips_imaginis, surr_maxT))
names(maxT_data) = c("time", "Thrips_imaginis", paste("T", as.character(seq(1, 1000)),
                                                      sep = ""))
rain_data = as.data.frame(cbind(seq(1:nrow(Thrips)), Thrips$Thrips_imaginis, surr_rain))
names(rain_data) = c("time", "Thrips_imaginis", paste("R", as.character(seq(1, 1000)),
                                                      sep = ""))
# Cross mapping
for (i in 1:1000) {
  targetCol = paste("T", i, sep = "") # as in maxT_data
  ccm_out = CCM(dataFrame = maxT_data, E = E, Tp = 0, columns = "Thrips_imaginis",
                target = targetCol, libSizes = "73 73 5", sample = 1)
  col = paste("Thrips_imaginis", ":", targetCol, sep = "")
  rho_surr$maxT[i] = ccm_out[1, col]
}
for (i in 1:1000) {
  targetCol = paste("R", i, sep = "") # as in rain_data
  ccm_out = CCM(dataFrame = rain_data, E = E, Tp = 0, columns = "Thrips_imaginis",
                target = targetCol, libSizes = "73 73 5", sample = 1)
  col = paste("Thrips_imaginis", ":", targetCol, sep = "")
  rho_surr$Rain[i] = ccm_out[1, col]
}
1 - ecdf(rho_surr$maxT)(ccm_matrix["maxT_degC", "Thrips_imaginis"])
1 - ecdf(rho_surr$Rain)(ccm_matrix["Rain_mm", "Thrips_imaginis"])



# Sunspot example from GitHub page (univariate) ---------------------------

df = data.frame(yr = as.numeric(time(sunspot.year)), 
                sunspot_count = as.numeric(sunspot.year))

plot(df$yr, df$sunspot_count, type = "l", xlab = "year", ylab = "sunspots") #time series
rho_E <- EmbedDimension(dataFrame = df, columns = "sunspot_count", target = "sunspot_count",lib = "1 280", pred = "1 280", showPlot = TRUE)
(E_r <-rho_E[which.max(rho_E$rho),"E"][1])# The optimal E of sunspots

simplex_out <- Simplex(dataFrame = df, lib = "1 190", pred = "191 287", columns = "sunspot_count", E = 3)

plot( df$yr, df$sunspot_count, type = "l", lwd = 2,
      xlab = "year", ylab = "sunspots")
lines( simplex_out$yr, simplex_out$Predictions, col = "red", lwd = 2)

ComputeError(simplex_out$Observations, simplex_out$Predictions)


# time series of red noise and logistic map - my way --------------------------------------------------------------------
dat2 <- read.csv('ESM2_Data_noise.csv',header=T)
dat2 <- dat2 %>% mutate(index=c(1:1000), .before=R)

# Data normalization
dat2 <- dat2 %>% mutate(
  Red2 = (R-mean(R))/sd(R),
  Logi2 = (L-mean(L))/sd(L), .keep="unused"
)

red_time_series<-ggplot(dat2, aes(x=index, y=Red2))+geom_line()+theme_classic() #mine
red_time_series #plot time series

logi_time_series<-ggplot(dat2, aes(x=index, y=Logi2))+geom_line()+theme_classic()+ylim(-2, 2) #mine
logi_time_series #plot time series


rho_E2 <- EmbedDimension(dataFrame = dat2, columns = "Red2", target = "Red2",lib = "1 500", pred = "501 1000", maxE=8, showPlot = TRUE)
optimalRhoRed<-as.ggplot(~plot(rho_E2, type="l"))
rho_E2[which.max(rho_E2$rho), "E"] #7 for red

rho_E3 <- EmbedDimension(dataFrame = dat2, columns = "Logi2", target = "Logi2",lib = "1 500", pred = "501 1000", maxE=8, showPlot = TRUE)
optimalRhoLogi<-as.ggplot(~plot(rho_E3, ylim=c(0.95,1.02), type="l"))
rho_E3[which.max(rho_E3$rho), "E"] #1 for logi

#test for nonlinearity
rho_thetaR = PredictNonlinear(dataFrame = dat2, columns = "Red2",
                              target = "Red2", lib = "1 500", pred = "501 1000", theta=seq(0,2,0.1), E = 7)
rho_vs_thetaR<-as.ggplot(~plot(rho_thetaR, type="l",ylim=c(0.4,0.5)))

rho_thetaL = PredictNonlinear(dataFrame = dat2, columns = "Logi2",
                              target = "Logi2", lib = "1 500", pred = "501 1000", theta=seq(0,2,0.1), E = 2)
rho_vs_thetaL<-as.ggplot(~plot(rho_thetaL, type="l",ylim=c(0.6,1)))

grid.arrange(red_time_series, logi_time_series, optimalRhoRed, optimalRhoLogi, rho_vs_thetaR, rho_vs_thetaL, ncol=2) #SUCCESS

# # time series of red noise and logistic map - error way --------------------------------------------------------------------
# dat <- read.csv('ESM2_Data_noise.csv',header=T)
# 
# # Data normalization 
# Red <- ((dat[,"R"]-mean(dat[,"R"]))/sd(dat[,"R"]))
# Logi <- ((dat[,"L"]-mean(dat[,"L"]))/sd(dat[,"L"]))
# 
# plot(Red, type="l") #plot time series
# 
# # Simplex projection for red noise and logistic map
# sim_r <- simplex(Red,lib=c(1,500),pred=c(501,1000),E=c(2:8))
# sim_l <- simplex(Logi,lib=c(1,500),pred=c(501,1000),E=c(2:8))
# 
# ## Plot predictive skill (rho) vs embedding dimension (E) 
# par(mfrow=c(2,1),mar=c(4,4,1,1))
# plot(rho~E,data=sim_r,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.3,0.4),col=2,main="Red noise")
# plot(rho~E,data=sim_l,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.95,1.02),col=4,main="Logistic map")
# 
# #S map for Red Noise & logistic map
# smap_r <- s_map(Red,E=7,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))
# smap_l <- s_map(Logi,E=2,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))
# 
# # Plot predictive skill (rho) vs state-dependency parameter (theta)
# plot(rho~theta,data=smap_r,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.4,0.5),col=2,main="Red noise")
# plot(rho~theta,data=smap_l,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.6,1),col=4,main="Logistic map")
# 
# # The optimal theta determined by maximizing rho
# (the_r <- smap_r[which.max(smap_r$rho),"theta"][1])
# (the_l <- smap_l[which.max(smap_l$rho),"theta"][1])


# # time series of Moran effect model - error way --------------------------------------------------------------------
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
# 
# ## Load package and data
# d <- read.csv("ESM5_Data_5spModel.csv")
# # Please reduce the number of data points if the calculation needs long time
# data_used <- 1:1000
# 
# # Specify the length of time series to be used to reconstruct state space (Library length)
# lib_point <- c(1,floor(max(data_used)/2))
# 
# # Specify which points will be predicted based on the reconstructed state space
# pred_point <- c(floor(max(data_used)/2)+1, max(data_used))
# 
# # Time series of C1 is normalized
# C1 <- as.numeric(scale(d[data_used,'C1']))
# 
# # Estimate the best embedding dimension
# simp_C1_tmp <- simplex(C1, E=1:10, silent = T)
# plot(simp_C1_tmp$E, simp_C1_tmp$mae, type="l", xlab="E", ylab="MAE")
# 
# # Best E = 3
# bestE_C1 <- simp_C1_tmp[which.min(simp_C1_tmp$mae),"E"]
# 
# # Perform univariate simplex projection
# # We need to specify time series (C1), embedding dimension (E), library length (lib), predictee (pred) and which output we need (stats_only). If you do not want to see warning message, "silent" option should be set as "T".
# simp_C1 <- simplex(C1, E=bestE_C1, lib=lib_point, pred=pred_point, stats_only = F, silent = T)
# C1_pred_uni <- na.omit(simp_C1$model_output$E3$Predictions)
# C1_obs_uni <- na.omit(simp_C1$model_output$E3$Observations)
# plot(C1_obs_uni, C1_pred_uni, xlab="Observed", ylab="Predicted")
# abline(0,1) # add 1:1 line
# 
# # Make multivariate embedding
# Embedding <- c("C1", "R", "P1")
# block <- d[,Embedding]
# 
# 
# block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))
# 
# # Do multivariate simplex projection using block_lnlp() function
# # We need to specify time series, method (simplex or s-map), library length (lib), predictee (pred) and which output we need (stats_only).
# mult_simp_C1 <-  block_lnlp(block[data_used,], method = "simplex", lib = lib_point, pred = pred_point,stats_only = F, silent = T)
# 
# C1_pred_mult <- na.omit(mult_simp_C1$model_output$Predictions)
# C1_obs_mult <- na.omit(mult_simp_C1$model_output$Observations)
# plot(C1_obs_mult, C1_pred_mult, xlab="Observed", ylab="Predicted")
# abline(0,1) # add 1:1 line
