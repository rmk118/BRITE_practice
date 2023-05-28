#Learning about EDM
#Ruby Krasnow
#5/28/23

library(rEDM)
library(Kendall)
library(dplyr)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(patchwork)


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

red_time_series<-ggplot(dat2, aes(x=index, y=Red2))+geom_line()+theme_classic()+ylab("Standardized density")+xlab("Time")+ggtitle("Red noise") #mine
red_time_series #plot time series

logi_time_series<-ggplot(dat2, aes(x=index, y=Logi2))+geom_line()+theme_classic()+ylim(-2, 2)+ylab("")+xlab("Time")+ggtitle("Logistic map") #mine
logi_time_series #plot time series

rho_E2 <- EmbedDimension(dataFrame = dat2, columns = "Red2", target = "Red2",lib = "1 500", pred = "501 1000", maxE=8, showPlot = FALSE)
optimalRhoRed<-(as.ggplot(~plot(rho_E2, type="l", ylab=expression(rho), xlab="E")))
optimalRhoRed
rho_E2[which.max(rho_E2$rho), "E"] #7 for red

rho_E3 <- EmbedDimension(dataFrame = dat2, columns = "Logi2", target = "Logi2",lib = "1 500", pred = "501 1000", maxE=8, showPlot = TRUE)
optimalRhoLogi<-as.ggplot(~plot(rho_E3, ylim=c(0.95,1.02), type="l", ylab=""))
rho_E3[which.max(rho_E3$rho), "E"] #1 for logi

#test for nonlinearity
rho_thetaR = PredictNonlinear(dataFrame = dat2, columns = "Red2",
                              target = "Red2", lib = "1 500", pred = "501 1000", theta=seq(0,2,0.1), E = 7)
rho_vs_thetaR<-as.ggplot(~plot(rho_thetaR, type="l",ylim=c(0.4,0.5), ylab=expression(rho), xlab=expression(theta)))

rho_thetaL = PredictNonlinear(dataFrame = dat2, columns = "Logi2",
                              target = "Logi2", lib = "1 500", pred = "501 1000", theta=seq(0,2,0.1), E = 2)
rho_vs_thetaL<-as.ggplot(~plot(rho_thetaL, type="l",ylim=c(0.6,1), ylab="",xlab=expression(theta)))

grid.arrange(red_time_series, logi_time_series, optimalRhoRed, optimalRhoLogi, rho_vs_thetaR, rho_vs_thetaL, ncol=2) #SUCCESS

allGraphs<- (red_time_series / optimalRhoRed / rho_vs_thetaR) | (logi_time_series / optimalRhoLogi / rho_vs_thetaL)
allGraphs

#time series of Moran effect model - my way --------------------------------------------------------------------
# Loading the time series for the Moran effect and mirage correlation models
dam2.pre <- read.csv('ESM3_Data_moran.csv',header=T) # Moran effect

dam2 <- scale(dam2[,-1], center = TRUE, scale = TRUE)# Data normalization

damShort<- dam2 %>% filter(Time<101)

#plot time series
moran_series<-ggplot(damShort, aes(x=Time))+
  geom_line(aes(y = N1, color = "N1"))+ theme_classic()+
  geom_line(aes(y = N2, color = "N2"))+
  ylab(expression("Density (capita*m" ^"-3"~")"))+xlab("Time")+ggtitle("Moran effect")+
  scale_color_manual(values=c("Black", "Red"))+ theme(legend.title = element_blank())
moran_series

#Determine optimal embedding dimension
rho_N1 <- EmbedDimension(dataFrame = dam2, columns = "N2", target = "N1",lib = "1 500", pred = "501 1000", maxE=8, showPlot = TRUE)
rho_N1[which.max(rho_N1$rho), "E"] #5 for N1

rho_N2 <- EmbedDimension(dataFrame = dam2, columns = "N1", target = "N2",lib = "1 500", pred = "501 1000", maxE=8, showPlot = TRUE)
rho_N2[which.max(rho_N2$rho), "E"] #6 for N2, but 5 is almost as good

# CCM analysis of the Moran effect model, N1 and N2
CCM(dataFrame = dam2, E = 5, Tp = 0, columns = "N1",
              target = "N2", libSizes = "100 1000 200", sample = 100, showPlot = TRUE)


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

#time series of Moran effect model - my way --------------------------------------------------------------------
# Loading the time series for the Moran effect and mirage correlation models
dam2.pre <- read.csv('ESM3_Data_moran.csv',header=T) # Moran effect

dam2 <- scale(dam2.pre[,-1], center = TRUE, scale = TRUE)# Data normalization

# Data normalization
dam3 <- dam2.pre %>% mutate(
  N1 = (N1-mean(N1))/sd(N1),
  R1 = (R1-mean(R1))/sd(R1),
  R2 = (R2-mean(R2))/sd(R2),
  N2 = (N2-mean(N2))/sd(N2)
)

damShort<- dam3 %>% filter(Time<101)

#plot time series
moran_series<-ggplot(damShort, aes(x=Time))+
  geom_line(aes(y = N1, color = "N1"))+ theme_classic()+
  geom_line(aes(y = N2, color = "N2"))+
  ylab(expression("Density (capita*m" ^"-3"~")"))+xlab("Time")+ggtitle("Moran effect")+
  scale_color_manual(values=c("Black", "Red"))+ theme(legend.title = element_blank())
moran_series
