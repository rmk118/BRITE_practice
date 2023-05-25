#Learning about EDM
#Ruby Krasnow
#5/25/23

library(rEDM)
library(Kendall)
## Data loading 


# time series of red noise and logistic map --------------------------------------------------------------------
dat <- read.csv('ESM2_Data_noise.csv',header=T)

plot(dat$R)
Red <- ((dat[,"R"]-mean(dat[,"R"]))/sd(dat[,"R"]))
plot(Red)


plot(dat$L)
Logi <- ((dat[,"L"]-mean(dat[,"L"]))/sd(dat[,"L"]))
plot(Logi)

## Simplex projection for red noise and logistic map
sim_r <- simplex(Red,lib=c(1,500),pred=c(501,1000),E=c(2:8))
sim_l <- simplex(Logi,lib=c(1,500),pred=c(501,1000),E=c(2:8))

i=c(2:8)
## Plot predictive skill (rho) vs embedding dimension (E) - gives error
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(rho[i]~i,data=sim_r,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.3,0.4),col=2,main="Red noise")
plot(rho~E,data=sim_l,type="l",xlab="Embedding dimension (E)",ylab=expression(rho),ylim=c(0.95,1.02),col=4,main="Logistic map")

## The optimal embedding dimension determined by maximizing rho
(E_r <-sim_r[which.max(sim_r$rho),"E"][1])# The optimal E of red noise
(E_l <-sim_l[which.max(sim_l$rho),"E"][1])# The optimal E of logistic map

# S map for Red Noise & logistic map
smap_r <- s_map(Red,E=E_r,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))
smap_l <- s_map(Logi,E=E_l,lib=c(1,500),pred=c(501,1000),theta=seq(0,2,0.1))

## Plot predictive skill (rho) vs state-dependency parameter (theta) 
plot(rho~theta,data=smap_r,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.4,0.5),col=2,main="Red noise")
plot(rho~theta,data=smap_l,type="l",xlab=expression(theta),ylab=expression(rho),ylim=c(0.6,1),col=4,main="Logistic map")

## The optimal theta determined by maximizing rho
(the_r <- smap_r[which.max(smap_r$rho),"theta"][1])
(the_l <- smap_l[which.max(smap_l$rho),"theta"][1])


# time series of Moran effect model --------------------------------------------------------------------
# Loading the time series for the Moran effect and mirage correlation models 
dam <- read.csv('ESM3_Data_moran.csv',header=T) # Moran effect
dac <- read.csv('ESM4_Data_competition.csv',header=T) # Mirage correlation


