# Classic SEIR model for transmission dynamics of 2019-nCoV
# β is the average number of infected individuals per infectious subject per unit time, 
# α is the reciprocal of average latent period, 
# v is the rate of recovery (reciprocal of duration of the infection).


# Wuhan
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
seir = function(N = 1000, SUS = .5*N,  EXP=0, INF = 1,IMM = N-SUS-EXP-INF, days = 200, beta=0.8, alpha=1/6, mu=1/10, 
                TMP = 0.04, CLS=1, IMP=10){
  M = data.frame(ID = 1:days,SUS = NA, EXP = NA, INF = NA, IMM = NA)
  M[1,] = c(1, SUS, EXP, INF, IMM)
  for(i in 2:days){
    if (i > 90)  beta <- beta * (1-TMP)^(i-90) # temperature effect Wuhan
    # if (i > 40)  beta <- beta * (1-TMP)^(i-40) # temperature effect Beijing
    if (i < (31+23))  { # city closure effect
    M$SUS[i] = M$SUS[i-1] - beta * M$INF[i-1] * M$SUS[i-1]/N 
    M$EXP[i] = M$EXP[i-1] + beta * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
    M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu # + IMP[i] # imported cases every day
    M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
    else { # city closure effect
      beta2 = beta * CLS
      M$SUS[i] = M$SUS[i-1] - beta2 * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta2 * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu # + IMP[i] # imported cases every day
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
  }
  return(invisible(M))
}
# WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW


# Beijing
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
seir = function(N = 1000, SUS = .5*N,  EXP=0, INF = 1,IMM = N-SUS-EXP-INF, days = 200, beta=0.8, alpha=1/6, mu=1/10, 
                TMP = 0.04, CLS=1, IMP=10){
  M = data.frame(ID = 1:days,SUS = NA, EXP = NA, INF = NA, IMM = NA)
  M[1,] = c(1, SUS, EXP, INF, IMM)
  for(i in 2:days){
    # if (i > 90)  beta <- beta * (1-TMP)^(i-90) # temperature effect Wuhan
    if (i > 40)  beta <- beta * (1-TMP)^(i-40) # temperature effect Beijing
    if (i < (31+23))  { # city closure effect
      M$SUS[i] = M$SUS[i-1] - beta * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu  + IMP[i] # imported cases every day
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
    else { # city closure effect
      beta2 = beta * CLS
      M$SUS[i] = M$SUS[i-1] - beta2 * M$INF[i-1] * M$SUS[i-1]/N 
      M$EXP[i] = M$EXP[i-1] + beta2 * M$INF[i-1] * M$SUS[i-1]/N - alpha*M$EXP[i-1] 
      M$INF[i] = M$INF[i-1] + alpha* M$EXP[i-1] - M$INF[i-1]*mu  + IMP[i] # imported cases every day
      M$IMM[i] = M$IMM[i-1] + M$INF[i-1]*mu
    }
  }
  return(invisible(M))
}
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB






city = read.csv("city_time_series.csv", header=T)
head(city)
WH = as.numeric(city[, c("Wuhan")])
WH = WH[!is.na(WH)]
WH = c(41,45,62,121,198,270,375,444,444,549,572,618,698,1590,1905,2261,3215,4109,5142, 6384, 8351, 10117) #Jan15-Feb5

# Match confirmed cases for Wuhan: infectious R0 = beta/mu/(S/N) = 2.3*5/2 = 5.75 
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir(N=10000000, beta=2.3, alpha=1/6, mu=1/5) # Wuhan
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,50), ylim=c(0, 9000),xlab="Days",ylab="Number of infected cases")
points(c(46:67)-22, WH, pch=16, col="blue")

beta=2.3
for (i in 1:200){
  B = rnorm(1, beta, beta*0.15)
  DD = seir(N=10000000, beta=B, alpha=1/6, mu=1/5) # Wuhan
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}
DD = seir(N=10000000, beta=2.3, alpha=1/6, mu=1/5) # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(46:67)-22, WH, pch=16, col="blue")
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# real situation for Wuhan (R0 = 5)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir(N=10000000, beta=10/5, alpha=1/6, mu=1/5, days=200, CLS=.5) # Wuhan
# plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,200), ylim=c(0, 10000000),xlab="Days since Dec. 1, 2019",ylab="Number of cases")
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(30, 150), ylim=c(0, 600000),xlab="Days since Dec. 1, 2019",ylab="Number of infected cases")
points(c(46:67) + 0, WH, pch=16, col="blue", cex=.5)
abline(v=60, lty=2)

beta=10/5
for (i in 1:200){
  B = rnorm(1, beta, beta*0.15)
  DD = seir(N=10000000, beta=B, alpha=1/6, mu=1/5, CLS=.5) # Wuhan
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
  lines(DD$ID, DD$SUS, col=adjustcolor("black", alpha.f=0.2))
  lines(DD$ID, DD$IMM, col=adjustcolor("green", alpha.f=0.2))
  lines(DD$ID, DD$EXP, col=adjustcolor("yellow",  alpha.f=0.2))
  peak[i] = DD[DD$INF == max(DD$INF), c('ID')]
}

mean(peak);sd(peak) # peak of the pandemic wave
DD = seir(N=10000000, beta=10/5, alpha=1/6, mu=1/5, CLS=.5) # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(46:67)+0, WH, pch=16, col="blue", cex=.8)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# temperature effect for Wuhan
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
DD = seir(N=10000000, beta=10/5, alpha=1/6, mu=1/5, days=200, CLS=.5, TMP=0.1) # Wuhan
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,160), ylim=c(0, 600000),xlab="Days since Dec. 1, 2019",ylab="Number of infected cases")
# points(c(46:67) + 0, WH, pch=16, col="blue", cex=.5)

TMP=0.1
for (i in 1:200){
  TMP = rnorm(1, TMP, TMP*0.3)
  DD = seir(N=10000000, beta=10/5, TMP=TMP, alpha=1/6, mu=1/5, CLS=.5) # Wuhan
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
}

DD = seir(N=10000000, beta=10/5, alpha=1/6, mu=1/5, CLS=.5, TMP=0) # Wuhan
lines(DD$ID, DD$INF, col="red", lwd=2)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE





city = read.csv("city_time_series.csv", header=T)
BJ = as.numeric(city[, c("Beijing")]) # Jan. 19
BJ = BJ[!is.na(BJ)]
BJ = BJ[!BJ==0]

# Matching confirmed cases for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; imp; mean(imp); sd(imp)
DD = seir(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,25), ylim=c(0, 400),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:21) + 0, BJ, pch=16, col="blue", cex=.5)

for (i in 1:200){
  imp = imp + rnorm(1, 0, 2) # SD 2 was based on the daily report of confirmed cases in Beijing
  DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing
  DD = DD[!DD$INF[5]<0,]
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.1))
}

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; imp; mean(imp); sd(imp)
DD = seir(N=22000000, beta=2.8/5, alpha=1/6, mu=1/5, days=30, IMP=imp) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:21) + 0, BJ, pch=16, col="blue", cex=1)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# Simulate imported cases for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 170))
DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,150), ylim=c(0, 2000),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:21) + 0, BJ, pch=16, col="blue", cex=.5)

for (i in 1:200){
  imp = imp + rnorm(1, 0, 3.85)
  DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
  DD = DD[!DD$INF[5]<0,]
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
}

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5;imp = c(imp, rep(0, 170))
DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:18) + 0, BJ, pch=16, col="blue", cex=.5)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


# Simulate R0 for Beijing
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5; 
imp = c(imp, rep(0, 170))
DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
plot(DD$ID, DD$INF, col="red", type="l", xlim=c(0,150), ylim=c(0, 10000),xlab="Days since Jan. 19, 2020",ylab="Number of infected cases")
points(c(1:length(BJ)) + 0, BJ, pch=16, col="blue", cex=.5)

peak.day = numeric(200)
peak = numeric(200)
beta=3/5
TMP=0.001
for (i in 1:200){
  # imp = imp + rnorm(1, 0, 3.85)
  B = rnorm(1, beta, beta*0.15)
  tmp = rnorm(1, TMP, TMP*0.4)
  DD = seir(N=22000000, beta=B, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=tmp) # Beijing
  peak.day[i] = DD[DD$INF == max(DD$INF), c('ID')]
  peak[i] = DD[DD$INF == max(DD$INF), c('INF')]
  DD$INF[DD$INF<0] <- 0
  lines(DD$ID, DD$INF, col=adjustcolor("red", alpha.f=0.2))
}
mean(peak);sd(peak) # peak of the pandemic wave in Beijing

imp = floor(dnorm(1:30, 15, 5) *300)*1.5 +5;imp = c(imp, rep(0, 170))
DD = seir(N=22000000, beta=3/5, alpha=1/6, mu=1/5, days=200, IMP=imp, TMP=0.001) # Beijing
lines(DD$ID, DD$INF, col="red", lwd=2)
points(c(1:18) + 0, BJ, pch=16, col="blue", cex=.5)
# EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE



