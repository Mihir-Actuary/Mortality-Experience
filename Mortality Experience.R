###############################################################################################
###########################          Mortality Experience Model Fixed       #################
###########################        LEE-CARTER FORECAST TO 2040                   ############
###############################################################################################

##########     Libraries

#install.packages("splines")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("fields")

library(fields)
library(splines)
library(ggplot2)
library(reshape2)


###############################################################################################
##########     Working Directory and Import Data
###############################################################################################

setwd("C:\\Users\\mihir\\Documents\\Mihir Documents\\Projects\\Mortality Excel")

Census_Data_C = read.csv("Age By Census Year and Gender CSV - 1996-2011.csv", skip=12, header=FALSE)
Age_Gender_Population_2022 = read.csv("Gender - Age Population CSV - 2022.csv", skip=11, header=FALSE)

Age_Gender_2006 = read.csv("Gender - Age - 2006 CSV.csv", skip=11, header=FALSE)
Age_Gender_2008 = read.csv("Gender - Age - 2008 CSV.csv", skip=11, header=FALSE)
Age_Gender_2009 = read.csv("Gender - Age - 2009 CSV.csv", skip=11, header=FALSE)
Age_Gender_2010 = read.csv("Gender - Age - 2010 CSV.csv", skip=11, header=FALSE)
Age_Gender_2011 = read.csv("Gender - Age - 2011 CSV.csv", skip=11, header=FALSE)
Age_Gender_2012 = read.csv("Gender - Age - 2012 CSV.csv", skip=11, header=FALSE)
Age_Gender_2013 = read.csv("Gender - Age - 2013 CSV.csv", skip=11, header=FALSE)
Age_Gender_2014 = read.csv("Gender - Age - 2014 CSV.csv", skip=11, header=FALSE)
Age_Gender_2015 = read.csv("Gender - Age - 2015 CSV.csv", skip=11, header=FALSE)
Age_Gender_2021 = read.csv("Gender - Age CSV - 2021.csv", skip=11, header=FALSE)

###############################################################################################
##########     Clean Census Data
###############################################################################################

Census_Data_C = Census_Data_C[-(122:124), -10]
colnames(Census_Data_C) = c("Age","Male1996","Female1996","Male2001","Female2001","Male2011","Female2011","MaleTotal","FemaleTotal")

Age_Gender_Population_2022 = Age_Gender_Population_2022[-(122:134), -5]
colnames(Age_Gender_Population_2022) = c("Age","Male2022","Female2022","Total2022")

Ages = 0:120

###############################################################################################
##########     Population Spline Interpolation
###############################################################################################

CensusYears = c(1996,2001,2011,2022)
TargetYears = 2006:2040

MalePopMatrix = cbind(Census_Data_C$Male1996,Census_Data_C$Male2001,Census_Data_C$Male2011,Age_Gender_Population_2022$Male2022)
FemalePopMatrix = cbind(Census_Data_C$Female1996,Census_Data_C$Female2001,Census_Data_C$Female2011,Age_Gender_Population_2022$Female2022)

# Replace zeros with small number only where necessary for log
MalePopMatrix[MalePopMatrix < 1] = 1
FemalePopMatrix[FemalePopMatrix < 1] = 1

MalePopInterp = matrix(NA,length(Ages),length(TargetYears))
FemalePopInterp = matrix(NA,length(Ages),length(TargetYears))

for(i in 1:length(Ages)) {
  model = lm(log(MalePopMatrix[i,]) ~ ns(CensusYears,df=3))
  MalePopInterp[i,] = exp(predict(model,newdata=data.frame(CensusYears=TargetYears)))
}
for(i in 1:length(Ages)) {
  model = lm(log(FemalePopMatrix[i,]) ~ ns(CensusYears,df=3))
  FemalePopInterp[i,] = exp(predict(model,newdata=data.frame(CensusYears=TargetYears)))
}

###############################################################################################
##########     Clean Death Data
###############################################################################################

clean_deaths = function(df) {
  df = df[-(122:134), -5]
  colnames(df) = c("Age","Male","Female","Total")
  df$Age = as.numeric(df$Age)
  Full = data.frame(Age=Ages)
  df = merge(Full,df,by="Age",all.x=TRUE)
  df$Male[is.na(df$Male)] = 0
  df$Female[is.na(df$Female)] = 0
  return(df)
}

DeathDataList = list(Age_Gender_2006, Age_Gender_2008, Age_Gender_2009, Age_Gender_2010,
                     Age_Gender_2011, Age_Gender_2012, Age_Gender_2013, Age_Gender_2014,
                     Age_Gender_2015, Age_Gender_2021)

DeathDataList = lapply(DeathDataList, clean_deaths)

MaleDeathMatrix = do.call(cbind, lapply(DeathDataList, function(x) x$Male))
FemaleDeathMatrix = do.call(cbind, lapply(DeathDataList, function(x) x$Female))

DeathYears = c(2006,2008,2009,2010,2011,2012,2013,2014,2015,2021)

###############################################################################################
##########     Lee-Carter Fit with Stabilization
###############################################################################################

lee_carter_fit = function(D,E) {
  mx = D/E
  mx[mx<=0] = NA  # only treat impossible/missing as NA
  logmx = log(mx)
  ax = rowMeans(logmx, na.rm=TRUE)
  centered = sweep(logmx,1,ax,"-")
  centered[!is.finite(centered)] = 0   # replace remaining NA/Inf with 0 for SVD
  sv = svd(centered)
  bx = sv$u[,1]
  kt = sv$d[1]*sv$v[,1]
  # normalize bx/kt
  bx = bx/sum(bx)
  kt = kt*sum(bx)
  return(list(ax=ax,bx=bx,kt=kt))
}

MaleExposure = MalePopInterp[,1:length(DeathYears)]
FemaleExposure = FemalePopInterp[,1:length(DeathYears)]

MaleLC = lee_carter_fit(MaleDeathMatrix,MaleExposure)
FemaleLC = lee_carter_fit(FemaleDeathMatrix,FemaleExposure)

###############################################################################################
##########     kt Forecast 2022-2040 with Linear + High Age Gompertz Smoothing
###############################################################################################

forecast_kt = function(kt,h) {
  t = 1:length(kt)
  model = lm(kt~t)
  future_t = (length(kt)+1):(length(kt)+h)
  forecast = predict(model,newdata=data.frame(t=future_t))
  return(forecast)
}

ForecastYears = 2022:2040
h = length(ForecastYears)

Male_kt_forecast = forecast_kt(MaleLC$kt,h)
Female_kt_forecast = forecast_kt(FemaleLC$kt,h)

# Construct mortality for each year
construct_mx = function(ax,bx,kt_forecast) {
  mx_forecast = sapply(1:length(kt_forecast), function(i) exp(ax + bx*kt_forecast[i]))
  # Gompertz smoothing for ages 60+
  ages_smooth = 60:length(ax)
  for(i in 1:ncol(mx_forecast)) {
    old_rates = mx_forecast[ages_smooth,i]
    # Fit simple Gompertz: log(mx) ~ age
    fit = lm(log(old_rates) ~ ages_smooth)
    mx_forecast[ages_smooth,i] = exp(predict(fit,newdata=data.frame(ages_smooth=ages_smooth)))
  }
  mx_forecast[mx_forecast>1] = 1  # cap at 1
  return(mx_forecast)
}

Male_mx_forecast = construct_mx(MaleLC$ax, MaleLC$bx, Male_kt_forecast)
Female_mx_forecast = construct_mx(FemaleLC$ax, FemaleLC$bx, Female_kt_forecast)

###############################################################################################
##########     Life Table Function
###############################################################################################

build_life_table = function(mx) {
  qx = mx/(1+0.5*mx)
  qx[qx>1] = 1
  lx = matrix(NA,nrow=nrow(mx),ncol=ncol(mx))
  dx = lx
  Lx = lx
  Tx = lx
  ex = lx
  for(j in 1:ncol(mx)) {
    lx[,j] = NA; lx[1,j] = 100000
    for(i in 2:nrow(mx)) lx[i,j] = lx[i-1,j]*(1-qx[i,j])
    dx[,j] = lx[,j]*qx[,j]
    Lx[,j] = (lx[,j]+c(lx[-1,j],0))/2
    Lx[nrow(mx),j] = lx[nrow(mx),j]/mx[nrow(mx),j]
    Tx[,j] = rev(cumsum(rev(Lx[,j])))
    ex[,j] = Tx[,j]/lx[,j]
  }
  df_list = lapply(1:ncol(mx), function(j) {
    data.frame(Age=Ages, mx=mx[,j], qx=qx[,j], lx=lx[,j], dx=dx[,j], Lx=Lx[,j], Tx=Tx[,j], ex=ex[,j])
  })
  names(df_list) = ForecastYears
  return(df_list)
}

MaleLifeTables = build_life_table(Male_mx_forecast)
FemaleLifeTables = build_life_table(Female_mx_forecast)

###############################################################################################
##########     Export Life Tables
###############################################################################################

setwd("C:\\Users\\mihir\\Documents\\Mihir Documents\\Projects\\Mortality Excel\\Trial 2 Output")

for(yr in ForecastYears) {
  write.csv(MaleLifeTables[[as.character(yr)]], paste0("Male_Life_Table_",yr,".csv"), row.names=FALSE)
  write.csv(FemaleLifeTables[[as.character(yr)]], paste0("Female_Life_Table_",yr,".csv"), row.names=FALSE)
}

###############################################################################################
##########     Plot Survival Curves for 2040
###############################################################################################

plot(Ages,MaleLifeTables[["2040"]]$lx/100000,type="l",col="blue",lwd=2,
     main="Survival Curve Forecast 2040",
     xlab="Age",ylab="Survival Probability")
lines(Ages,FemaleLifeTables[["2040"]]$lx/100000,col="red",lwd=2)
legend("topright",legend=c("Male","Female"),col=c("blue","red"),lwd=2)

###############################################################################################
##########     Plot All Survival Curves 2022-2040
###############################################################################################

# Set up colors: gradient from blue (2022) to red (2040)
library(grDevices)
num_years = length(ForecastYears)
colors = colorRampPalette(c("blue","red"))(num_years)

plot(Ages, MaleLifeTables[[1]]$lx/100000, type="n", ylim=c(0,1),
     xlab="Age", ylab="Survival Probability",
     main="Survival Curves 2022–2040")

# Plot male curves
for(i in 1:num_years) {
  lines(Ages, MaleLifeTables[[i]]$lx/100000, col=colors[i], lwd=2)
}

# Plot female curves as dashed lines
for(i in 1:num_years) {
  lines(Ages, FemaleLifeTables[[i]]$lx/100000, col=colors[i], lwd=2, lty=2)
}

legend("topright", 
       legend=c("Male","Female"),
       col=c("black","black"),
       lty=c(1,2),
       lwd=2,
       bty="n")

# Optional: add color legend to show years

#image.plot(legend.only=TRUE, zlim=range(ForecastYears),
#           col=colors, legend.lab="Year", legend.line=2.5)

###############################################################################################
##########     Create Average Life Table 2022-2040
###############################################################################################

# Initialize matrices to store values
mx_mat = qx_mat = lx_mat = dx_mat = Lx_mat = Tx_mat = ex_mat = matrix(NA, nrow=length(Ages), ncol=length(ForecastYears))

# Fill matrices with values from all tables
for(i in 1:length(ForecastYears)) {
  mx_mat[,i] = MaleLifeTables[[i]]$mx
  qx_mat[,i] = MaleLifeTables[[i]]$qx
  lx_mat[,i] = MaleLifeTables[[i]]$lx
  dx_mat[,i] = MaleLifeTables[[i]]$dx
  Lx_mat[,i] = MaleLifeTables[[i]]$Lx
  Tx_mat[,i] = MaleLifeTables[[i]]$Tx
  ex_mat[,i] = MaleLifeTables[[i]]$ex
}

# Compute row-wise averages across all years
AvgMaleTable = data.frame(
  Age = Ages,
  mx = rowMeans(mx_mat),
  qx = rowMeans(qx_mat),
  lx = rowMeans(lx_mat),
  dx = rowMeans(dx_mat),
  Lx = rowMeans(Lx_mat),
  Tx = rowMeans(Tx_mat),
  ex = rowMeans(ex_mat)
)

# Repeat for females
mx_mat_f = qx_mat_f = lx_mat_f = dx_mat_f = Lx_mat_f = Tx_mat_f = ex_mat_f = matrix(NA, nrow=length(Ages), ncol=length(ForecastYears))

for(i in 1:length(ForecastYears)) {
  mx_mat_f[,i] = FemaleLifeTables[[i]]$mx
  qx_mat_f[,i] = FemaleLifeTables[[i]]$qx
  lx_mat_f[,i] = FemaleLifeTables[[i]]$lx
  dx_mat_f[,i] = FemaleLifeTables[[i]]$dx
  Lx_mat_f[,i] = FemaleLifeTables[[i]]$Lx
  Tx_mat_f[,i] = FemaleLifeTables[[i]]$Tx
  ex_mat_f[,i] = FemaleLifeTables[[i]]$ex
}

AvgFemaleTable = data.frame(
  Age = Ages,
  mx = rowMeans(mx_mat_f),
  qx = rowMeans(qx_mat_f),
  lx = rowMeans(lx_mat_f),
  dx = rowMeans(dx_mat_f),
  Lx = rowMeans(Lx_mat_f),
  Tx = rowMeans(Tx_mat_f),
  ex = rowMeans(ex_mat_f)
)

###############################################################################################
##########     Save Average Life Tables (Optional)
###############################################################################################

write.csv(AvgMaleTable, "Avg_Male_LifeTable_2022_2040.csv", row.names = FALSE)
write.csv(AvgFemaleTable, "Avg_Female_LifeTable_2022_2040.csv", row.names = FALSE)

cat("Average Male Life Table saved to 'Avg_Male_LifeTable_2022_2040.csv'\n")
cat("Average Female Life Table saved to 'Avg_Female_LifeTable_2022_2040.csv'\n")


###############################################################################################
########################################   END   ##############################################
###############################################################################################
