library(rpart)
library(ggplot2)
library(rpart.plot)
library(plyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
XtremLow <- function(Tcl, Lat, Lon, Elev){
  pacificsouth <- 1/((((Lat - -22.7)/13)^2 + ((Lon - -82.3)/14)^2)^2+1)
  amazon2 <- 1/((((Lat - -10.2)/5)^2 + ((Lon - -59.9)/10)^2)^2+1)
  amazon1 <- 1/((((Lat - -2.8)/14)^2 + ((Lon - -61.3)/19)^2)^2+1)
  pacificcent <- 1/((((Lat - 4.1)/21)^2 + ((Lon - -122.4)/41)^2)^2+1)
  mexico <- 1/((((Lat - 26)/6)^2 + ((Lon - -98.4)/12)^2)^2+1)
  florida <- 1/((((Lat - 27.5)/4)^2 + ((Lon - -81.1)/8)^2)^2+1)
  pacificnorth <- 1/((((Lat - 32.9)/26)^2 + ((Lon - -145)/27)^2)^2+1)
  oklahoma <- 1/((((Lat - 33.6)/4)^2 + ((Lon - -98.4)/8)^2)^2+1)
  arizona <- 1/((((Lat - 34)/12)^2 + ((Lon - -113.1)/8)^2)^2+1)
  atlantic <- 1/((((Lat - 34)/15)^2 + ((Lon - -60.7)/19)^2)^2+1)
  himalayas <- 1/((((Lat - 35.3)/6)^2 + ((Lon - 91.3)/13)^2)^2+1)
  kentucky <- 1/((((Lat - 38.5)/3)^2 + ((Lon - -87.6)/9)^2)^2+1)
  detroit <- 1/((((Lat - 41.8)/3)^2 + ((Lon - -82.6)/4)^2)^2+1)
  ontario <- 1/((((Lat - 44.6)/2)^2 + ((Lon - -79.2)/6)^2)^2+1)
  montana <- 1/((((Lat - 45.4)/5)^2 + ((Lon - -111.8)/10)^2)^2+1)
  minn <- 1/((((Lat - 47.6)/6)^2 + ((Lon - -92.6)/12)^2)^2+1)
  hudson <- 1/((((Lat - 60)/7)^2 + ((Lon - -87)/34)^2)^2+1)
  siberia <- 1/((((Lat - 61.2)/20)^2 + ((Lon - 105.7)/39)^2)^2+1)
  california <- 1/((((Lat - 34.8)/9)^2 + ((Lon - -128.2)/9)^2)^2+1)
  washington <- 1/((((Lat - 46)/5)^2 + ((Lon - -126.6)/5)^2)^2+1)
  colorado <- 1/((((Lat - 38.3)/2)^2 + ((Lon - -108.8)/3)^2)^2+1)
  hawaii <- 1/((((Lat - 21.3)/7)^2 + ((Lon - -157.5)/11)^2)^2+1)
  chess <- 1/((((Lat - 37)/3)^2 + ((Lon - -74)/3)^2)^2+1)
  
  Tclx<-	-9.171	+
    Tcl *	1.202	+
    Lat *	-0.04149	+
    Elev *	0.0008691	+
    Lat * Elev *	-0.00002455	+
    pacificsouth *	-1.792	+
    amazon2 *	2.573	+
    amazon1 *	-1.014	+
    pacificcent *	-0.749	+
    mexico *	-0.8227	+
    florida *	-3.557	+
    pacificnorth *	-1.246	+
    oklahoma *	0.1758	+
    arizona *	2.605	+
    chess *	0.8347	+
    atlantic *	0.2967	+
    himalayas *	-1.814	+
    kentucky *	-2.644	+
    detroit *	0	+
    ontario *	-2.314	+
    montana *	-4.415	+
    minn *	1.136	+
    hudson *	-5.154	+
    siberia *	-3.797	+
    california *	4.48	+
    washington *	3.597	+
    colorado *	1.458	+
    hawaii *	6.673	
  return(Tclx)}

Days <- c(31.00, 28.25, 31.00, 30.00, 31.00, 30.00, 31.00, 31.00, 30.00, 31.00, 30.00, 31.00)
DayNumber <- c(16.000,45.625,75.250,106.125,136.250,166.750,197.250,228.250,258.750,289.250,319.750,350.250)
dcl <- 0.409*sin(2*3.141592*DayNumber/365-1.39)

GetSolarRad <- function(Month, Lat){
  declination <- dcl[Month]
  
  hs <- acos(pmin(pmax(-tan(Lat/360*2*3.141592) * tan(declination),-1),1))
  Ra <- 117.5 * (hs*sin(Lat/360*2*3.141592)*sin(declination) +
                   cos(Lat/360*2*3.141592)*cos(declination)*sin(hs)) / 3.141592
  return(Ra)
}

GetDayLength<- function(Month, Lat){
  declination <- dcl[Month]
  
  Dl <- ifelse(Lat + declination*360/2/3.141592 > 89.16924, 24, ifelse(Lat - declination*360/2/3.141592 >= 90, 0, (atan(-((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))/(-((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))*((sin(-0.83/360*2*3.141592)-sin(declination)*sin(Lat/360*2*3.141592))/(cos(declination)*cos(Lat/360*2*3.141592)))+1)^0.5)+2*atan(1))/3.141592*24))
  return(Dl)}


GetVp  <- function(p,th,tl) {#Based on linear regression using 10 minute WorldClim 2.0 data with vapor pressure estimates
  Vpmax = 0.6108*exp(17.27*th/(th+237.3)) #saturation vapor pressure kPa
  Vpmin = 0.6108*exp(17.27*tl/(tl+237.3)) #saturation vapor pressure kPa
  Vp0 <- (Vpmin*7.976e-01+
            Vpmin*log(p+1)*9.499e-02+
            Vpmin*Vpmax*-6.599e-02)
  Vp <- pmax(0,pmin(Vpmin,Vp0))
  return(Vp)}

GetSolar <- function(Ra, Elev, th, tl, p) {#Based on linear regression using 10 minute WorldClim 2.0 data with solar radiation estimates
  Rso <- (0.75+2*10^-5*Elev)*Ra
  Rs0 <- (Rso*9.521e-01+
            Rso*log(p+1)*-9.087e-02+
            Rso*tl*-3.644e-03+
            Rso*log(p+1)*th*1.335e-03)
  Rs <- pmax(0.3*Rso,pmin(Rso,Rs0))
  return(Rs)}

GetPET <- function(Ra, th, tl, p){
  Vpmax = 0.6108*exp(17.27*th/(th+237.3)) #saturation vapor pressure kPa
  Vpmin = 0.6108*exp(17.27*tl/(tl+237.3)) #saturation vapor pressure kPa
  logp <- log(p+1)
  e0 <- Ra*0.0508780  +
          Vpmax*0.7893714  +
          Vpmin*-0.5589255  +
          logp*-0.1309403  +
          Ra*Vpmax*0.0049383
  e <- pmax(0,e0)
    return(e)}



GetNetSolar <- function(Ra, Elev, th, tl, p){
  Vp = GetVp(p,th,tl)
  Rso <- (0.75+2*10^-5*Elev)*Ra
  Rs <- GetSolar(Ra, Elev, th, tl, p)
  Rnl <- 4.901*10^-9 * (1.35*Rs/(Rso+0.000001)-0.35) * (0.34 - 0.14 * Vp^0.5) * ((th+273.16)^4 + (tl+273.16)^4)/2
  Rns <- (1-0.23)*Rs
  Rn <- pmax(0,Rns - Rnl)
  return(Rn)}

GetTransGrow <- function(th, tl) {#Adjust to reduction in transpiration due to cold, with evaporation only outside growing season
  ts = 0.8 #assumed T/ET ratio during growing season
  tw = 0.0 #assumed T/ET ratio during freezing season
  t <- (th+tl)/2
  tr <- 10 #generally as mean temperatures get below 10 transpiration shuts down, regardless of warm daytime temperatures
  G0 <- (t-0)/(tr) 
  evmin = (tw)+(1-ts)
  G = G1*(1-evmin)+evmin
  return(G)}

Biomeclimate <- readRDS(file='data/RadBiomeclimate.RDS')


t. = grep("^t01$", colnames(Biomeclimate)):grep("^t12$", colnames(Biomeclimate))
th. = grep("^th01$", colnames(Biomeclimate)):grep("^th12$", colnames(Biomeclimate))
tl. = grep("^tl01$", colnames(Biomeclimate)):grep("^tl12$", colnames(Biomeclimate))
p. = grep("^p01$", colnames(Biomeclimate)):grep("^p12$", colnames(Biomeclimate))
e. = grep("^e01$", colnames(Biomeclimate)):grep("^e12$", colnames(Biomeclimate))

#PET and derivatives
for (i in 1:12){#i=1
  Biomeclimate[,e.[i]] <- 0.85*GetTransGrow(Biomeclimate[,th.[i]], Biomeclimate[,tl.[i]])*
    GetPET(GetSolarRad(i,Biomeclimate$Latitude), Biomeclimate[,th.[i]], Biomeclimate[,tl.[i]], Biomeclimate[,p.[i]])*
    Days[i]
}
Biomeclimate$pAET <- pmin(Biomeclimate[,e.[1]], Biomeclimate[,p.[1]])
for(i in 1:12){
  Biomeclimate$var <- pmin(Biomeclimate[,e.[i]], Biomeclimate[,p.[i]])
  Biomeclimate$pAET <- pmax(Biomeclimate$pAET, Biomeclimate$var)
};Biomeclimate$var = NULL

Biomeclimate$Surplus <- 0
for(i in 1:12){
  Biomeclimate$var <- pmax(0,Biomeclimate[,p.[i]]-Biomeclimate[,e.[i]])
  Biomeclimate$Surplus <- Biomeclimate$Surplus+Biomeclimate$var
};Biomeclimate$var = NULL

Biomeclimate$Deficit <- 0
for(i in 1:12){
  Biomeclimate$var <- pmax(0,Biomeclimate[,e.[i]]-Biomeclimate[,p.[i]])
  Biomeclimate$Deficit <- Biomeclimate$Deficit+Biomeclimate$var
};Biomeclimate$var = NULL

#Tg
Biomeclimate$TgA <- 0
Biomeclimate$TgB <- 0
monthsA <- c(1,2,3,4,11,12)
monthsB <- c(5,6,7,8,9,10)
for(i in 1:6){
  Biomeclimate$varA <- pmax(0, Biomeclimate[,t.[monthsA[i]]])
  Biomeclimate$TgA <- Biomeclimate$TgA+Biomeclimate$varA
  Biomeclimate$varB <- pmax(0, Biomeclimate[,t.[monthsB[i]]])
  Biomeclimate$TgB <- Biomeclimate$TgB+Biomeclimate$varB
}
Biomeclimate$Tg <- pmax(Biomeclimate$TgA,Biomeclimate$TgB)/6
Biomeclimate$varA=NULL;Biomeclimate$TgA=NULL;Biomeclimate$varB=NULL;Biomeclimate$TgB=NULL

#Twh
Biomeclimate$Twh <- apply(Biomeclimate[,th.[1:12]], MARGIN = 1, FUN = 'max') 
#Tw
Biomeclimate$Tw <- apply(Biomeclimate[,t.[1:12]], MARGIN = 1, FUN = 'max') 
#Tc
Biomeclimate$Tc <- apply(Biomeclimate[,t.[1:12]], MARGIN = 1, FUN = 'min') 
#Tcl
Biomeclimate$Tcl <- apply(Biomeclimate[,tl.[1:12]], MARGIN = 1, FUN = 'min') 
#Tclx
Biomeclimate$Tclx <- XtremLow(Biomeclimate$Tcl, Biomeclimate$Latitude, Biomeclimate$Longitude, Biomeclimate$Elevation)

saveRDS(Biomeclimate, 'data/RadBiomeclimate.RDS')


#####
Biomeclimate <- readRDS( 'data/RadBiomeclimate.RDS')



unique(Biomeclimate[,c('BIOME', 'biomname')])
Biomeclimate$PET <- apply(Biomeclimate[,which(colnames(Biomeclimate)=='e01'):which(colnames(Biomeclimate)=='e12')], MARGIN = 1, FUN='sum')

Biomeclimate$MAP <- apply(Biomeclimate[,which(colnames(Biomeclimate)=='p01'):which(colnames(Biomeclimate)=='p12')], MARGIN = 1, FUN='sum')

Biomeclimate$MAAT <- apply(Biomeclimate[,which(colnames(Biomeclimate)=='t01'):which(colnames(Biomeclimate)=='t12')], MARGIN = 1, FUN='mean')

Biomeclimate$M <- (Biomeclimate$MAP/Biomeclimate$PET +0.0001)
Biomeclimate$Cindex <- pmin(Biomeclimate$Tc,Biomeclimate$Tclx+15)
#Biomeclimate$Deficit <- pmin(200, Biomeclimate$Deficit)
#Biomeclimate$Surplus <- pmin(200, Biomeclimate$Surplus)
corep <- function(x){
  e=x[1:12]
  p=x[13:24]
  cor(e,p)
}
Biomeclimate$corep <- apply(Biomeclimate[,c(e.,p.)], 1, corep)
fp3aet <- function(x){
  qmon = c(11,12,1,2,3,4,5,6,7,8,9,10,11,12)
  e=x[1:12]
  p=x[13:24]
  a = c(1:12)*0
  for(i in 1:12){
    a[i] = min(e[qmon[i]],p[qmon[i]])+
      min(e[qmon[i+1]],p[qmon[i+1]])+
      min(e[qmon[i+2]],p[qmon[i+2]])
  }
  p3aet = max(a)
  return(p3aet)
}
Biomeclimate$p3aet <- apply(Biomeclimate[,c(e.,p.)], 1, fp3aet)



selectBiome<-subset(Biomeclimate, !is.na(corep))
selectBiome$biom2 <- 'not'
selectBiome$biom2 <- ifelse(selectBiome$Tg <12, 'boreal',
                            ifelse(selectBiome$Cindex >= 15, 'tropical',
                                   ifelse(selectBiome$Deficit < 150 & selectBiome$M >= 1, 'isopluvial',
                                          ifelse(selectBiome$Surplus < 25 & selectBiome$M < 0.5, 'isoxeric',
                                                 ifelse(selectBiome$corep < 0, 'xerothermic','pluviothermic')))))
# selectBiome$biom2 <- ifelse(selectBiome$BIOME %in% c(1,3,5,6), 'tree',
#                             ifelse(selectBiome$BIOME %in% c(2,4), 'tree',
#                                    ifelse(selectBiome$BIOME %in% c(12), 'shrub',
#                                           ifelse(selectBiome$BIOME %in% c(7,8), 'grass',
#                                                  ifelse(selectBiome$BIOME %in% c(13), 'desert',
#                                                         ifelse(selectBiome$BIOME %in% c(10,11), 'tundra',
#                                                  'not'))))))
# 
# selectBiome <-  subset(selectBiome, !biom2 %in% 'not' & Tg > 12 & Cindex < 15 & Cindex > 0 & biom2 %in% c('grass','shrub'))
selectBiomecount<-aggregate(selectBiome$ECO_NAME, by=list(selectBiome$biom2),FUN=length)
colnames(selectBiomecount)<-c("biom2","x")
selectBiomecount$wt <- 1/(selectBiomecount$x+1)
selectBiome<-merge(selectBiome,selectBiomecount, by="biom2")
colnames(selectBiome)

biomeclass <- rpart(biom2 ~  p3aet + Tg + Cindex +  M + Deficit + Surplus, 
                    data = selectBiome,weights=selectBiome$wt, method="class", 
                    control = list(maxdepth = 5, cp=0.0005, minsplit=100))
# Make plot
biomeclass
png(filename="biome2class.png",width = 10, height = 3, units = 'in', res = 600)

rpart.plot(biomeclass, extra=108) # Make plot

dev.off()




selectBiome<-subset(Biomeclimate, BIOME %in% c(1,2,4,5,6,7,8,10,11,12,13))

selectBiome<-subset(Biomeclimate, BIOME %in% c(7,12,13))
selectBiome <- subset(selectBiome, !is.na(corep))
selectBiomecount<-aggregate(selectBiome$ECO_NAME, by=list(selectBiome$biomname),FUN=length)
colnames(selectBiomecount)<-c("biomname","x")
selectBiomecount$wt <- 1/(selectBiomecount$x+1)
selectBiome<-merge(selectBiome,selectBiomecount, by="biomname")
colnames(selectBiome)
selectBiome$corep
biomeclass <- rpart(biomname ~  p3aet + Tg + Cindex +  M + Deficit + Surplus + pAET +corep, data = selectBiome,weights=selectBiome$wt, method="class", control = list(maxdepth = 2, cp=0.0005, minsplit=100))
# Make plot
biomeclass

biomeclass <- rpart(biomname ~  p3aet + Tg + Cindex +  M + Deficit + Surplus + pAET, data = selectBiome,weights=selectBiome$wt, method="class", control = list(maxdepth = 5, cp=0.0005, minsplit=100))
# Make plot

png(filename="biomeclass.png",width = 10, height = 3, units = 'in', res = 600)

rpart.plot(biomeclass, extra=108) # Make plot

dev.off()



t. = grep("^t01$", colnames(Biomeclimate)):grep("^t12$", colnames(Biomeclimate))
th. = grep("^th01$", colnames(Biomeclimate)):grep("^th12$", colnames(Biomeclimate))
tl. = grep("^tl01$", colnames(Biomeclimate)):grep("^tl12$", colnames(Biomeclimate))
p. = grep("^p01$", colnames(Biomeclimate)):grep("^p12$", colnames(Biomeclimate))
e. = grep("^e01$", colnames(Biomeclimate)):grep("^e12$", colnames(Biomeclimate))





