library(foreign)
library(raster)
library(readr)
library(sp)
library(rgdal)
library(scales) #for transparency
library(mgcv)
library(ggplot2)

rm(list = ls())

q = 0.08

coast = readOGR(dsn='/Users/ktanaka/COG/GIS/Gulf of Maine 15 m contour lines/', layer='gom15ctr')
zones = readOGR(dsn='/Users/ktanaka/COG/GIS/Maine Lobster Management Zones/', layer='lob_zone_lines_dd')

a = list("sp.polygons", coast)

load("/Users/ktanaka/COG/Data/fall_spring.RData")
load("/Users/ktanaka/COG/Data/TweedieGAMResult_w_FVCOM.RData")
load("/Users/ktanaka/COG/Data/TweedieGAMResult.RData")

lobster1 = "std_female_adult"
lobster2 = "std_male_adult"
lobster3 = "std_female_juvenile"
lobster4 = "std_male_juvenile"

data = spring

l1 = function(lobster){
  
  # lobster = "std_male_juvenile"
  
  ts2001 = data[ which(data$year == 2001),]
  ts2002 = data[ which(data$year == 2002),]
  ts2003 = data[ which(data$year == 2003),]
  ts2004 = data[ which(data$year == 2004),]
  ts2005 = data[ which(data$year == 2005),]
  ts2006 = data[ which(data$year == 2006),]
  ts2007 = data[ which(data$year == 2007),]
  ts2008 = data[ which(data$year == 2008),]
  ts2009 = data[ which(data$year == 2009),]
  ts2010 = data[ which(data$year == 2010),]
  ts2011 = data[ which(data$year == 2011),]
  ts2012 = data[ which(data$year == 2012),]
  ts2013 = data[ which(data$year == 2013),]
  ts2014 = data[ which(data$year == 2014),]
  
  x = list(ts2001=ts2001, ts2002=ts2002, ts2003=ts2003, ts2004=ts2004, 
           ts2005=ts2005, ts2006=ts2006, ts2007=ts2007, ts2008=ts2008, 
           ts2009=ts2009, ts2010=ts2010, ts2011=ts2011, ts2012=ts2012,
           ts2013=ts2013, ts2014=ts2014)
  
  gravity <- matrix(0, length(x), ncol = 4)
  
  for (i in 1:length(x)){
    
    x[[i]] = na.exclude(x[[i]][, c("La", "Lo", "D", lobster)])
    tstemp = data.frame(x[[i]])
    
    X=sum(tstemp[, 4]*tstemp[, 2])/sum(tstemp[, 4]) #for Lon
    Y=sum(tstemp[, 4]*tstemp[, 1])/sum(tstemp[, 4]) #for Lat
    Z=sum(tstemp[, 4]*tstemp[, 3])/sum(tstemp[, 4]) #for Depth
    
    gravity[i, 1] = X
    gravity[i, 2] = Y
    gravity[i, 3] = Z
    
  }
  
  gravity[,4] = c(2001:2014)
  
  gt = as.data.frame(gravity)
  
  list1 <- cbind(1:5, rep("2001-2005"))
  list2 <- cbind(6:10, rep("2006-2010"))
  list3 <- cbind(11:14, rep("2011-2014"))
  list = rbind(list1, list2, list3)
  
  gt = cbind(gt, list[,2])
  
  colnames(gt) = c("Lon", "Lat", "Depth", "Year", "Period")
  
  st.err <- function(x) {
    sd(x)/sqrt(length(x))
  }
  
  lon_se <- aggregate(Lon ~ Period, gt, st.err)
  lat_se <- aggregate(Lat ~ Period, gt, st.err)
  se = merge(lon_se, lat_se)
  colnames(se) = c("Period","Lon_se","Lat_se")
  
  rm(list, list1, list2, list3, gt, lon_se, lat_se)
  
  y1 = gravity[c(1:5),] #2001-2005
  y2 = gravity[c(6:10),] #2006-2010
  y3 = gravity[c(11:14),] #2011-2014
  
  n = 5
  z1 = aggregate(y1, list(rep(1:(nrow(y1)%/%n+1), each = n, len = nrow(y1))), mean, na.rm=TRUE)[-1]
  z1[,4] = c("2001-2005")
  
  z2 = aggregate(y2, list(rep(1:(nrow(y2)%/%n+1), each = n, len = nrow(y2))), mean, na.rm=TRUE)[-1]
  z2[,4] = c("2006-2010")
  
  n = 4
  z3 = aggregate(y3, list(rep(1:(nrow(y3)%/%n+1), each = n, len = nrow(y3))), mean, na.rm=TRUE)[-1]
  z3[,4] = c("2011-2014")
  
  gravity = rbind(z1, z2, z3)
  colnames(gravity) = c("Lon","Lat","Depth","Period")
  
  n = 5
  x1 = aggregate(y1, list(rep(1:(nrow(y1)%/%n+1), each = n, len = nrow(y1))), sd, na.rm=TRUE)[-1]
  x1[,4] = c("2001-2005")
  
  x2 = aggregate(y2, list(rep(1:(nrow(y2)%/%n+1), each = n, len = nrow(y2))), sd, na.rm=TRUE)[-1]
  x2[,4] = c("2006-2010")
  
  n = 4
  x3 = aggregate(y3, list(rep(1:(nrow(y3)%/%n+1), each = n, len = nrow(y3))), sd, na.rm=TRUE)[-1]
  x3[,4] = c("2011-2014")
  
  gravity_sd = rbind(x1,x2,x3)
  colnames(gravity_sd) = c("Lon_sd","Lat_sd","Depth_sd","Period")
  
  gravity = cbind(gravity[,c(1:2)], gravity_sd[,c(1:2)], se[,c(2:3,1)])
  
  points(gravity[, 1], gravity[, 2], pch = 16, col = 4, cex = 1)
  #text(gravity[,1], gravity[,2], gravity[,4], cex = 2, pos = 1, offset = 2.5)
  
  s <- seq(length(x)-1)  # one shorter than data
  arrows(gravity[s,1], gravity[s,2], gravity[s+1,1], gravity[s+1,2], 
         length = 0.2, angle=20, lwd = 2, col = alpha("blue", 1))
  scale = gravity
  coordinates(scale)<-~Lon+Lat
  scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))
  
  return(gravity)
  
  
}#observed COGs
l2 = function(lobster){
  
  ts2001 = data[ which(data$year == 2001),]
  ts2002 = data[ which(data$year == 2002),]
  ts2003 = data[ which(data$year == 2003),]
  ts2004 = data[ which(data$year == 2004),]
  ts2005 = data[ which(data$year == 2005),]
  ts2006 = data[ which(data$year == 2006),]
  ts2007 = data[ which(data$year == 2007),]
  ts2008 = data[ which(data$year == 2008),]
  ts2009 = data[ which(data$year == 2009),]
  ts2010 = data[ which(data$year == 2010),]
  ts2011 = data[ which(data$year == 2011),]
  ts2012 = data[ which(data$year == 2012),]
  ts2013 = data[ which(data$year == 2013),]
  ts2014 = data[ which(data$year == 2014),]
  
  x = list(ts2001=ts2001, ts2002=ts2002, ts2003=ts2003, ts2004=ts2004, 
           ts2005=ts2005, ts2006=ts2006, ts2007=ts2007, ts2008=ts2008, 
           ts2009=ts2009, ts2010=ts2010, ts2011=ts2011, ts2012=ts2012,
           ts2013=ts2013, ts2014=ts2014)
  
  gravity <- matrix(0, length(x), ncol = 4)
  
  if (lobster == "std_female_adult") {
    model = sp_f_adu
  }
  if (lobster == "std_male_adult") {
    model = sp_m_adu
  }
  if (lobster == "std_female_juvenile") {
    model = sp_f_juv
  }
  if (lobster == "std_male_juvenile") {
    model = sp_m_juv
  }
  
  for (i in 1:length(x)){
    
    x[[i]] = na.exclude(x[[i]][, c("La", "Lo", "D","Te", "S", "Do", lobster)])
    tstemp = data.frame(x[[i]])
    tstemp$prediction = predict.gam(model, tstemp, se.fit = F, type = "response")
    
    
    X=sum(tstemp[, 8]*tstemp[, 2])/sum(tstemp[, 8]) #for Lon
    Y=sum(tstemp[, 8]*tstemp[, 1])/sum(tstemp[, 8]) #for Lat
    Z=sum(tstemp[, 8]*tstemp[, 3])/sum(tstemp[, 8]) #for Depth
    
    gravity[i, 1] = X
    gravity[i, 2] = Y
    gravity[i, 3] = Z
    
  }
  
  
  gravity[,4] = c(2001:2014)
  
  y1 = gravity[c(1:5),] #2001-2005
  y2 = gravity[c(6:10),] #2006-2010
  y3 = gravity[c(11:14),] #2011-2014
  
  
  n = 5
  y1 = aggregate(y1, list(rep(1:(nrow(y1)%/%n+1), each = n, len = nrow(y1))), mean, na.rm=TRUE)[-1]
  y1[,4] = c("2001-2005")
  
  y2 = aggregate(y2, list(rep(1:(nrow(y2)%/%n+1), each = n, len = nrow(y2))), mean, na.rm=TRUE)[-1]
  y2[,4] = c("2006-2010")
  
  n = 4
  y3 = aggregate(y3, list(rep(1:(nrow(y3)%/%n+1), each = n, len = nrow(y3))), mean, na.rm=TRUE)[-1]
  y3[,4] = c("2011-2014")
  
  gravity = rbind(y1, y2, y3)
  
  points(gravity[, 1], gravity[, 2], pch = 16, col = 2, cex = 1)
  #text(gravity[,1], gravity[,2], gravity[,4], cex = 2, pos = 1, offset = 2.5)
  
  s <- seq(length(x)-1)  # one shorter than data
  arrows(gravity[s,1], gravity[s,2], gravity[s+1, 1], gravity[s+1, 2], length = 0.2, angle=20, lwd = 2, col = alpha("red", 1))
  coordinates(gravity)<-~V1+V2
  scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))
}#modeled COGs

xlim = c(-69.4, -68.7)
ylim = c(43.8, 44.15)

plot(coast,
     cex.main = 2,
     xlim = xlim,
     ylim = ylim,
     xlab = "", ylab = "",
     col = "gray", xaxt = "n", yaxt = "n", lwd = 0.5, cex.lab = 1)
lines(zones)
maps::map(add = T, fill = T)

degAxis(1, cex.axis = 1)
degAxis(2, las = 2, cex.axis = 1)
scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))

rm(mean)

l1(lobster1)#FA
legend("bottomright", "SP_ADU_F", lty = 0.5, lwd = 0.5, bty = "n", cex = 4)
l2(lobster1)#FA

plot(NULL,
     cex.main = 2,
     xlim = xlim,
     ylim = ylim,
     xlab = "", ylab = "",
     col = "gray", xaxt = "n", yaxt = "n", lwd = 0.5, cex.lab = 1)
lines(zones)
maps::map(add = T, fill = T)

degAxis(1, cex.axis = 1)
degAxis(2, las = 2, cex.axis = 1)
scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))

l1(lobster2)#MA
legend("bottomright", "SP_ADU_M", lty = 0.5, lwd = 0.5, bty = "n", cex = 4)
l2(lobster2)#MA

plot(NULL,
     cex.main = 2,
     xlim = xlim,
     ylim = ylim,
     xlab = "", ylab = "",
     col = "gray", xaxt = "n", yaxt = "n", lwd = 0.5, cex.lab = 1)
lines(zones)
maps::map(add = T, fill = T)

degAxis(1, cex.axis = 1)
degAxis(2, las = 2, cex.axis = 1)
scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))

l1(lobster3)#FJ
legend("bottomright", "SP_JUV_F", lty = 0.5, lwd = 0.5, bty = "n", cex = 4)
l2(lobster3)#FJ

plot(NULL,
     cex.main = 2,
     xlim = xlim,
     ylim = ylim,
     xlab = "", ylab = "",
     col = "gray", xaxt = "n", yaxt = "n", lwd = 0.5, cex.lab = 1)
lines(zones)
maps::map(add = T, fill = T)

degAxis(1, cex.axis = 1)
degAxis(2, las = 2, cex.axis = 1)
scalebar(d = 20, xy = NULL, type = 'bar', divs = 4, below = 'km', cex = 2, adj = c(0,-1))

l1(lobster4)#MJ
legend("bottomright", "SP_JUV_M", lty = 0.5, lwd = 0.5, bty = "n", cex = 4)
l2(lobster4)#MJ
