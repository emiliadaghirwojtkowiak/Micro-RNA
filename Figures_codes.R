# miRNA project ##
# Figures #


#########
#Fig.1 ##
#########

# Raw data (centered and standardized miRNA levels) summarized as boxplots for 49 miRNAs 
# in patients and controls. The box and whiskers plots depict mean, 
# 25th and 75th percentiles. Blue dots overlaid are individual data points. 

rm(list=ls())

data1 = read.csv("MyDataEDW.csv")
colnames(data1) <- c("SUBID", "MIRNAID", "casecont","DV", "hist", "age" , "weight","group" , "study_id", "Ycensord", "REP"); head(data1);
data1 <-as.data.frame(data1)
data <- subset(data1, data1$study_id =="1") # RT1 
data <- subset(data, data$REP =="1" & data$MIRNAID<50) # RT1 i bez kontrolnego miRNA
data$DV =  data$DV - 40*log(2)
dim(data)


font.settings <- list( font = 0.7, cex = 0.7)
my.theme <- list(
  box.umbrella = list(col="salmon",alpha=0.4),
  box.rectangle = list(col="salmon",fill="salmon",alpha=0.4),
  #box.dot = list(col = "black", cex=.5),
  plot.symbol = list(cex = .5, col = 1), #outlier size and color
  par.xlab.text = font.settings,
  par.ylab.text = font.settings,
  axis.text = font.settings,
  par.sub=font.settings)

bwplot(scale(DV) ~ as.factor(casecont) | mirna.names, data=data, layout=c(10,5),par.settings = my.theme,
       par.strip.text=list(cex=0.7), cex=0.5,lty=1.3,alpha=.3,do.out = F,ylab="Centered and standardized miRNA",
       panel=function(x,y,...){
         panel.bwplot(x,y,...)
         panel.stripplot(x,y,col="dodgerblue",do.out=T,jitter.data=TRUE,factor=0.5,...)
       })

###########
##Fig. 2###
###########
# The summary of a marginal posterior distribution representing fold change between disease and control subjects for 49 miRNA. The distribution was summarized as a boxplot with 5th, 25th, 50th, 75th and 95th percentile. 
# Grey line denotes no effect for miRNAs. 


boxplot((lapply(as.data.frame(exp(betaCANCER)), quantile, probs = c(0.95, 0.5, 0.05))), outline = FALSE,names = mirna.names, las=2,cex.axis=0.45, medcol = "red", whisklty = 1,staplelwd = 4, outpch = 8,
        cex.lab=1.75, mgp=c(1,1,0),mgp= c(3, 0, 0))
title(xlab="", ylab = expression(paste("(",beta[j], ")"), c(1,2,2,3) ,cex.lab=3,xaxt="n"),mgp= c(2, 1,1))
abline(h=1 ,col = "darkgray", lwd=2)
abline(h = 1.5, col = "darkred", lwd=2) 

#########
#Fig.3###
#########

# The summary of a marginal posterior distribution of an effect size for 49 miRNA. The distribution was summarized as a boxplot with 5th, 25th, 50th, 75th and 95th percentile. 
# Grey line describes no effect. 

EffSize1 <- select(samples, starts_with("EffSize1["));dim(EffSize1);
EffSize2 <- select(samples, starts_with("EffSize2["));dim(EffSize2);

setwd("C:/dane/prof_Limon/rak_jajnika/codes/model with covariares/models_3/Results_Model006")
tiff("effSize_M006.tiff", width = 10, height = 8, units = 'in', res = 300)


boxplot((lapply(as.data.frame((EffSize2)), quantile, probs = c(0.95, 0.5, 0.05))),lty=1, medcol = "white", boxlty = 0, whisklty = 4, staplelwd = 4, outpch = 8, outcex = 3,
        ylim=c(-4,4), at = 1:ncol(EffSize1) + 0.5, boxfill="red", boxwex=0.25, staplewex = 1, outwex = 1, boxfill = "grey34",pch=18, medlty = 1,medlwd = 2.5, lwd=.5,cex.lab=1.75, mgp=c(1,1,0),mgp= c(3, 0, 0),
        xaxt="n")

boxplot((lapply(as.data.frame((EffSize1)), quantile, probs = c(0.95, 0.5, 0.05))),lty=1, medcol = "white", boxlty = 0, whisklty = 4, staplelwd = 4, outpch = 8, outcex = 3,
        add = TRUE,boxwex = 0.25, staplewex = 1, outwex = 1, boxfill = "gray20",pch=18, medlty = 1,medlwd = 2.5, lwd=.5,cex.lab=1.75, mgp=c(1,1,0),mgp= c(3, 0, 0),
        names= mirna.names[-50],las=3)
abline(h= 0 ,col = "darkgray", lwd=2);
abline(h= 0.138 ,col = "blue", lwd=2);



# Fig.4 #
#########

layout(matrix(c(1,2,3,4,5,5), 2,1, byrow = TRUE), widths=c(1,1), heights=c(1,1))

## OMEGA
c <- data.frame(omega3, omega1)
boxplot((lapply(as.data.frame((omega3)), quantile, probs = c(0.95, 0.5, 0.05))), outline = FALSE,las=2,cex.axis=0.45,medcol = "red", whisklty = 1,staplelwd = 4, outpch = 8,
        cex.lab=1.75, mgp=c(1,1,0),xaxt="n",
        col=names(c)=="omega1",rgb(0.8,0.1,0.3,0.6))
title(ylab = expression(paste(" ",omega[j],""),c(5, 8, 9, 5) + 0.1 ,cex.lab=3),mgp= c(2, 1,1))
abline(h = 0, col = "darkgray", lwd=2)

## SIGMA
boxplot((lapply(as.data.frame((sigma)), quantile, probs = c(0.95, 0.5, 0.05))),names = mirna.names, outline = FALSE,las=2,cex.axis=0.45, medcol = "red", whisklty = 1,staplelwd = 4, outpch = 8,
        cex.lab=1.75, mgp=c(1,1,0),mgp= c(3, 0, 0))
title(xlab="", ylab = expression(paste("",sigma[j], ""), c(1,2,2,3) ,cex.lab=3),mgp= c(2, 1,1))
abline(h= 0 ,col = "darkgray", lwd=2)


########
#Fig.5#
#######

#Weighted residuals (top) and weighted residuals versus fitted values (bottom) for 49 miRNAs.
#At the top figure the dots are distributed across the line of identity. 
#At the bottom figure the dots are quite evenly distributed across zero with no visible pattern or trend.


#predicted values
Ypred.1 <- select(samples, starts_with("Ycond"));dim(Ypred.1)
d11 = as.data.frame(Ypred.1); dim(d11)
tr11 <- t(d11); dim(tr11)
xsplit <- rep(1:49)
tmp11 <- split.data.frame(tr11,xsplit); class(tmp11)
str(tmp11)
library("abind")
all.matrix.RT11 <- abind(tmp11, along=3); dim(all.matrix.RT11) # 183 / 3000/ 50, Combine multi-dimensional arrays
xx.RT11 <- apply(all.matrix.RT11, c(1,2), mean); dim(xx.RT11) # 183 / 1002 , c(1,2) means keep the 1st and 2nd dimensions and apply the function (mean) over the 3rd
xx11 <- apply(all.matrix.RT11, 1, colMeans); dim(xx11) #49/161

#observed
setwd("C:/dane/prof_Limon/rak_jajnika/codes/model with covariares/data_file/")
data <- read.csv("MyDataEDW.csv")
colnames(data) <- c("SUBID", "MIRNAID", "casecont","DV", "hist", "age" , "weight","group" , "study_id", "Ycensord", "REP"); head(data)
data <- subset(data, data$study_id =="1" & data$REP=="1" & data$MIRNAID<50)


Yc = (40-CT)/log(2);
Y = (Yc - mean(Yc, na.rm=T))/ sd(Yc, na.rm=T);
Y <- t(matrix(Y,nrow = 178, ncol = 49))
predicted <- as.numeric(xx11);  
observed <- as.matrix(Y); 
observed<- as.numeric(Y); length(observed);

plot(predicted , observed, col=topo.colors(222), xlab="Predicted CT", ylab="Observed CT",
     cex = 1.3, xlim=c(-4,4), ylim=c(-4,4) ); segments(-50,-50,50,50,col = "darkgray", lwd=2);

# residuals
sigma <- select(samples, starts_with("sigma"));dim(sigma)
colnames(sigma) <- NULL; rownames(sigma) <- NULL; 

# WRES gruped by MIRNAID
layout(matrix(c(1,2,3,4,5,5), 2,1, byrow = TRUE), widths=c(1,1), heights=c(1,1));


a <- xyplot(observed-predicted/length(sigma[MIRNAID]) ~ MIRNAID, ylab="WRES",xlab="miRNA",grid = TRUE,
            scales = list(x = list(at = "", equispaced.log = FALSE)),
            type = c("p", "smooth"), lwd = 2,col.line = "darkorange")


b <- xyplot(observed - predicted ~ predicted[ID]/length(sigma[MIRNAID]),grid = TRUE, 
            ylab="WRES",xlab="predicted",do.out=F,
            scales = list(x = list(at = c(-2,0,2), equispaced.log = T)),
            type = c("p", "smooth"), lwd = 2,col.line = "darkorange")


grid.arrange(a,b, nrow=2)




