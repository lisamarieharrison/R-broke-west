#20170214: conceptual figure for Lisa's krill hurdle paper
library(grImport)
library(png)
library(ggplot2)
library(gridGraphics)
require(gridExtra)

wd='C:/Users/martin_cox/Documents/projects/Lisa/'

img1 <- readPNG(paste(wd,"krill.png",sep=''))
g1 <- rasterGrob(img1, interpolate=FALSE)

#Example environment
waterT=seq(1,0.9,length.out=30) ###################################
waterT=append(waterT,seq(waterT[length(waterT)],0.6,length.out=20)) ###################################
waterT=append(waterT,seq(waterT[length(waterT)],0.5,length.out=190)) ###################################
dat=data.frame(obs=1:length(waterT),waterT=waterT)
rm(waterT)
p1=qplot(obs,waterT,data=dat,geom='line')+
  labs(x='Observation number',y=expression(paste("Water temperature ,",degree,"C")),
       title='Example environment - water temperature')

#environmental pref
dat$drift=1/nrow(dat)
#assume depth pref is Gaussian 
dpMean=0.75;dpSD=0.05 ###################################
dat$pref=dnorm(dat$waterT,0.75,0.05)/sum(dnorm(dat$waterT,dpMean,dpSD))
maxy=max(c(dat$drift,dat$pref)) #for plotting
p2=qplot(obs,drift,data=dat,geom='line')+
  labs(x='Observation number',y='Environmental preference - drifting',
       title='Example preference - drifting')+
  scale_y_continuous(limits=c(0, maxy))

p3=qplot(obs,pref,data=dat,geom='line')+
  labs(x='Observation number',y='Environmental preference - swimming',
       title='Example preference - swimming')+
  scale_y_continuous(limits=c(0, maxy))


#observed distribution
#uniform (drifting):
nswarms=10   ###################################
xxu=runif(nswarms,1,nrow(dat))
yyu=dat$waterT[xxu]

#add little krill to plot (drifting):
xoffset=14 ###################################
yoffset=0.2 ###################################

p4=qplot(obs,waterT,data=dat,geom='line')+
  labs(x='Observation number',y=expression(paste("Water temperature ,",degree,"C")),
       title='As observed distribution - drifting')

for(i in 1:nswarms)
  p4=p4+annotation_custom(g1, xmin=xxu[i]-xoffset, xmax=xxu[i]+xoffset, ymin=yyu[i]-yoffset, ymax=yyu[i]+yoffset)


#Gaussian (swimming):
xxg=sample(x=1:nrow(dat),size=nswarms,replace=TRUE,prob=dat$pref)
yyg=dat$waterT[xxg]

#add little krill to plot (swimming):
p5=qplot(obs,waterT,data=dat,geom='line')+
  labs(x='Observation number',y=expression(paste("Water temperature ,",degree,"C")),
       title='As observed distribution - swimming')

for(i in 1:nswarms)
  p5=p5+annotation_custom(g1, xmin=xxg[i]-xoffset, xmax=xxg[i]+xoffset, ymin=yyg[i]-yoffset, ymax=yyg[i]+yoffset)

#plot
windows(width=8,height=12)
grid.arrange(p1,p1,p2,p3,p4,p5,nrow=3,ncol=2)
      


