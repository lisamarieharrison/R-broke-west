#Conceptual figure for krill hurdle paper
#author: Martin Cox, with modifications by Lisa-Marie Harrison
#date: 14/02/2017


library(grImport)
library(png)
library(ggplot2)
library(gridGraphics)
require(gridExtra)

wd <- 'C:/Users/43439535/Dropbox/uni/hurdle paper/figures/'

img1 <- readPNG(paste0(wd, "krill_new.png"))
transparent <- img1[,,4] == 0
img1 <- as.raster(img1[,,1:3])
img1[transparent] <- NA
img1[img1 == "#000606"] <- NA

g1   <- rasterGrob(img1, interpolate=FALSE)




#Example environment
waterT <- seq(1, 0.9, length.out=30) 
waterT <- append(waterT, seq(waterT[length(waterT)], 0.6, length.out=20)) 
waterT <- append(waterT, seq(waterT[length(waterT)], 0.5, length.out=190)) 
dat <- data.frame(obs=1:length(waterT), waterT=waterT)

p1 <- qplot(obs, waterT, data=dat, geom='line') +
  labs(x='Observation number',y=expression(paste("Water temperature (",degree,"C)"))) +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())

p1b <- p1 +  labs(x='Observation number', y="")

#environmental pref
dat$drift <- 1/nrow(dat)
#assume depth pref is Gaussian 
dpMean <- 0.75
dpSD <- 0.05 
dat$pref <- dnorm(dat$waterT, 0.75, 0.05)/sum(dnorm(dat$waterT, dpMean, dpSD))
maxy <- max(c(dat$drift,dat$pref))

p2 <- qplot(waterT, drift, data=dat, geom='line')+
  labs(x = expression(paste("Water temperature (",degree,"C)")), y = 'Probability of presence',
       title='Drifting - No preference') +
  scale_y_continuous(limits=c(0, maxy)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())

p3 <- qplot(waterT, pref, data=dat, geom='line') +
  labs(x = expression(paste("Water temperature (",degree,"C)")), y = '',
       title=expression(paste("Swimming - Preference of ~0.7",degree,"C"))) +
  scale_y_continuous(limits=c(0, maxy)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())


#observed distribution
#uniform (drifting):
nswarms <- 10  
xxu <- runif(nswarms, 1, nrow(dat))
yyu <- dat$waterT[xxu]

#add little krill to plot (drifting):
xoffset <- 14 
yoffset <- 0.2

p4 <- qplot(obs, waterT, data=dat, geom='line')+
  labs(x='Observation number', y=expression(paste("Water temperature (",degree,"C)")),
       title='Drifting') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())

for (i in 1:nswarms) {
  p4 <- p4 + annotation_custom(g1, xmin=xxu[i]-xoffset, xmax=xxu[i]+xoffset, ymin=yyu[i]-yoffset, ymax=yyu[i]+yoffset)
}


#Gaussian (swimming):
xxg <- sample(x=1:nrow(dat), size=nswarms, replace=TRUE, prob=dat$pref)
yyg <- dat$waterT[xxg]

#add little krill to plot (swimming):
p5 <- qplot(obs, waterT, data=dat, geom='line')+
  labs(x='Observation number', y="",
       title='Swimming') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())

for (i in 1:nswarms) {
  p5 <- p5 + annotation_custom(g1, xmin=xxg[i]-xoffset, xmax=xxg[i]+xoffset, ymin=yyg[i]-yoffset, ymax=yyg[i]+yoffset)
}

#plot
plus <- textGrob("+", gp=gpar(fontface="bold", fontsize = 40))
equals <- textGrob("=", gp=gpar(fontface="bold", fontsize = 40))
windows(width=9, height=12)
grid.arrange(arrangeGrob(p1, p1b, ncol = 2, top = textGrob("Example Environment - Water Temperature", gp=gpar(fontface="bold"))), 
             plus, arrangeGrob(p2, p3, ncol = 2, top = textGrob("Example Preference", gp=gpar(fontface="bold"))), equals,
             arrangeGrob(p4, p5, ncol = 2, top = textGrob("Observed Distribution", gp=gpar(fontface="bold"))), nrow=5, ncol=1, heights=c(3, 1, 3, 1, 3))

