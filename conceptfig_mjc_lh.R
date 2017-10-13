library(grImport)
library(png)
library(ggplot2)
library(gridGraphics)
require(gridExtra)

wd <- 'C:/Users/Lisa/Dropbox/uni/hurdle paper/figures/'


img1 <- readPNG(paste0(wd, "krill_new.png"))
transparent <- img1[,,4] == 0
img1 <- as.raster(img1[,,1:3])
img1[transparent] <- NA
img1[img1 == "#000606"] <- NA

g1   <- rasterGrob(img1, interpolate=FALSE)




#Example environment
waterT <- seq(0, 1, length.out=30) 
pdf <- dunif(waterT, min = 0, max = 1)
dat <- data.frame(pdf=pdf, waterT=waterT)

p1 <- qplot(waterT, pdf, data=dat, geom='line') +
  labs(y='PDF',x=expression(paste("Water temperature (",degree,"C)")), title="") +
  geom_text(size = 5, x = 0.9, y = 1.4, label = "(a)") +
  theme(legend.title=element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        panel.background = element_blank())


#environmental pref
dat$drift <- 1/nrow(dat)
#assume depth pref is Gaussian 
dpMean <- 0.5
dpSD <- 0.05 
dat$pref <- dnorm(dat$waterT, 0.5, 0.05)/sum(dnorm(dat$waterT, dpMean, dpSD))
maxy <- max(c(dat$drift,dat$pref))

p2 <- qplot(waterT, drift, data=dat, geom='line')+
  labs(x = expression(paste("Water temperature (",degree,"C)")), y = 'Probability of presence',
       title='Drifting - No preference') +
  scale_y_continuous(limits=c(0, maxy)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(size = 5, x = 0.9, y = 0.2, label = "(b)") +
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
       title=expression(paste("Swimming - Preference of ~0.5",degree,"C"))) +
  geom_text(size = 5, x = 0.9, y = 0.2, label = "(c)") +
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
temp_unif <- runif(nswarms, 0, 1)

#add little krill to plot (drifting):
p4 <- qplot(waterT, pdf, data=dat, geom='line')+
  labs(y='PDF', x=expression(paste("Water temperature (",degree,"C)")),
       title='Drifting') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(size = 5, x = 0.9, y = 1.4, label = "(d)") +
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
  p4 <- p4 + annotation_custom(g1, xmin=temp_unif[i]-0.2, xmax=temp_unif[i], ymin=0.9, ymax=1.1)
}


#Gaussian (swimming):
x_swim <- rnorm(nswarms, 0.5, 0.05)

#add little krill to plot (swimming):
p5 <- qplot(waterT, pdf, data=dat, geom='line')+
  labs(y="", x=expression(paste("Water temperature (",degree,"C)")),
       title='Swimming') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(size = 5, x = 0.9, y = 1.4, label = "(e)") +
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
  p5 <- p5 + annotation_custom(g1, xmin=x_swim[i]-0.1, xmax=x_swim[i]+0.1, ymin=0.9, ymax=1.1)
}

png("C:/Users/Lisa/Dropbox/uni/paperwork/examination/plots/concept_fig.png", width = 9, height = 12, units = "in", res = 100)

#plot
plus <- textGrob("+", gp=gpar(fontface="bold", fontsize = 40))
equals <- textGrob("=", gp=gpar(fontface="bold", fontsize = 40))
blank <- grid.rect(gp=gpar(col="white"))

grid.arrange(arrangeGrob(blank, p1, ncol = 3, widths=c(1, 2, 1), top = textGrob("Example Environment - Water Temperature", gp=gpar(fontface="bold"))), 
             plus, arrangeGrob(p2, p3, ncol = 2, top = textGrob("Example Preference", gp=gpar(fontface="bold"))), equals,
             arrangeGrob(p4, p5, ncol = 2, top = textGrob("Observed Distribution", gp=gpar(fontface="bold"))), nrow=5, ncol=1, heights=c(3, 1, 3, 1, 3))

dev.off()
