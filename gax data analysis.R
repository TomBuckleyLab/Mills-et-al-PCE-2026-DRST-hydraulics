library(ggplot2)
library(dplyr)
library(viridis)
library(tidyr)

options("scipen"=100, "digits"=8)

## load in data
x=read.delim("gas exchange data.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
colnames(x)[1] <- "rep" ## something weird about how it loaded in this column name; fix


## correct D to deltaw (kPa to mmol/mol)
x$d <- 10*x$d


## get median D and T by experiment
x <- x 	%>% group_by(spp, rep) %>% 
  mutate(	medianD = median(d),
          medianT = median(t))


## add columns to say whether a given row was the high or low D or T for that expt
x$lohiD <- ifelse(x$d>x$medianD, "D", "d") ## D=high deltaw; d=low
x$lohiT <- ifelse(x$t>x$medianT, "T", "t") ## T=high T; t=low
x$trt <- paste0(x$lohiD, x$lohiT)


## create solumn with unique identifier for each plant/expt
x$spp_rep <- paste0(x$spp, "_", x$rep)

## discard extraneous columns
x <- x[,c("g", "d", "t", "trt", "spp", "spp_rep")]


## convert to wide format so that the four experimental values
##  for each day (dt, dT, etc) are in separate columns
y <- pivot_wider(data=x,
                 names_from="trt",
                 values_from=c("g", "d", "t"))


## calculations 
y$dT.d <- y$t_dT - y$t_dt ## T change at low D
y$dT.D <- y$t_DT - y$t_Dt ## T change at high D
y$dT <- 0.5*(y$dT.d + y$dT.D) ## avg T change

y$d <- 0.5*(y$d_dt + y$d_dT) ## average of low D between low and high T
y$D <- 0.5*(y$d_Dt + y$d_DT) ## average of high D between low and high T

y$g.d <- 0.5*(y$g_dt + y$g_dT) ## avg g at low D
y$g.D <- 0.5*(y$g_Dt + y$g_DT) ## avg g at high D

y$dlng.d <- (y$g_dT - y$g_dt)/y$g.d ## relative change in g at low D
y$dlng.D <- (y$g_DT - y$g_Dt)/y$g.D ## rel change in g at high D

y$dlngdt.d <- y$dlng.d/y$dT.d ## dlng/dT at low D
y$dlngdt.D <- y$dlng.D/y$dT.D ## dlng/dT at high D


## difference in dlng/dT between high and low D. this is the main output
##  multiply by 100 to express it as a percent DRST
##  multiply by 10 (as a standard 10C shift in T)
y$ddg <- 100*(y$dlngdt.D - y$dlngdt.d)*10  


### now do theoretical predictions based on K
q=read.delim("q10_distribution.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
m=read.delim("oren_distribution.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)


### make df with randomly resampled values of Oren slope and resulting K/b (stom sensitivity parameter)
nq <- 10000
qm <- data.frame(matrix(nrow=nq, ncol=0))
qm$m <- sample(m$m, nrow(qm), replace=T)
## K/b = 10*(1/m - 1) where m is the unitless version of the Oren slope
qm$kb <- 10*(1/qm$m - 1)


## get mean low and hi delta-w and lo and hi T across expts for each spp
sppavg <- y %>% group_by(spp) %>% summarize(dhi = mean(D), dlo = mean(d),
                                            tlo = mean(mean(t_dt),mean(t_Dt)),
                                            thi = mean(mean(t_dT),mean(t_DT)))

sppavg$dT <- sppavg$thi - sppavg$tlo ## mean T shift for each spp

us <- unique(x$spp) ## list of species codes

## stack qm on itself 6 times (one for each spp); we will fill w/ randomly sampled Q10s
qm2 <- rbind(qm,qm)
qm4 <- rbind(qm2,qm2)
qm <- rbind(qm4,qm2)

## create column w/ spp names
qm$spp <- c(rep(us[1],nq),rep(us[2],nq),rep(us[3],nq),rep(us[4],nq),rep(us[5],nq),rep(us[6],nq))

## pull in the avg deltaw and T values from earlier
qm <- left_join(qm, sppavg, by="spp")

## make three copies, to hold values based on root, leaf and stem Q10 values
qmr <- qm
qml <- qm
qms <- qm

## isolate literature Q10s for roots and stems
q.r <- subset(q, type=="r")
q.l <- subset(q, type=="l")

## randomly sample Q10s
qmr$q10 <- sample(q.r$Q10, nq, replace=T); qmr$organ <- "root"
qml$q10 <- sample(q.l$Q10, nq, replace=T); qml$organ <- "leaf"
qms$q10 <- 1.24; qms$organ <- "stem" ## 1.24 is assumed stem value (for viscosity)

## calculate predicted ddg
##   (difference in dlng/dT between high and low D) * standard 10C T shift
qmr$ddg.pred <- 100*((qmr$dhi/(qmr$kb + qmr$dhi) - qmr$dlo/(qmr$kb + qmr$dlo))*((1/(qmr$q10^(-qmr$dT/10)) - 1)/qmr$dT))*10
qml$ddg.pred <- 100*((qml$dhi/(qml$kb + qml$dhi) - qml$dlo/(qml$kb + qml$dlo))*((1/(qml$q10^(-qml$dT/10)) - 1)/qml$dT))*10
qms$ddg.pred <- 100*((qms$dhi/(qms$kb + qms$dhi) - qms$dlo/(qms$kb + qms$dlo))*((1/(qms$q10^(-qms$dT/10)) - 1)/qms$dT))*10

## pull results back together
qm <- rbind(qmr, qml)
qm <- rbind(qm, qms)


### calculate means and medians by spp and Q10 type
qmspp <- qm %>% group_by(spp, organ) %>% summarize(ddg.mean = mean(ddg.pred),
                                          ddg.median = median(ddg.pred))

## order spp for presentation
qmspp$spp <- factor(qmspp$spp, levels=c("vr", "vvpn", "ep", "hv", "zm", "pv"))
qm$spp <- factor(qm$spp, levels=c("vr", "vvpn", "ep", "hv", "zm", "pv"))


## quantiles
qmq <- qm %>% group_by(spp, organ) %>% summarize(q90 = quantile(ddg.pred, 0.9),
                                                 q95 = quantile(ddg.pred, 0.95),
                                                 q99 = quantile(ddg.pred, 0.99),
                                                 q995 = quantile(ddg.pred, 0.995))

## overall avg by organ
qmsppavg <- qm %>% group_by(organ) %>% summarize(md = mean(ddg.pred))



#############
############
#### figure with generic continuous predictions vs Dw

### need dlng/dT per unit dlnK/dT as function of Dw
### this is Dw/(K/b + Dw)
### repeat for three values of K/b (3.3, 6.7 and 12.6 mmol/mol)
### and three values of dlnK/dT (Q10 = 1.24, 1.525, 1.955; medians for stem, leaf, root respectively)

z <- data.frame(matrix(nrow=100, ncol=0))
z$deltaw <- seq(from=5, to=30, by=25/99)
zo <- z
z <- rbind(z,z); z <- rbind(z,z); z <- rbind(z,z); z <- rbind(z,zo);
z$kb <- c(rep(3.3,300), rep(6.7,300), rep(12.6,300))

#Q10s: means in each dist
q10r <- round(mean(q.r$Q10),2); q10l <- round(mean(q.l$Q10),2); q10s = 1.24
q10seq <- c(rep(q10s,100), rep(q10l,100), rep(q10r,100))
z$Q10 <- c(q10seq, q10seq, q10seq)
z$dlnK.dT <- (z$Q10 - 1)/10 ## K @ 30C/K @ 20C = Q10; dlnK = dK/K @ 20C = Q10 - 1; dlnK/dT = this/10
z$dlng.dT <- z$dlnK.dT*z$deltaw/(z$deltaw + z$kb)
z$dlng10 <- 100*10*z$dlng.dT ## 10C standardized T shift
z$sensitivity <- ifelse(z$kb==3.3, "high", ifelse(z$kb==6.7, "medium", "low")) ## text tag for legend
z$Q10 <- as.factor(z$Q10)
z$Q10 <- factor(z$Q10, levels=c(q10r, q10l, q10s))
z$sensitivity <- factor(z$sensitivity, levels=c("high", "medium", "low"))


ggplot() +
   geom_line(data=z, aes(x=deltaw, y=dlng10, linetype=sensitivity, color=Q10), size=1) +
   theme_bw() + 
   scale_color_manual(values=c("#2a788e", "#7ad151", "#440154")) +
   xlab(expression(paste(Delta,"w (mmol ",mol^-1,")"))) +
   ylab(expression(paste("predicted % change in ",italic("g")["s"]," with ",10^"o","C warming"))) +
   geom_vline(xintercept=13, alpha=0.15, color="black", size=2)+
   geom_vline(xintercept=22, alpha=0.15, color="black", size=2) +
   geom_point(aes(x=13, y=36.2), color="red", size=3, alpha=0.3)+
   geom_point(aes(x=22, y=42.15), color="red", size=3, alpha=0.3)+
   guides(linetype=guide_legend(title="stomatal\nsensitivity"))


ggsave("fig 1 (theoretical predictions vs dw).png", device="png", dpi=600)



###############
###############
####  main results figure

## SE function, needed for summaries
SE <- function(x, na.rm=FALSE) {
   if (na.rm) x <- na.omit(x)
   sqrt(var(x)/length(x))
}


## generate df with means and SEs of diff (between high and low Dw) in % change in gs w/ 10C shift ('ddg') by species
yy <- y %>% group_by(spp) %>%
   summarize(ddg.mean = mean(ddg),
             ddg.se = SE(ddg),
             ddg.median = median(ddg))

yy$spp <- factor(yy$spp, levels=c("vr", "vvpn", "ep", "hv", "zm", "pv"))


## replace spp codes with full names
(oldnames = unique(yy$spp))  ## enclosing code in parens forces output to display

newnames = c(expression(paste(italic("V. riparia"))),
             expression(paste(italic("V. vinifera"))),
             expression(paste(italic("E. polyanthemos"))), ## \n is line break
             expression(paste(italic("H. vulgare"))),
             expression(paste(italic("Z. mays"))),
             expression(paste(italic("P. vulgaris"))))


ggplot() +
   scale_color_viridis(begin=0.8, end=0, discrete=T) +
   scale_fill_viridis(begin=0.8, end=0, discrete=T) +
   # violin plots of theoretical predictions, by Q10 source (organ)
   geom_violin(data=qm, aes(x=ddg.pred, fill=organ, y=spp), 
               position=position_dodge(0.8), color="transparent",
               alpha=0.5, linewidth=0, trim=T, scale="width") +
   # use error bars w/ negligible SEs to show means and medians separately for each spp
   geom_errorbar(data=qmspp, aes(xmin=ddg.mean-0.01, xmax=ddg.mean+0.01, y=spp, color=organ), 
                 position=position_dodge(0.8), width=0.5)+
   geom_errorbar(data=qmspp, aes(xmin=ddg.median-0.01, xmax=ddg.median+0.01, y=spp, color=organ), 
                 position=position_dodge(0.8), width=0.5, linetype="dashed")+
   geom_errorbar(data=yy, aes(y=spp, xmin=ddg.mean-0.01, xmax=ddg.mean+0.01), 
                 stat="identity", position=position_dodge(), 
                 color="black", alpha=0.8, linewidth=1, width=0.5) +
   # use actual error bars to show +- SE of observations for each spp
   geom_errorbar(data=yy, aes(y=spp, xmin=ddg.mean-ddg.se, xmax=ddg.mean+ddg.se), 
                 width=0, alpha=0.3, linewidth=15, color="black")+
   # reference line for zero
   geom_vline(xintercept=0, color="black", linewidth=0.5)+
   xlab(expression(atop(paste("difference (between high and low ",Delta,italic("w"),")"),
                       paste("in % increase in ",italic("g")["s"]," following 10"^"o","C warming"))
                  )) +
   scale_y_discrete(labels=newnames) + 
   theme(axis.text.y=element_text(face="italic", hjust=0))+
   guides(fill=guide_legend(title="predictions\nby organ"),
          color=guide_legend(title="predictions\nby organ"))+
   ylab("species") +
   xlim(-11,45)+
   theme_bw()


ggsave("fig 4 (main results).png", device="png", dpi=600)



#############
##### figure showing distributions of leaf and root Q10s, and K/b

q$type <- ifelse(q$type=="l", "leaf", "root")
qs <- data.frame(matrix(nrow=2, ncol=0))
qs$Q10 <- 1.24
qs$type <- "stem"
q <- rbind(q, qs)
q <- subset(q, type != "stem")
pq <- ggplot() +
   geom_density(data=q, aes(x=Q10, fill=type), color="transparent", alpha=0.5) +
   scale_fill_viridis(begin=0.8, end=0.4, discrete=T) +
   guides(fill=guide_legend(title="organ"))+
   theme_bw() + theme(legend.position = c(0.7, 0.7))+
   xlab(expression(paste(italic("Q")["10"]," (unitless)"))); pq

m$kb <- 10*(1/m$m - 1)
pkb <- ggplot() +
   geom_density(data=m, aes(x=kb), fill="black", alpha = 0.3, color="transparent") +
   theme_bw() + xlim(0, 25)+
   xlab(expression(paste(italic("K/b")))); pkb

library(ggpubr)
ggarrange(pkb, pq, nrow=1, ncol=2, labels = c("(a)", "(b)"),
          label.x = c(0.15, 0.14), label.y=c(0.97, 0.97))

ggsave("fig 2 (q10 and Kb distributions).png", device="png", dpi=600)



