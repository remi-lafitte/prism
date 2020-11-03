#McIntosh, Brown & Young (2018)--------------
#"Meta-analysis of the visuospatial aftereffects of prism adaptation, with two novel experiments"
#ANALYSIS CODE, META_ANALYSIS
#LOAD REQUIRED LIBRARIES--------------
library(ggplot2)
library(gridExtra)
library(metafor)
library(ggpmisc)
library(tidyverse)
library(magrittr)
#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv---------
#load the data
meta <- read.csv(here::here("McINTOSH_19_PA_META","PSA_META.csv"), sep = ";", dec = ",")

#see the data structure
str(meta)

#calculate standard error per study
meta$SE <- 1/sqrt(meta$n)

#CORRELATIONS to check possible predictors--------
cor(meta[c("duration_c","prism", "number","deg", "d")], method="p", use="complete.obs")
# Moderate corrrelations between d and prism/number of movement.

#formula for scattergrams
my.formula <- y ~ x

#LANDMARK: SCATTERGRAM RAW_STANDARD ES-----------
ggplot(meta, aes(x=deg, y=d)) +
  geom_point(size = 4, alpha=.6) +
  ggtitle("(a)") +
  geom_smooth(method="lm", colour="black", se = FALSE, fullrange = TRUE) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +     
  scale_x_continuous(limits = c(0, 7), breaks = seq(from = 0, to=7, by=1)) +
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(from = 0, to=2.5, by=.5)) +
  labs(x = 'Raw effect size (cm)', y= 'Standardised effect size (d)') +
  geom_hline(yintercept=0, linetype="dotted") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major =element_blank(), panel.grid.minor =element_blank(),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.23, .9),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank()) -> PA_META_Fig1
PA_META_Fig1

#META_ANALYSIS--------------

#simple random effects models
ma <- rma.uni(yi=d, sei=SE, data=meta, method= "REML")

#UNMODERATED FUNNEL PLOT, LANDMARK--------------

#GET META_ESTIMATE
estimate<- as.numeric(summary(ma)[[1]])
estimate
#get Tau_squared
T <- as.numeric(summary(ma)[9])
T <- paste("tau^2 == ", round(T, 2))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(meta[!is.na(meta$d), "SE"]), 0.001)

#95% CI vectors
ll95 <- estimate - (1.96*se.seq)
ul95 <- estimate + (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#LANDMARK FUNNEL PLOT----------
ggplot(meta[!is.na(meta$d), ], aes(x = d, y = SE)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=estimate) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=number, shape = as.factor(prism))) +
  # scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  # scale_colour_gradient(low="grey80", high="black", name="Prism") +
  labs(x= 'Standardised effect size (d)', y= 'Standard error') +
  scale_y_reverse() +
  scale_x_continuous(limits = c(0,1.6), breaks=seq(0,1.6,0.5)) +
  ggtitle("(b)") +
  annotate("text", x = -.5, y = 0, label = "Total heterogeneity:", hjust = "left") +
  annotate("text", x = -.5, y = 0.025, label = T, parse = TRUE, hjust = "left") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank(),
        panel.grid.minor.x =element_blank(), axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.88, .78),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig2

PA_META_Fig2


#random effects models, with moderators-------
ma <- rma.uni(yi=d, sei=SE, data=meta, mods=~prism+number, method= "REML"); summary(ma)

#MODERATED FUNNEL PLOT LANDMARK------------
#meta-analytic estimate and CIs
fun <- funnel(ma)
fun <- cbind(fun, meta[!is.na(meta$d), c("prism", "number")])

#get Tau_squared
T <- as.numeric(summary(ma)[9])
T <- paste("tau^2 == ", round(T, 2))

#get asymmetry
EGG <- paste0("z = ", round(regtest(ma)$zval, 2), ", p = ", round(regtest(ma)$pval, 2))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(fun$y), 0.001)

#95% CI vectors
ll95 <- - (1.96*se.seq)
ul95 <- (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#LANDMARK FUNNEL PLOT
ggplot(fun, aes(x = x, y = y)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=0) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=number, shape=as.factor(prism))) +
  # scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  # scale_colour_gradient(low="grey80", high="black", name="Prism") +
  labs(x= 'Residual effect size (d)', y= 'Standard error') +
  scale_y_reverse() +
  scale_x_continuous(limits = c(-1.25,1.25), breaks=seq(-1,1,0.5)) +
  ggtitle("(e)") +
  annotate("text", x = -1.25, y = 0, label = "Residual heterogeneity:", hjust = "left") +
  annotate("text", x = -1.25, y = 0.025, label = T, parse = TRUE, hjust = "left") +
  annotate("text", x = -1.25, y = 0.06, label = "Asymmetry:", hjust = "left") +
  annotate("text", x = -1.25, y = 0.085, label = EGG, hjust = "left") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank(),
        panel.grid.minor.x =element_blank(), axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.88, .78),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig3

PA_META_Fig3

#UNPACK MODEL PREDICTIONS-----------------
#make df for results
predict_MA<-as.data.frame(predict.rma(ma, level = 95, addx = TRUE))

#plot meta-analytic models
ggplot(predict_MA, aes(x=X.number, y= pred, col = as.factor(X.prism))) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_point(size = 5, alpha = 0.3)+
  scale_color_manual(name= "Prism (°)", values = c("8.5", "10", "12.4",
                                                   "15", "17")) +
  # geom_line(size=1)+
  # geom_ribbon(aes(ymin=ci.lb,ymax=ci.ub), fill="grey", alpha="0.5") +
  geom_smooth(colour="black", method="lm", fill = "grey") +
  geom_smooth(aes(y=cr.lb), linetype="dotted", colour="black", method="lm",
              fullrange = TRUE) +
  geom_smooth(aes(y=cr.ub), linetype="dotted", colour="black", method="lm",
              fullrange = TRUE) +
  scale_y_continuous(limits=c(0,2), breaks = seq(from=-0, to=2, by=.5)) +
  labs(x="Number of pointing", y = "Predicted average effect size (d)") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
        axis.title = element_text(size=18), axis.text = element_text(size=12),
        strip.text = element_text(size=16))-> PA_META_Fig4
PA_META_Fig4



#prediction for 10degrees PA and 250 movements-------
ma$beta[1,1] + (ma$beta[2,1]*10) + (ma$beta[3,1]*250)
# 0.74 
ma$beta[1,1] + (ma$beta[2,1]*10) + (ma$beta[3,1]*200)

ma$beta[1,1] + (ma$beta[2,1]*10) + (ma$beta[3,1]*200)

(predict_MA %>% arrange(X.number))
newmods <- cbind(number=seq(from = 200, to = 400, by=20), prism = 10)
p<-as.data.frame(predict.rma(ma, newmods=newmods, level = 95, addx = TRUE))
ggplot(data = p, aes(x = X.number, y = pred))+
  geom_line(col = "red", size=2)+
  geom_smooth(aes(y=cr.lb), linetype="dotted", colour="black", method="lm") +
  geom_smooth(aes(y=cr.ub), linetype="dotted", colour="black", method="lm")+
  geom_smooth(aes(y=ci.lb), linetype="longdash", colour="blue", method="lm") +
  geom_smooth(aes(y=ci.ub), linetype="longdash", colour="blue", method="lm")+
  scale_y_continuous(breaks=seq(0,2,0.1))+
  scale_x_continuous(breaks=seq(200,400,50))+
  labs(x = "Number of movements", y = "Cohen's d prediction")+
  # geom_vline(xintercept = c(200,250,300))+
  theme_bw(base_size = 20)-> PA_META_Fig5
p %>% filter(X.number %in% c(200,260,300))
  
#forest plot------------------
forest(ma)

#WRITE FILES-----------
png("PA_META_Fig.png", units="in", width=10.5, height=15.75, res=200)
grid.arrange(PA_META_Fig1, PA_META_Fig2, PA_META_Fig3, PA_META_Fig4,
             PA_META_Fig5, ncol=2)
dev.off()

#tidy environment------------------
rm(dfCI, fun, newmods,predict_MA, EGG, estimate, ll95, ma, my.formula, se.seq, T, t, ul95)
