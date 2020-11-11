#McIntosh, Brown & Young (2018)
#"Meta-analysis of the visuospatial aftereffects of prism adaptation, with two novel experiments"
#ANALYSIS CODE, META_ANALYSIS
#LOAD REQUIRED LIBRARIES
library(ggplot2)
library(gridExtra)
library(metafor)
library(ggpmisc)

#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv
#load the data
meta <- read.csv("PA_META.csv")
#restrict to leftward prisms
meta <- meta[meta$direction == "L", ]

#dichotomise exposure duration
meta$expose <- meta$duration >= 10

#impute tabletop distance
meta[is.na(meta$distance), "distance"] <- 375

#encode line length as degrees of visual angle (dva)
meta$dva_LM <- atan(meta$length_LM/meta$distance)*(180/pi)
meta$dva_LB <- atan(meta$length_LB/meta$distance)*(180/pi)

#calculate standard error per study
meta$SE <- 1/sqrt(meta$n)

#CORRELATIONS to check possible predictors
cor(meta[c(5,6,16,17,11)], method="p", use="complete.obs")
cor(meta[c(5,6,16,18,14)], method="p", use="complete.obs")

#####################
#formula for scattergrams
my.formula <- y ~ x

#LANDMARK: SCATTERGRAM RAW_STANDARD ES
ggplot(meta, aes(x=raw_LM, y=d_LM)) +
  geom_point(size = 4, alpha=.6) +
  ggtitle("(a)") +
  geom_smooth(method="lm", colour="black", se = FALSE, fullrange = TRUE) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +     
  scale_x_continuous(limits = c(-.9, 4.75), breaks = seq(from = 0, to=4, by=1)) +
  scale_y_continuous(limits = c(-.5, 2.4), breaks = seq(from = -.5, to=2.5, by=.5)) +
  labs(x = 'Raw effect size (% of half-length)', y= 'Standardised effect size (d)') +
  geom_hline(yintercept=0, linetype="dotted") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major =element_blank(), panel.grid.minor =element_blank(),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.23, .9),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank()) -> PA_META_Fig4a

#BISECTION: SCATTERGRAM RAW_STANDARD ES
ggplot(meta, aes(x=raw_LB, y=d_LB)) +
  geom_point(size = 4, alpha=.6) +
  ggtitle("(b)") +
  geom_smooth(method="lm", colour="black", se = FALSE, fullrange = TRUE) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +     
  scale_x_continuous(limits = c(-.9, 4.75), breaks = seq(from = 0, to=4, by=1)) +
  scale_y_continuous(limits = c(-.5, 2.4), breaks = seq(from = -.5, to=2.5, by=.5)) +
  labs(x = 'Raw effect size (% of half-length)', y= 'Standardised effect size (d)') +
  geom_hline(yintercept=0, linetype="dotted") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major =element_blank(), panel.grid.minor =element_blank(),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.23, .9),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank()) -> PA_META_Fig4b

##############
#META_ANALYSIS

#simple random effects models
ma_LM <- rma.uni(yi=d_LM, sei=SE, data=meta, method= "REML")
ma_LB <- rma.uni(yi=d_LB, sei=SE, data=meta, method= "REML")

###############################
#UNMODERATED FUNNEL PLOT, LANDMARK

#GET META_ESTIMATE
estimate<- as.numeric(summary(ma_LM)[[1]])

#get Tau_squared
T <- as.numeric(summary(ma_LM)[9])
T <- paste("tau^2 == ", round(T, 2))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(meta[!is.na(meta$d_LM), "SE"]), 0.001)

#95% CI vectors
ll95 <- estimate - (1.96*se.seq)
ul95 <- estimate + (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#LANDMARK FUNNEL PLOT
ggplot(meta[!is.na(meta$d_LM), ], aes(x = d_LM, y = SE)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=estimate) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=prism, shape=expose)) +
  scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  scale_colour_gradient(low="grey80", high="black", name="Prism") +
  labs(x= 'Standardised effect size (d)', y= 'Standard error') +
  scale_y_reverse() +
  scale_x_continuous(limits = c(-.5,1.5), breaks=seq(-.5,1.5,0.5)) +
  ggtitle("(c)") +
  annotate("text", x = -.5, y = 0, label = "Total heterogeneity:", hjust = "left") +
  annotate("text", x = -.5, y = 0.025, label = T, parse = TRUE, hjust = "left") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank(),
        panel.grid.minor.x =element_blank(), axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.88, .78),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig4c

#UNMODERATED FUNNEL PLOT, BISECTION

#GET META_ESTIMATE
estimate<- as.numeric(summary(ma_LB)[[1]])

#get Tau_squared
T <- as.numeric(summary(ma_LB)[9])
T <- paste("tau^2 == ", round(T, 2))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(meta[!is.na(meta$d_LB), "SE"]), 0.001)

#95% CI vectors
ll95 <- estimate - (1.96*se.seq)
ul95 <- estimate + (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#BISECT FUNNEL PLOT
ggplot(meta[!is.na(meta$d_LB), ], aes(x = d_LB, y = SE)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=estimate) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=prism, shape=expose)) +
  scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  scale_colour_gradient(low="grey80", high="black", name="Prism") +
  labs(x= 'Standardised effect size (d)', y= 'Standard error') +
  scale_y_reverse() +
  scale_x_continuous(limits = c(-0.5,2.5), breaks=seq(-0.5,2.5,0.5)) +
  ggtitle("(d)") +
  annotate("text", x = -.5, y = 0, label = "Total heterogeneity:", hjust = "left") +
  annotate("text", x = -.5, y = 0.025, label = T, parse = TRUE, hjust = "left") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank(),
        panel.grid.minor.x =element_blank(), axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.88, .78),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig4d

##########################
#random effects models, with moderators
ma_LM <- rma.uni(yi=d_LM, sei=SE, data=meta, mods=~prism+expose, method= "REML"); summary(ma_LM)
ma_LB <- rma.uni(yi=d_LB, sei=SE, data=meta, mods=~prism+expose, method= "REML"); summary(ma_LB)
##########################
#for exploratory comparison, removing Experiment 2
summary(rma.uni(yi=d_LB, sei=SE, data=meta[meta$STUDY != "McIntosh2018", ], mods=~prism+expose, method= "REML"))

################
#MODERATED FUNNEL PLOT LANDMARK
#meta-analytic estimate and CIs
fun_LM <- funnel(ma_LM)
fun_LM <- cbind(fun_LM, meta[!is.na(meta$d_LM), c("prism", "expose")])

#get Tau_squared
T <- as.numeric(summary(ma_LM)[9])
T <- paste("tau^2 == ", round(T, 2))

#get asymmetry
EGG <- paste0("z = ", round(regtest(ma_LM)$zval, 2), ", p = ", round(regtest(ma_LM)$pval, 2))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(fun_LM$y), 0.001)

#95% CI vectors
ll95 <- - (1.96*se.seq)
ul95 <- (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#LANDMARK FUNNEL PLOT
ggplot(fun_LM, aes(x = x, y = y)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=0) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=prism, shape=expose)) +
  scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  scale_colour_gradient(low="grey80", high="black", name="Prism") +
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
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig4e

##############
#MODERATED FUNNEL PLOT BISECT
#meta-analytic estimate and CIs
fun_LB <- funnel(ma_LB)
fun_LB <- cbind(fun_LB, meta[!is.na(meta$d_LB), c("prism", "expose")])

#get Tau_squared
T <- as.numeric(summary(ma_LB)[9])
T <- paste("tau^2 == ", round(T, 2))

#get asymmetry
EGG <- paste0("z = ", round(regtest(ma_LB)$zval, 2), ", p = ", round(regtest(ma_LB)$pval, 3))

#SE vector from zero to max SE (min n), with small increment
se.seq <- seq(0, max(fun_LB$y), 0.001)

#95% CI vectors
ll95 <- - (1.96*se.seq)
ul95 <- (1.96*se.seq)

#data frame of CIs
dfCI = data.frame(ll95, ul95, se.seq)

#BISECT FUNNEL PLOT
ggplot(fun_LB, aes(x = x, y = y)) +
  geom_line(data = dfCI, aes(x = ll95, y = se.seq), colour = 'grey') +
  geom_line(data = dfCI, aes(x = ul95, y = se.seq), colour = 'grey') +
  geom_vline(xintercept=0) +
  #plot the data
  geom_point(size = 4, alpha=.9, aes(colour=prism, shape=expose)) +
  scale_shape_manual(values=c(19,17), name="Exposure", labels=c("short", "long")) +
  scale_colour_gradient(low="grey80", high="black", name="Prism") +
  labs(x= 'Residual effect size (d)', y= 'Standard error') +
  scale_y_reverse(limits = c(.55, 0), breaks=seq(0,5,.1)) +
  scale_x_continuous(limits = c(-1.25,1.25), breaks=seq(-1,1,0.5)) +
  ggtitle("(f)") +
  annotate("text", x = -1.25, y = 0, label = "Residual heterogeneity:", hjust = "left") +
  annotate("text", x = -1.25, y = 0.028, label = T, parse = TRUE, hjust = "left") +
  annotate("text", x = -1.25, y = 0.073, label = "Asymmetry:", hjust = "left") +
  annotate("text", x = -1.25, y = 0.105, label = EGG, hjust = "left") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major.y =element_blank(), panel.grid.minor.y =element_blank(),
        panel.grid.minor.x =element_blank(), axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.88, .78),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) -> PA_META_Fig4f

###############
#UNPACK MODEL PREDICTIONS
#make df for results
predict_MA <- read.csv(text="pred,se,ci.lb,ci.ub,cr.lb,cr.ub,X.intrcpt,X.prism,X.time,TASK,EXPOSE")

for(t in c(FALSE, TRUE)){
  newmods <- cbind(prism=seq(from = 8, to = 17, by=.5), expose=t)
  pLM <- as.data.frame(predict.rma(ma_LM, newmods=newmods, level = 95, addx = TRUE))
  pLB <- as.data.frame(predict.rma(ma_LB, newmods=newmods, level = 95, addx = TRUE))
  pLM$TASK <- "Landmark task"
  pLB$TASK <- "Line bisection"
  pLM$EXPOSE <- t
  pLB$EXPOSE <- t
  predict_MA <- rbind(predict_MA, (rbind(pLM, pLB)))
}

#relevel & rename exposure factor
predict_MA$EXPOSE <- factor(predict_MA$EXPOSE, levels=c(FALSE, TRUE),
                            labels = c("Short exposure (< 10 mins)", "Long exposure (10 mins +)"))

#plot meta-analytic models
ggplot(predict_MA, aes(x=X.prism, y= pred)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_line(size=1)+
  geom_ribbon(aes(ymin=ci.lb,ymax=ci.ub), fill="grey", alpha="0.5") +
  geom_smooth(aes(y=cr.lb), linetype="dotted", colour="black") +
  geom_smooth(aes(y=cr.ub), linetype="dotted", colour="black") +
  scale_y_continuous(limits=c(-1.5,2.5), breaks = seq(from=-1.5, to=2.5, by=.5)) +
  labs(x="Prism strength (degrees)", y = "Predicted average effect size (d)") +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(),
        axis.title = element_text(size=18), axis.text = element_text(size=12),
        strip.text = element_text(size=16)) +
  facet_grid(EXPOSE~TASK) -> PA_META_Fig5

############
#WRITE FILES
png("PA_META_Fig4.png", units="in", width=10.5, height=15.75, res=200)
grid.arrange(PA_META_Fig4a, PA_META_Fig4b, PA_META_Fig4c, PA_META_Fig4d, PA_META_Fig4e, PA_META_Fig4f, ncol=2)
dev.off()

png("PA_META_Fig5.png", units="in", width=10, height=10, res=200)
PA_META_Fig5
dev.off()

############
#tidy environment
rm(dfCI, fun_LB, fun_LM, newmods, pLB, pLM, predict_MA, EGG, estimate, ll95, ma_LB, ma_LM, my.formula, se.seq, T, t, ul95)
