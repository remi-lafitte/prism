#McIntosh, Brown & Young (2018)--------------
#"Meta-analysis of the visuospatial aftereffects of prism adaptation, with two novel experiments"
#ANALYSIS CODE, META_ANALYSIS
#LOAD REQUIRED LIBRARIES--------------
library(ggplot2)
library(gridExtra)
library(metafor)
library(ggpmisc)
library(tidyverse)
library(tidybayes)
library(magrittr)
library(brms)
library(BEST)
#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv---------
#load the data
meta <- read.csv(here::here("McINTOSH_19_PA_META","PSA_META.csv"), sep = ";", dec = ",")

#see the data structure
str(meta)
# View(meta)


#PRIOR PREDICTIVE CHECKING
sample_mu <-rnorm(1e4,0,10)
sample_sigma <-rcauchy(1000, 0, 1)
ppc<-data.frame(x = rnorm(1000, sample_mu, sample_sigma)) %>% 
  na.omit()

ggplot(data = ppc, aes(x) ) +
  scale_x_continuous(limits = c(-2,4), breaks = seq(-2,4,0.5) )+
  geom_histogram(col="white") +
  labs(x = "D", y = "Nombre d'échantillons") +
  theme_bw(base_size = 20)

# meta$d %>% range()
# median(meta$d)
# mean(meta$d)

BEST::plotPost(ppc$x, credMass = 0.95, cex = 1.5, 
               xlab = expression(theta), xlim = c(-4, 4) )

# meta$STUDY----
get_prior()
p <-c(
  prior(normal(0, 10), coef = intercept),
  prior(cauchy(0, 1), class = sd)
)
p2 <-c(
  prior(normal(0, 10), class = Intercept),
  prior(cauchy(0, 1), class = sd)
)
m <-brm(d | se(SD) ~ 0 + intercept + (1 | ID),
        family = gaussian(),data = meta,
        prior = p,
         save_all_pars = TRUE,
         iter = 2000, warmup = 1000, chains = 4,seed = 14,
         cores = parallel::detectCores(),
         control = list(adapt_delta = .95))

m2 <-brm(d | se(SD) ~ 1 + (1 | ID),
        family = gaussian(),data = meta,
        prior = p2,
        save_all_pars = TRUE,
        iter = 2000, warmup = 1000, chains = 4,seed = 14,
        cores = parallel::detectCores(),
        control = list(adapt_delta = .95))

mcmc_plot(m,type = "hist")
mcmc_plot(m, type = "rhat")
summary(m)
summary(m2)
pairs(m)


cor(meta[c("prism","prism_c","duration_raw", "duration_c",
           "number","deg", "d")], method="p", use="complete.obs")
# prism strenght
# duration raw
# number = missing obs +++
# meta$number


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
  scale_x_continuous(limits = c(0, 15), breaks = seq(from = 0, to=15, by=1)) +
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(from = 0, to=2.5, by=.5)) +
  labs(x = 'Raw effect size (cm)', y= 'Standardised effect size (d)') +
  geom_hline(yintercept=0, linetype="dotted") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major =element_blank(), panel.grid.minor =element_blank(),
        axis.title = element_text(size=16),
        plot.title = element_text(size=16), legend.position = c(0.23, .9),
        legend.key.size =  unit(0.1, "in"), legend.background = element_blank()) 


