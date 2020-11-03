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
library(brms)
#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv---------
#load the data
meta <- read.csv(here::here("McINTOSH_19_PA_META","PSA_META.csv"), sep = ";", dec = ",")

#see the data structure
str(meta)
# View(meta)

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

prior <- c(
  prior(normal(0, 1.5), coef = intercept),
  # conservative prior
  prior(cauchy(0, 1), class = sd)
)
# meta$STUDY
m1 <- brm(data = meta, family = gaussian,
            d | se(SD) ~ 1 + (1 | ID) + (1 | STUDY),
            prior = c(prior(normal(0, 1.5), class = Intercept),
                      prior(cauchy(0, 1), class = sd)),
  save_all_pars = TRUE,
  iter = 2000, warmup = 1000, chains = 4,seed = 14,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
print(m1)
# good fit
m1 %>%
  plot(
    pars = c("^b_", "^sd_"),
    combo = c("dens_overlay", "trace"),
    theme = theme_bw(base_size = 16)
  )

# forest plot------------------
# load tidybayes
library(tidybayes)
get_variables(m1)

m1 %>%
  spread_draws(m1, b_Intercept[Intercept]) %>%
  head(10)

