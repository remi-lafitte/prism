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
  prior(normal(0, 2), coef = intercept),
  prior(cauchy(0, 1), class = sd)
)
# meta$STUDY
m1 <- brm(data = meta, family = gaussian,
            d | se(SD) ~ 1 + (1 | ID),
            prior = c(prior(normal(0, 1.5), class = Intercept),
                      prior(cauchy(0, 1), class = sd)),
  save_all_pars = TRUE,
  iter = 2000, warmup = 1000, chains = 4,seed = 14,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# simple model, without moderators, such as prism strenght, or duration.
# it will be the next step...
print(m1)
# good fit
m1 %>%
  plot(
    pars = c("^b_", "^sd_"),
    combo = c("dens_overlay", "trace"),
    theme = theme_bw(base_size = 16)
  )
print(m1)
# forest plot------------------
# load tidybayes
library(tidybayes)
get_variables(m1)

m1 %>%
  spread_draws(b_Intercept,r_ID[ID,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_Intercept + r_ID) %>%
  ungroup() %>%
  mutate(outcome = str_replace_all(ID, "[.]", " ")) %>% 
  # plot
  ggplot(aes(x = mu, y = reorder(outcome))) +
  geom_vline(xintercept = fixef(m1)[1, 1], color = "white", size = 1) +
  geom_vline(xintercept = fixef(m1)[1, 3:4], color = "white", linetype = 2) +
  geom_halfeyeh(.width = .95, size = 2/3) +
  labs(x = expression(italic("Cohen's d")),
       y = NULL) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0))+
  theme_bw(base_size=20)





# Study-specific effects are deviations + average
out_r <- m1 %>%
  spread_draws(b_Intercept,r_ID[ID,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_Intercept + r_ID)


# Average effect
out_f <- spread_draws(fit_rem, b_Intercept) %>% 
  mutate(study = "Average")
# Combine average and study-specific effects' data frames
out_all <- bind_rows(out_r, out_f) %>% 
  ungroup() %>%
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(study = fct_relevel(study, "Average"))
# Data frame of summary numbers
out_all_sum <- group_by(out_all, study) %>% 
  mean_qi(b_Intercept)
#> Warning: unnest() has a new interface. See ?unnest for details.
#> Try `cols = c(.lower, .upper)`, with `mutate()` needed
# Draw plot
out_all %>%   
  ggplot(aes(b_Intercept, study)) +
  geom_density_ridges(
    rel_min_height = 0.01, 
    col = NA,
    scale = 1
  ) +
  geom_pointintervalh(
    data = out_all_sum, size = 1
  ) +
  geom_text(
    data = mutate_if(out_all_sum, is.numeric, round, 2),
    # Use glue package to combine strings
    aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
    hjust = "inward"
  )
#> Picking joint bandwidth of 0.0226