library(ggplot2)
library(gridExtra)
library(metafor)
library(ggpmisc)
library(tidyverse)
library(glue)
library(tidybayes)
library(magrittr)
library(bayesplot)
library(brms)
library(BEST)
library(broom.mixed)
library(ggridges)
library(brmstools)
#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv---------
#load the data
meta <- read.csv(here::here("McINTOSH_19_PA_META","PSA_META.csv"), sep = ";", dec = ",")
hist(meta$SD, breaks=10)
# setting the seed for reproducibility
set.seed(123)
meta$degz<-(meta$deg - mean(meta$deg))/sd(meta$deg)
#standardize degree = simplify the prior spec

cor(meta[c("prism","prism_c","duration_raw", "duration_c",
           "number","deg", "d")], method="p", use="complete.obs")

# defining the priors
prior1 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd)
)

range(meta$deg)
hist(meta$deg)
# hist(rcauchy(1e4, scale = 2, location = 3),xlim = c(- 200, 200), breaks = 10000)   

# PRIOR CHECK-----------
bmod0 <- brm(
  degz | se(SD) ~ 1 + (1|ID),
  data = meta,
  prior = prior1,
  sample_prior = "only",
  save_all_pars = TRUE,
  chains = 4,
  warmup = 5000,
  iter = 20000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# get_variables(bmod0)

mcmc_areas(
  as.array(bmod0), 
  pars = c("sd_ID__Intercept",
           "sd_ID__Intercept","Intercept"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)

bmod1 <- brm(
  deg | se(SD) ~ 1 + (1|ID),
  data = meta,
  prior = prior1,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

mcmc_areas(
  as.array(bmod1), 
  pars = c("sd_ID__Intercept",
           "sd_ID__Intercept","Intercept"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)

prior_summary(bmod1)
tidy(bmod1)
tidy(bmod1,effects = "ran_vals")

f<-forest(bmod1, show_data = T, sort = F)
f+
  scale_x_continuous(limits = c(-5,15), breaks = seq(-5,15,1))
(-0.23*sd(meta$deg))+mean(meta$deg)
 meta %>%
  group_by(ID) %>%
  summarise(Observed = mean(d) ) %>%
  mutate(Estimated = coef(bmod1)$ID[, , ] %>% data.frame %>% pull(Estimate)) %>%
  gather(type, Observed, Observed:Estimated) %>%
  ggplot(aes(x = ID, y = Observed, fill = type) ) +
  # geom_hline(yintercept = mean(post$b_Intercept), linetype = 2) +
  geom_point(pch = 21, size = 5, alpha = 0.8, colour = "white", show.legend = TRUE) +
  # scale_color_manual(values = rev(wes_palette(n = 2, name = "Chevalier1") ) )  +
  # scale_fill_manual(values = rev(wes_palette(n = 2, name = "Chevalier1") ) )  +
  theme_bw(base_size = 20) +
  scale_x_discrete(name=NULL)+
  ylab("D") +
  theme(legend.title = element_blank()) +
  coord_flip()

study.draws <- spread_draws(bmod1, r_ID[ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_ID + b_Intercept)
pooled.effect.draws <- spread_draws(bmod1, b_Intercept) %>% 
  mutate(ID = "Pooled Effect")
forest.data <- bind_rows(study.draws, pooled.effect.draws) %>% 
  ungroup() %>%
  mutate(ID = str_replace_all(ID, "[.]", " ")) %>% 
  mutate(ID = reorder(ID, b_Intercept))
forest.data.summary <- group_by(forest.data, ID) %>% 
  mean_qi(b_Intercept)

ggplot(aes(b_Intercept, relevel(ID, "Pooled Effect", after = Inf)), 
       data = forest.data) +
  geom_vline(xintercept = fixef(bmod1)[1, 1], color = "grey", size = 1) +
  geom_vline(xintercept = fixef(bmod1)[1, 3:4], color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  geom_density_ridges(fill = "grey", rel_min_height = 0.01, col = NA, scale = 1,
                      alpha = 0.8) +
  geom_pointintervalh(data = forest.data.summary, size = 1) +
  geom_text(data = mutate_if(forest.data.summary, is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf), hjust = "inward") +
  labs(x = "Standardized Mean Difference",
       y = element_blank()) +
  theme_minimal()

#more complex model------------
bmod1 <- brm(
  deg | se(SD) ~ 1 + prism_c + duration_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior1,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
