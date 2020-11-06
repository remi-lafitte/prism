#LIBRARY---------------
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
library("hrbrthemes") # for the forest plot

#MANUALLY SET WORKING DIRECTORY TO DIRECTORY CONTAINING DATAFILE PA_META.csv---------

#load the data
meta <- read.csv(here::here("PSA_META","PSA_META.csv"), sep = ";", dec = ",")

# setting the seed for reproducibility
set.seed(123)

#standardize degree = simplify the prior spec (not used)
# meta$degz<-(meta$deg - mean(meta$deg))/sd(meta$deg)
#CORRELATION---------
# identify moderators = prism degree and duration
cor(meta[c("prism_c", "duration_c","number"
           ,"deg", "d")], method="p", use="complete.obs")
str(meta)
# defining the priors--------------
prior1 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd)
  # cauchy distribution
  # hist(rcauchy(1e4, scale = 2, location = 3),xlim = c(- 200, 200), breaks = 10000) 
)

# PRIOR CHECKING-----------
# model with prior only--------------
bmod_prior <- brm(
  deg | se(SD) ~ 1 + (1+STUDY) + (1|ID),
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

# plot distribution of the parameters-----------
mcmc_areas(
  as.array(bmod_prior), 
  pars = c("sd_ID__Intercept",
           "sd_ID__Intercept","Intercept"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "mean"
) + ggplot2::labs(
  title = "Prior parameter distributions",
  subtitle = "with medians and 80% intervals"
)

# Zero model-------
prior0 <- c(
  prior(normal(4, 2), class = Intercept) # degree' prior
)
bmod0 <- brm(
  deg | se(SD) ~ 1,
  data = meta,
  prior = prior0,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

# model1 - no moderator---------------
bmod1 <- brm(
  deg | se(SD) ~ 1 + (1|STUDY) + (1|ID),
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

# mcmc_areas(
#   as.array(bmod1), 
#   pars = c("sd_ID__Intercept",
#            "sd_ID__Intercept","Intercept"),
#   prob = 0.8, # 80% intervals
#   prob_outer = 0.95, # 95%
#   point_est = "mean"
# ) + ggplot2::labs(
#   title = "Prior parameter distributions",
#   subtitle = "with medians and 80% intervals"
# )
# prior_summary(bmod1)
tidy(bmod1)
VarCorr(bmod1, robust = F)
# tidy(bmod1,effects = "ran_vals")

#leave-one-out-----------
article_names <- sort(unique(meta$STUDY) )
bmods <- setNames(vector("list", length(article_names) ), article_names)
for (i in seq_along(article_names) ) {
  print(article_names[i])
  subdata <- droplevels(subset(meta, STUDY != article_names[i]) )
  capture.output({bmods[[i]] <- update(bmod1, newdata = subdata)})
}

intercepts_LOO <- as.numeric(unlist(lapply(bmods, function(x) brms::fixef(x)[1]) ) )
range(intercepts_LOO)
# the estimate is robust

# forest plot------------------
f<-
  forest(bmod1, grouping = "STUDY",
         fill_ridge = "dodgerblue", show_data = T, sort =F)+
  xlab("Effect size (Degree)") 
f2<-f+
  theme_bw()+
  scale_x_continuous(limits = c(-4,20), breaks = seq(-4,20,2))+
  geom_vline(xintercept = 0, lty = "dashed")+
  labs(y = "Article", x = "Overall effect size (?)")
f2
png("PSA_META_forest.png", units="in", width=10, height=10, res=200)
f2
dev.off()
summary(bmod1)

#PREDICTION-----------
# get_variables(bmod1)
s<-posterior_samples(bmod1, pars = c("b_Intercept","sd_STUDY__Intercept"))
BEST::plotPost(s$sd_STUDY__Intercept, credMass = 0.95)
predictive_interval(bmod1, prob = 0.95)

# prism degree effect------------
# we want to assess the effect of  prism degres
prior2 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd),
  prior(normal(0, 10), class = b)) 
# slope of ... between small and big optic deviation = non informative prior
bmod2 <- brm(
  deg | se(SD) ~ 1 + prism_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior2,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# check mcmc chain
plot(bmod2, combo = c("dens_overlay", "trace"), 
     theme = theme_bw(base_size = 16))
# retrieving the posterior samples
post_prism <- posterior_samples(bmod2, pars = "^b_")

# BF (long)
# (bf_prism <- 
#   bayes_factor(bmod2, bmod1,
#                repetitions = 1e2, cores = parallel::detectCores()))

prism_small<- post_prism[,1] - (post_prism[,2]*0.5)
prism_big<- post_prism[,1] + (post_prism[,2]*0.5)

png("PSA_META_small_prism.png", units="in", width=10, height=10, res=200)
BEST::plotPost(prism_small)
dev.off()

png("PSA_META_big_prism.png", units="in", width=10, height=10, res=200)
BEST::plotPost(prism_big)
dev.off()

BEST::plotPost(prism_small)
sd(prism_small)
BEST::plotPost(prism_big)
sd(prism_big)
hdi <- (5.88-0.102)/2
hdi
(bmod_prism_est <- 
    tidy(bmod2,parameters = c("^b_", "^sd_"), prob = 0.95))


# There is no reliable evidence for difference between the two conditions
# (beta = 2.32, 95% CrI [-1.58, 6.1], BF01 = 0.83).
post_prism %>%
  sample_n(size = 1e2) %>%
  rownames_to_column("draw") %>%
  expand(nesting(draw, b_Intercept, b_prism_c), a = c(-0.5,0.5)) %>%
  mutate(d = b_Intercept + b_prism_c * a) %>%
  ggplot(aes(x = a, y = d) ) +
  geom_point(data = meta, aes(x = prism_c, y = deg), size = 2) +
  geom_line(aes(group = draw), color = "purple", size = 0.5, alpha = 0.5) +
  labs(x = "Prism strenght", y = "Degree") +
  theme_bw(base_size = 20)

# difference between small and big
# c1 <- (post_prism[, 1] + post_prism[, 2]) - (post_prism[, 1] + post_prism[, 3])

# prism * number of movements------------
meta$number_c<- (meta$number - mean(meta$number))/sd(meta$number)
# we remove the study of GUINET13, a potential outlier.
# meta2<-meta
# meta2 %<>% filter(STUDY != "Guinet (2013)")

prior3 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd),
  prior(normal(0, 10), class = b)) 
# slope of ... between small and big optic deviation = non informative prior

bmod3 <- brm(
  deg | se(SD) ~ 1 + prism_c*duration_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior3,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# check mcmc chain
plot(bmod3, combo = c("dens_overlay", "trace"), 
     theme = theme_bw(base_size = 12))
tidy(bmod3)

bmod4 <- brm(
  deg | se(SD) ~ 1 + duration_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior3,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# check mcmc chain
plot(bmod4, combo = c("dens_overlay", "trace"), 
     theme = theme_bw(base_size = 12))
tidy(bmod4)
# model comparison------------
# calcul du WAIC et ajout du WAIC à chaque modèle
bmod0 <-add_criterion(bmod0, "waic")
bmod1 <-add_criterion(bmod1, "waic")
bmod2 <-add_criterion(bmod2, "waic")
bmod3 <-add_criterion(bmod3, "waic")
bmod4<-add_criterion(bmod4, "waic")

# comparaison des WAIC de chaque modèle
w <-loo_compare(bmod0,bmod1,bmod2,bmod3,bmod4, criterion = "waic")
print(w, simplify = FALSE)
model_weights(bmod0,bmod1,bmod2,bmod3,bmod4, weights = "waic") %>% 
  round(digits = 3)
# this model comparison shows that the best model contains the estimates of
# prism strenght.
# we keep model 2 = 
post_prism <- posterior_samples(bmod2, pars = "^b_")
prism_small<- post_prism[,1] - (post_prism[,2]*0.5)
prism_big<- post_prism[,1] + (post_prism[,2]*0.5)

BEST::plotPost(prism_small)
sd(prism_small)
BEST::plotPost(prism_big)
sd(prism_big)
hdi <- (5.88-0.102)/2
hdi

# small forest---------
meta_small<-meta %>% filter(prism == "small")
b_small <- brm(
  deg | se(SD) ~ 1 + (1|STUDY) + (1|ID),
  data = meta_small,
  prior = prior1,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
tidy(b_small)
f_small<-
  forest(b_small, grouping = "STUDY",
         fill_ridge = "dodgerblue", show_data = T, sort =F)+
  xlab("Effect size (Degree)") +
  theme_bw()+
  scale_x_continuous(limits = c(-4,20), breaks = seq(-4,20,2))+
  geom_vline(xintercept = 0, lty = "dashed")+
  labs(y = "Article", x = "Overall effect size (°)")
f_small
png("PSA_META_small_forest.png", units="in", width=10, height=10, res=200)
f_small
dev.off()

# COHEN D------------
prior_d <- c(
  prior(normal(0, 5), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd),
  prior(normal(0, 10), class = b)) 
# slope of ... between small and big optic deviation = non informative prior
bmod_d <- brm(
  d | se(SD) ~ 1 + prism_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior_d,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
# check mcmc chain
plot(bmod_d, combo = c("dens_overlay", "trace"), 
     theme = theme_bw(base_size = 16))
# retrieving the posterior samples
post_prism <- posterior_samples(bmod_d, pars = "^b_")
tidy(bmod_d)
prism_small<- post_prism[,1] - (post_prism[,2]*0.5)
prism_big<- post_prism[,1] + (post_prism[,2]*0.5)

png("PSA_META_small_prism_d.png", units="in", width=10, height=10, res=200)
BEST::plotPost(prism_small)
dev.off()
sd(prism_small)
png("PSA_META_big_prism_d.png", units="in", width=10, height=10, res=200)
BEST::plotPost(prism_big)
dev.off()

#mod1----------
prior1 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd))
bmod1 <- brm(
  deg | se(SE) ~ 1 + (1|STUDY) + (1|ID),
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
tidy(bmod1)
f3<-
  forest(bmod1, grouping = "STUDY",
         fill_ridge = "dodgerblue", show_data = T, sort =F)+
  xlab("Effect size (Degree)") +
  theme_bw()+
  scale_x_continuous(limits = c(-4,20), breaks = seq(-4,20,2))+
  geom_vline(xintercept = 0, lty = "dashed")+
  labs(y = "Article", x = "Overall effect size (d)")
f3
png("PSA_META_SE.png", units="in", width=10, height=10, res=200)
f3
dev.off()

#mod2-----------------
prior2 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd),
  prior(normal(0, 10), class = b)) 
# slope of ... between small and big optic deviation = non informative prior
bmod2 <- brm(
  deg | se(SE) ~ 1 + prism_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior2,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

post_prism <- posterior_samples(bmod2, pars = "^b_")
prism_small<- post_prism[,1] - (post_prism[,2]*0.5)
prism_big<- post_prism[,1] + (post_prism[,2]*0.5)
BEST::plotPost(prism_small)
BEST::plotPost(prism_big)
sd(prism_small)

#mod3----------------
bmod3 <- brm(
  d | se(SE) ~ 1 + prism_c + (1|ID) + (1|STUDY),
  data = meta,
  prior = prior2,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 4000,
  iter = 16000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)
tidy(bmod3)
post_prism <- posterior_samples(bmod3, pars = "^b_")
prism_small<- post_prism[,1] - (post_prism[,2]*0.5)
BEST::plotPost(prism_small)
#prediction------------

conditions <- make_conditions(bmod3b, "prism_d")
plot(conditional_effects(
  bmod3b, effects = "number",conditions=conditions,
  spaghetti = F, nsamples = 1e2),points = TRUE)

coef<-fixef(bmod3b)
# intercept = small prism + 0 pointing !
post_prism_nb2 <- posterior_samples(bmod3b, pars = "^b_")
m_mpg = brm(
  mpg ~ hp * cyl,
  data = mtcars,
  
  file = "models/tidy-brms_m_mpg.rds"  # cache model (can be removed)
)


# meta2 %>%
#   ggplot(aes(x = number, y = deg)) +
#   geom_abline(intercept = fixef(bmod3b)[1], 
#               slope     = fixef(bmod3b)[3]) +
#   geom_point(shape = 1, size = 2, color = "royalblue") +
#   theme_bw() +
#   theme(panel.grid = element_blank())
# mcmc_plot(bmod3b)
# 
# # we need new `nd` data
# nd <- expand.grid(number =  seq(0,350,20),
#                   prism_d = c(0,1),
#                   SD = 0)


# Data frame to evaluate average effects predictions on
# fit<-fitted(bmod3b, newdata = nd,re_formula = NA)
# pred<-predict(bmod3b, newdata = nd,re_formula = NA)
# 
# # re_formula is needed
# fitavg<-cbind(nd, fit) %>% 
#   filter(prism_d == 0) %>% select(-prism_d , -SD)
# 
# predavg<-cbind(nd, pred) %>% 
#   filter(prism_d == 0) %>% select(-prism_d , -SD,-Est.Error)
# 
# p1 <- ggplot(meta2, aes(x = number, y = deg)) +
#   geom_point(shape = 1) +
#   # geom_smooth(method = "lm", fill = "dodgerblue", level = .95) +
#   coord_cartesian(ylim = c(-15, 15))
# 
# colnames(fitavg)<-c("number", "deg", "SD","lower", "upper")
# colnames(predavg)<-c("number", "deg","lower", "upper")
# p1
# p1 +
#   geom_line(data = fitavg,group = "number", col = "black", size = 1) +
#   geom_line(data = fitavg, aes(y = lower), col = "black", lty = 2) +
#   geom_line(data = fitavg, aes(y = upper), col = "black", lty = 2) +
#   geom_ribbon(data = predavg, aes(ymin = lower, ymax = upper), 
#               col = "NA", alpha = .2,outline.type="full")+
#   theme_bw(base_size = 20)


coef[1,1] + (coef[3,1]*300)
# 4.2 degrees
post_prism_nb2 %>% 
  sample_n(size = 1e2) %>%
  rownames_to_column("draw") %>%
  expand(nesting(draw, b_Intercept, b_number), a = seq(50,400,10)) %>% 
  mutate(deg = b_Intercept + (b_number*a)) %>% 
  ggplot(aes(x = a, y = deg))+
  geom_line(aes(group = draw), color = "grey", 
            size = 0.5, alpha = 0.5)+
  geom_smooth(method = "lm", se =F)+
  theme_bw()

f <-
  fitted(b8.3, 
         newdata = nd,
         probs = c(.015, .985)) %>%
  data.frame() %>%
  bind_cols(nd) %>%
  mutate(cont_africa = ifelse(cid == 1, "African nations", "Non-African nations"))

post_prism_nb %>%
  sample_n(size = 1e2) %>%
  rownames_to_column("draw") %>%
  expand(nesting(draw, b_Intercept, b_number_c, b_prism_c)) %>%
  mutate(d1 = b_Intercept + (b_prism_c * -0.5) +  b_number_c) %>%
  mutate(d2 = b_Intercept + (b_prism_c * +0.5) +  b_number_c) %>%
  ggplot(aes(x = b_number_c, y = d1) ) +
  geom_point(data = meta, aes(x = number_c, y = deg), size = 2) +
  geom_line(aes(group = draw), color = "purple", size = 0.5, alpha = 0.5) +
  labs(x = "Number of movements (prism < 12.5°)", y = "PSA A-E (°)") +
  theme_bw(base_size = 20)

## shrinkage----------
meta %>%
  group_by(STUDY) %>%
  summarise(Observed = mean(deg) ) %>%
  mutate(Estimated = coef(bmod1)$ID[, , ] %>% data.frame %>% pull(Estimate)) %>%
  gather(type, Observed, Observed:Estimated) %>%
  ggplot(aes(x = STUDY, y = Observed, fill = type) ) +
  # geom_hline(yintercept = mean(post$b_Intercept), linetype = 2) +
  geom_point(pch = 21, size = 5, alpha = 0.8, colour = "white", show.legend = TRUE) +
  # scale_color_manual(values = rev(wes_palette(n = 2, name = "Chevalier1") ) )  +
  # scale_fill_manual(values = rev(wes_palette(n = 2, name = "Chevalier1") ) )  +
  theme_bw(base_size = 20) +
  scale_x_discrete(name=NULL)+
  ylab("D") +
  theme(legend.title = element_blank()) +
  coord_flip()

#other plot scripts------------------
study.draws <- spread_draws(bmod1, r_STUDY[ID,], b_Intercept) %>% 
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
