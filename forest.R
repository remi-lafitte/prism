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
meta$degz<-(meta$deg - mean(meta$deg))/sd(meta$deg)
#CORRELATION---------
# identify moderators = prism degree and duration
cor(meta[c("prism","prism_c","duration_raw", "duration_c",
           "number","deg", "d")], method="p", use="complete.obs")

# defining the priors--------------
prior1 <- c(
  prior(normal(4, 2), class = Intercept), # degree' prior
  prior(cauchy(0, 2), class = sd)
  # cauchy distribution
  # hist(rcauchy(1e4, scale = 2, location = 3),xlim = c(- 200, 200), breaks = 10000) 
)

# PRIOR CHECKING-----------
# model with prior only--------------
bmod0 <- brm(
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
bmod2 <- brm(
  deg | se(SD) ~ 1 + (1|STUDY),
  data = meta,
  prior = prior1,
  sample_prior = FALSE,
  save_all_pars = TRUE,
  chains = 4,
  warmup = 3000,
  iter = 12000,
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

# forest plot------------------
f<-
  forest(bmod1, grouping = "STUDY",
         fill_ridge = "dodgerblue", show_data = T, sort =F)+
  xlab("Effect size (Degree)") 
f2<-f+
  theme_bw()+
  scale_x_continuous(limits = c(-4,20), breaks = seq(-4,20,2))+
  geom_vline(xintercept = 0, lty = "dashed")+
  labs(y = "Article", x = "Overall effect size (°)")
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

# shrinkage----------
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
# other plot scripts------------------
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

#more complex models------------
# we want to assess the effect of 
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
#forest fx--------
# model <- bmod1_forest
# 
# grouping = "article"; pars = NA; level = 0.95; av_name = "Average"; 
# sort = FALSE; show_data = TRUE; col_ridge = NA; fill_ridge = "grey75"; 
# density = TRUE; text = TRUE; rel_min_height = 0.01; scale = 0.9;
# digits = 2; theme_forest = FALSE;

forest_ll <- function(
  model, grouping = NA, pars = NA, level = 0.95, av_name = "Average", 
  sort = TRUE, show_data = FALSE, col_ridge = NA, fill_ridge = "grey75", 
  density = TRUE, text = TRUE, rel_min_height = 0.01, scale = 0.9,
  digits = 2, theme_forest = TRUE
) {
  
  if (!requireNamespace("ggridges", quietly = TRUE) ) {
    
    stop("ggridges package needed for this function. Please install it.", call. = FALSE)
    
  }
  
  grouping <- brmstools:::get_grouping(model, grouping)
  probs <- c(0.5 - level / 2, 0.5 + level / 2)
  
  lwr <- paste0(probs[1] * 100, "%ile")
  upr <- paste0(probs[2] * 100, "%ile")
  
  samples <- tidycoef(model, pars = pars, grouping = grouping)
  
  samples_sum <- tidycoef(
    model, pars = pars, grouping = grouping, 
    summary = TRUE, level = level
  )
  
  samples[[grouping]] <- ifelse(is.na(samples[[grouping]]), av_name, samples[[grouping]])
  
  samples_sum[[grouping]] <- ifelse(is.na(samples_sum[[grouping]]), 
                                    av_name, samples_sum[[grouping]])
  
  samples_sum[["Interval"]] <- paste0(round(samples_sum[["Estimate"]], 
                                            digits), " [", round(samples_sum[[lwr]], digits), ", ", 
                                      round(samples_sum[[upr]], digits), "]")
  
  if (sort) {
    
    samples_sum <- dplyr::arrange_(samples_sum, "type", "Parameter", "Estimate")
    
  } else {
    
    samples_sum <- dplyr::arrange(samples_sum, type, desc(article) )
    
  }
  
  samples_sum[["order"]] <- forcats::fct_inorder(paste0(samples_sum[["type"]], 
                                                        samples_sum[[grouping]], samples_sum[["Parameter"]]) )
  
  samples <- dplyr::left_join(samples, samples_sum[, c(grouping, 
                                                       "Parameter", "order")], by = c(grouping, "Parameter") )
  
  g <- ggplot(samples_sum, aes_string("Estimate", "order") ) + 
    scale_y_discrete(labels = samples_sum[[grouping]], breaks = samples_sum[["order"]]) +
    geom_vline(
      xintercept = samples_sum %>% filter(type == "b") %>% pull(Estimate),
      linetype = 3
    ) +
    geom_point()
  
  if (density) {
    
    g <- g +
      ggridges::geom_density_ridges(
        data = samples, aes_string(x = "value"),
        rel_min_height = rel_min_height, 
        scale = scale, col = col_ridge, fill = fill_ridge
      ) +
      geom_point()
    
  }
  
  g <- g +
    geom_segment(
      aes_(y = ~order, yend = ~order, x = as.name(lwr), xend = as.name(upr) )
    )
  
  if (text) {
    
    g <- g +
      geom_text(
        data = samples_sum[samples_sum[["type"]] == "b", ],
        aes_string(label = "Interval", x = "Inf"),
        hjust = "inward", fontface = "bold", size = 3
      ) + 
      geom_text(
        data = samples_sum[samples_sum[["type"]] == "r", ],
        aes_string(label = "Interval", x = "Inf"), 
        hjust = "inward", size = 3
      )
    
  }
  
  if (show_data & length(unique(samples_sum[["Parameter"]]) ) == 1) {
    
    tmp <- dplyr::left_join(model$data, samples_sum[, c(grouping, "order")])
    
    g <- g +
      geom_point(
        data = tmp, aes_string(
          attr(attr(model$data, "terms"), "term.labels")[1], "order"),
        shape = 8
      )
    
  }
  
  # g <- g + facet_wrap("Parameter", scales = "free", strip.position = "bottom")
  
  if (theme_forest) g <- g + theme_forest()
  
  g
  
}
