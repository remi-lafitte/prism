# Simulation of data -------------------
library(truncnorm) 
library(tidyverse)
library(dplyr)

simu<-function(x){
# VV----------
# VV in pretest--------
vv_pre_prism_left <- rnorm(n=10, mean = 0, sd =1)
vv_pre_prism_control <- rnorm(n=10, mean = 0, sd =1)
vv_pre_vr_left <- rnorm(n=10, mean = 0, sd =1)
vv_pre_vr_control <- rnorm(n=10, mean = 0, sd =1)
# VV in posttest--------
vv_post_prism_left <- rnorm(n=10, mean = -0.5, sd =1)
vv_post_prism_control <- rnorm(n=10, mean = 0, sd =1)

vr_left<-rtruncnorm(n=1, a=3, b=17, mean=10, sd=7) # rnorm with lower and upper limits
vv_post_vr_left <- rnorm(n=10, mean = vr_left, sd =2)
vv_post_vr_control <- rnorm(n=10, mean = 0, sd =1)

vv_all<-cbind(vv_pre_prism_control, vv_pre_prism_left, vv_pre_vr_control, vv_pre_vr_left,
              vv_post_prism_control, vv_post_prism_left, vv_post_vr_control, vv_post_vr_left)

# SSA----------
# SSA in pretest--------

ssa_pre<-rtruncnorm(n=1, a=-5, b=5, mean=0, sd=3) # rnorm with lower and upper limits

ssa_pre_prism_left <- rnorm(n=10, mean = ssa_pre, sd =3)
ssa_pre_prism_control <- rnorm(n=10, mean = ssa_pre, sd =3)
ssa_pre_vr_left <- rnorm(n=10, mean = ssa_pre, sd =1)
ssa_pre_vr_control <- rnorm(n=10, mean = ssa_pre, sd =1)
# ssa in posttest--------
ssa_post<-rtruncnorm(n=1, a=-2, b=6, mean= 3, sd= 3) # rnorm with lower and upper limits

ssa_post_prism_left <- rnorm(n=10, mean = ssa_post, sd =1)
ssa_post_prism_control <- rnorm(n=10, mean = ssa_pre, sd =1)

ssa_post_vr<-rtruncnorm(n=1, a=-2, b=6, mean= 1, sd= 3) # rnorm with lower and upper limits
ssa_post_vr_left <- rnorm(n=10, mean = ssa_post_vr, sd =2)
ssa_post_vr_control <- rnorm(n=10, mean = ssa_pre, sd =1)

ssa_all<-cbind(ssa_pre_prism_control, ssa_pre_prism_left, ssa_pre_vr_control, ssa_pre_vr_left,
              ssa_post_prism_control, ssa_post_prism_left, ssa_post_vr_control, ssa_post_vr_left)

#Multi VV------------------
corMat <- matrix(c
            (1,0,0,0,0,0,0,0,
              0,1,0,0,0,0,0,0,
              0,0,1,0,0,0,0,0,
              0,0,0,1,0,0,0,0,
              0,0,0,0,1,0,0,0,
              0,0,0,0,0,1,0,0,
              0,0,0,0,0,0,1,0,
              0,0,0,0,0,0,0,1), 
            nrow = 8, ncol = 8)


mu_vv <- c(apply(vv_all, 2, mean)) 
stddev_vv <- c(apply(vv_all, 2, sd))
covMat_vv <- stddev_vv %*% t(stddev_vv) * corMat
df_vv<-as.data.frame(MASS::mvrnorm(n=10, mu = mu_vv, Sigma = covMat_vv))

mu_ssa <- c(apply(ssa_all, 2, mean)) 
stddev_ssa <- c(apply(ssa_all, 2, sd))
covMat_ssa <- stddev_ssa %*% t(stddev_ssa) * corMat
df_ssa<-as.data.frame(MASS::mvrnorm(n=10, mu = mu_ssa, Sigma = covMat_ssa))

df_vv_long<-
  df_vv %>% 
  gather(var, vv) %>%
  separate(var, sep ="_", into = c("measure", "time", "train", "side")) %>% 
  select(-measure) %>% 
  as_tibble()

df_ssa_long<-
  df_ssa %>% 
  gather(var, ssa) %>%
  separate(var, sep ="_", into = c("measure", "time", "train", "side")) %>% 
  select(-measure) %>% 
  as_tibble()

df_final<-df_vv_long %>% 
  inner_join(.,df_ssa_long, by = c("train", "time", "side")) %>% 
  mutate(id = x)


return(df_final)
}

 
ls<-lapply(seq(101,130,1), function(x) simu(x=x))
df_simu<-bind_rows(ls)
df_simu %>% 
  # filter(train == "prism") %>% 
  group_by(id, side,train, time) %>% 
  summarise(vv = mean(vv)) %>% 
ggplot(aes(x = side, y = vv,fill = time, col = time))+
  geom_point(width = 0.3, size=2,alpha=0.5,
             position = position_dodge(0.5))+
  geom_hline(yintercept = 0)+
  theme_bw(base_size=20)+
  facet_wrap(~ train, scale="free")+
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 1, shape=21, col = "black",
               stroke= 1.5)

p<-
  df_simu %>% 
  rename("Visual Vertical" = "vv", "Straight Ahead"=ssa) %>% 
  mutate(train = fct_recode(train, 
             "Prism" = "prism", 
             "Tilted Room" = "vr"),
         side = fct_recode(side, 
                           "Control" = "control", 
                           "Leftward" = "left")) %>% 
  gather(., key = "task", value = "degree", 
         c(`Straight Ahead`, `Visual Vertical`)) %>% 
  group_by(id, side,train, time, task) %>% 
  summarise(degree = mean(degree)) %>% 
  ggplot(aes(x = side, y = degree,fill = time, col = time))+
  geom_point(width = 0.3, size=2,alpha=0.5,
             position = position_dodge(0.5))+
  geom_hline(yintercept = 0)+
  theme_bw(base_size=20)+
  facet_wrap(~ train*task, scale="free")+
  labs(y = "Degree", x = "Deviation")+
  scale_color_discrete(name = "Time", labels =c("Post-test", "Pre-test"))+
  scale_fill_discrete(name = "Time", labels =c("Post-test", "Pre-test"))+
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 1, shape=21, col = "black",
               stroke= 1.5)
  
png("Expected_Data.png", units="in", width=14, height=10, res=200)
p
dev.off()
p

#STATISTICS----------
# T-TEST----------
# aggregate
df_simu_mean<-df_simu %>% 
  group_by(id, time, side, train) %>% 
  summarize(vv = mean(vv), ssa= mean(ssa)) %>% 
  ungroup()
# dcast = wide format
library(data.table)
setDT(df_simu_mean)
 df_simu_wide<-
  dcast(df_simu_mean, id ~ time+side+train,
                    value.var = c("ssa", "vv"))
df_t<-
  as.data.frame(df_simu_wide) %>% 
   select(id, contains("prism")) %>% 
   mutate(
     time =
     (vv_post_control_prism + vv_post_left_prism)-
     (vv_pre_control_prism + vv_pre_left_prism),
     side =
     (vv_post_left_prism + vv_pre_left_prism)-
   (vv_post_control_prism + vv_pre_control_prism),
      time_side = 
     (vv_post_left_prism + vv_pre_control_prism)-
     (vv_pre_left_prism + vv_post_control_prism))

lm(time_side ~ 1, df_t) %>% summary()
lm(time ~ 1, df_t) %>% summary()
lm(side ~ 1, df_t) %>% summary()

sort(rstudent(lm(time_side ~ 1, df_t)))
x<-out(lm(time_side ~ 1, df_t), df_t)
x$outlier
# xtabs(vv ~ id+time, data = df_simu_mean)
# aggregate(vv ~ id, df_simu_mean,I)

# ANOVA----------
library(afex)
m1<-afex::aov_4(vv ~ time*side + (time*side|id), 
                data =df_simu[df_simu$train=="prism",])
                                                              
summary(m1)
# le scénario avec shift VV de -0.5° et SD = 1 fonctionne pour train = prism!

m2<-afex::aov_4(ssa ~ time*side + (time*side|id), 
                data =df_simu[df_simu$train=="vr",])
summary(m2)
# le scénario avec shift PSA de +1° et SD = 3 ne 
# fonctionne pas pour train = VR !

# LMER--------------
# Interaction between the variables Deviation and Time
library(lme4)
library(lmerTest)
str(df_simu)
lmer(formula = vv ~ time*side + (time*side|id), data =df_simu)
# warning = model too complex
# BRMS---------------------
# library(tidybayes)
# library(brms)
# 
# df_simu$time_c<-ifelse(df_simu$time == "pre", -0.5, 0.5)
# df_simu$side_c<-ifelse(df_simu$side == "control", -0.5, 0.5)
# 
# prior_vv <- c(
#   prior(normal(0, 10), class = Intercept), 
#   # uninformative prior
#   prior(cauchy(0, 10), class = sd),
#   prior(normal(0, 10), class = b),
#   prior(lkj(2), class = cor)
# )
# # modelisation
# m1 <- brm(
#   vv ~ 1 + time_c*side_c + (1 + time_c*side_c|id),
#   data = df_simu,
#   prior = prior_vv,
#   sample_prior = TRUE,
#   save_all_pars = TRUE,
#   chains = 4,
#   warmup = 4000,
#   iter = 16000,
#   cores = parallel::detectCores(),
#   control = list(adapt_delta = .99)
# )
# result
# tidy(mod_degree)