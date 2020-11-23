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
  filter(train == "prism") %>% 
  group_by(id, side, time) %>% 
  summarise(vv = mean(vv)) %>% 
ggplot(aes(x = side, y = vv,fill = time, col = time))+
  geom_point(width = 0.3, size=2, position = position_dodge(0.5))+
  geom_hline(yintercept = 0)+
  theme_bw(base_size=20)+
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 1, shape=21, col = "black",
               stroke= 1.5)
