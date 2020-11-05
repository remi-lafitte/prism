#SIMULATION OF DATA----------------------
library(tibble)
library(ggplot2)
library(tidyverse)
d<-read.table(here::here("SIMULATION", "psa_2019.txt"), sep=",", dec=",",
           header=T)
d %>% 
  filter(time == "pre", hand== "right") %>% 
  summarise(sd(bias))
levels(d2$group) <-c("sham", "lpa")
# SD = 4.33
# FUNCTION----------
msd <- function(df, ivs = "", dv, unit=""){
  
  ivs<-syms(ivs)
  
  msd<-
    df %>%
    select(!!!ivs, !!sym(dv)) %>%
    group_by(!!! ivs) %>%
    dplyr::summarise(mean = mean(!!sym(dv), na.rm = TRUE),
                     sd = sd(!!sym(dv), na.rm = TRUE),
                     lower.sd = mean - sd,
                     upper.sd = mean + sd,
                     n = n()) %>%
    mutate(se = sd / sqrt(n),
           lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
           upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se) %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    mutate(txt = paste("*M* = ",mean,unit,", *SD* = ",sd,
                       unit,sep = ""))
  
  
  
  return(msd)
}
#function that returns mean, sd, ci, and (M, SD) with APA style
#DATA--------------
v_lpa<-seq(101,130,1)
v_sham<-seq(131,160,1)

simu<-function(id, mu=0, sigma=1){
psa<-as_tibble(rnorm(n=10, mean = mu, sd = sigma))
psa$id<-rep(id, 10)
return(psa)
}
# simu(id = 56)

lpa<-lapply(v_lpa, function (x) simu(id= x, mu = 3, sigma=3))
sham<-lapply(v_sham, function (x) simu(id= x, mu = 0, sigma=3))

df<-bind_rows(lpa, sham)

df_avg<-
  df %>% 
  mutate(group = c(rep("LPA",300), rep("SHAM", 300))) %>% 
  group_by(id,group) %>% 
  rename(ae = "value") %>% 
  summarise(ae = mean(ae)) %>% ungroup() %>% 
  mutate(group_c = ifelse(group == "SHAM", -0.5, 0.5))
  
#PLOT-----------------
msd<-msd(df = df_avg, ivs = "group", dv = "ae")
pd<-position_dodge(0.5)

ggplot(data = df_avg, aes(x = group_c, y = ae, group = 1))+
  geom_jitter(width =0.02)+
  geom_smooth(method = "lm", se=T)+
  geom_errorbar(data  = msd, 
                aes(y=mean,ymin = lower.ci, ymax = upper.ci), 
                width = 0.3, position = pd, size = 0.5,
                col = "black")+
  geom_hline(yintercept = 0, lty = "dashed")+
  theme_bw(base_size = 20)+
  

fig1<- ggplot(data = msd, aes(x = group, y = mean, 
                       col = group, fill = group)) +
  labs(y = "After-effect (°)")+
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_jitter(data = df_avg, 
              aes(group, y = ae, col = group, fill = group),
              size = 3, 
              alpha = 0.3,
              position = position_jitterdodge(jitter.width = 0.05,  
                                              dodge.width = 0.5)
  ) +
  geom_errorbar(data  = msd, aes(ymin = lower.ci, ymax = upper.ci), 
                width = 0.1, position = position_dodge(0.1), size = 1.2,
                col = "black") +
  scale_x_discrete(name="Group",labels = c("SHAM", "LPA")) +
  geom_point(shape = 21,size = 4, 
             stroke = 1.3, 
             col = "black",
             position = position_dodge(0.1) )+
  theme_bw(base_size = 20)+
  guides(fill = F, color = F)
  
png("PA_SIMU_Fig1.png", units="in", width=10, height=10, res=200)
fig1
dev.off()

fig1

# p1<-c(2,4,4,2,2,1,1,2,1,2)
# p2<-c(3,0,-1,1.5,2,2.5,0,-1.5,0,-1)
# cbind(p1,p2)


