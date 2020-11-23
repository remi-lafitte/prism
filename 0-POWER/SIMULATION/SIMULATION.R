#SIMULATION OF DATA----------------------
library(tibble)
library(ggplot2)
library(tidyverse)
library(MASS) 
library(tidyr)
library(magrittr)
library(purrr)
library(dplyr)
library(rethinking)
library(ggpubr)
select<- dplyr::select

# 1) Covariance between pre and post test--------------
old<-read.table(here::here("SIMULATION", "psa_2019.txt"),
                header = T, sep= ",", dec = ",")
d<-old %>%
  filter(group == "control") %>% 
  dplyr::select(time, subject, bias)  
pr<-d %>% filter(., time == "pre") %>% rename(pre = "bias") %>% select(-time)
po<-d %>% filter(., time == "post")%>% rename(post = "bias")%>% 
  select(-time, -subject)

df<-cbind(pr, po) %>% na.omit()
v<-seq(1:max(df$subject))
ls<-lapply(v, function(x) as.data.frame(cor(df[df$subject == x,]$pre, 
                                            df[df$subject == x,]$post)))
df2<-bind_rows(ls)
colnames(df2)<-"cor"
df2<-df2 %>% na.omit

mean(df2$cor)
sd(df2$cor)

plot(df2$cor)

# no evidence of covariance between pre test and post test.
df_t<-df %>% mutate(w = post-pre) %>%
  group_by(subject) %>% 
  summarise(m = mean(w))  
t.test(m ~ 1,df_t, mu = 0, y = NULL)

# no link between pre and post test.
# simplify data creation = no covariance
# covariance can range between 0 and 0.2 (cor)

# data frame simulation-----------------
## parameters-------
set.seed(5) # used to replicate example
pre <- 0 # average psa AE
post <- 3 # average psa AE
sigma <- 4 # std 
N <- seq(1:30)



sham_pre <- lapply(N, function (x) as.data.frame(rnorm(10 , pre , sigma)))
sham_post <- lapply(N, function (x) as.data.frame(rnorm(10 , pre , sigma)))
lpa_pre <- lapply(N, function (x) as.data.frame(rnorm(10 , pre, sigma)))
lpa_post <- lapply(N, function (x) as.data.frame(cbind(rnorm(10 ,post, sigma),x)))

spr<-bind_rows(sham_pre)
spo<-bind_rows(sham_post)
lpr<-bind_rows(lpa_pre)
lpo<-bind_rows(lpa_post)
x<-cbind(spr, spo, lpr, lpo)
colnames(x)<-c("s1", "s2", "l1", "l2", "id")
head(x)

# within subject variable = short format
wis <-x  
# within subject variable = long format
wil <-x %>% 
  gather(key = 'var', value = 'psa', -id) %>% 
  separate(var, sep = -1, into = c("group", "session")) 
wil_avg<-wil %>% 
  group_by(id, group, session) %>% 
  dplyr::summarise(psa=mean(psa))

#2)WITHIN SUBJECT approach------
library(afex)
m1<-afex::aov_4(psa ~ group*session + (1 + session*group | id), wil_avg)
m1
p<-ggplot(data = wil_avg, aes(x = group, y = psa,col = session))+
  geom_hline(yintercept=0, lty ="dashed")+
  geom_boxplot(size=0.2,position = position_dodge(0.7),
               width=0.4)+
  # stat_summary(fun= "mean", geom = "point",
  #              position = position_dodge(0.7), size = 6)
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 0.3)+
  scale_x_discrete(name = "Group", labels =c("LPA", "Control"))+
  scale_color_discrete(name = "Session", labels =c("Pre-test", "Post-test"))+
  labs(y="PSA (degree)")+
  theme_bw(base_size = 20)
  
png("PSA_SIMU_WI.png", units="in", width=10, height=10, res=200)
p
dev.off()

#3)BETWEEN SUBJECT APPROACH-------------
N1<-seq(100,130,1)
N2<-seq(131,160,1)
sham_pre <- lapply(N1, function (x) as.data.frame(rnorm(10 , pre , sigma)))
sham_post <- lapply(N1, function (x) as.data.frame(cbind(rnorm(10 , pre , sigma),x)))
lpa_pre <- lapply(N2, function (x) as.data.frame(rnorm(10 , pre, sigma)))
lpa_post <- lapply(N2, function (x) as.data.frame(cbind(rnorm(10 ,post, sigma),x)))

spr<-bind_rows(sham_pre)
spo<-bind_rows(sham_post)
lpr<-bind_rows(lpa_pre)
lpo<-bind_rows(lpa_post)
x1<-cbind(spr, spo)
x2<-cbind(lpr, lpo)
colnames(x1)<-c( "post", "pre",  "id")
colnames(x2)<-c("post","pre",  "id")
x11<-x1 %>% gather(key = "session", value = "psa",-id) %>% mutate(group="sham")
x22<-x2 %>% gather(key = "session", value = "psa",-id) %>% mutate(group="lpa")
x3<-bind_rows(x11,x22)

# between subject variable = long format
bwl <-x3  

bw_avg<-bwl %>% 
  group_by(id, group, session) %>% 
  dplyr::summarise(psa=mean(psa))
library(afex)
m2<-afex::aov_4(psa ~ group*session + (1 + session | id), bw_avg)
m2

p2<-ggplot(data = bw_avg, aes(x = group, y = psa,col = session))+
  geom_hline(yintercept=0, lty ="dashed")+
  geom_boxplot(size=0.2,position = position_dodge(0.7),
               width=0.4)+
  # stat_summary(fun= "mean", geom = "point",
  #              position = position_dodge(0.7), size = 6)
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 0.3)+
  scale_x_discrete(name = "Group", labels =c("LPA", "Control"))+
  scale_color_discrete(name = "Session", labels =c("Pre-test", "Post-test"))+
  labs(y="PSA (degree)")+
  theme_bw(base_size = 20)

png("PSA_SIMU_BW.png", units="in", width=10, height=10, res=200)
p2
dev.off()

#COV----------
mu <- c(0,0,0,3) 
stddev <- c(rep(4,4))
corMat <- matrix(c(1, 0.3, 0.3,0.3,
                   0.3, 1, 0.3,0.3,
                   0.3, 0.3, 1,0.3,
                   0.3, 0.3, 0.3,1),
                 ncol = 4)
corMat
covMat <- stddev %*% t(stddev) * corMat
covMat
ls_cov <- 
  lapply(N, function (x)
  as.data.frame(cbind(MASS::mvrnorm(n = 10, mu = mu, 
                              Sigma = covMat, empirical = FALSE),x)))
d_cov<-bind_rows(ls_cov)
  
  
# further form
colnames(d_cov)<-c("s1", "s2", "l1", "l2", "id")
# within subject variable = short format
wis <-d_cov
# within subject variable = long format
wil <-d_cov %>% 
  gather(key = 'var', value = 'psa', -id) %>% 
  separate(var, sep = -1, into = c("group", "session")) 
wil_avg<-wil %>% 
  group_by(id, group, session) %>% 
  dplyr::summarise(psa=mean(psa))

library(afex)
m3<-afex::aov_4(psa ~ group*session + (1 + session*group | id), wil_avg)
m3
p3<-ggplot(data = wil_avg, aes(x = group, y = psa,col = session))+
  geom_hline(yintercept=0, lty ="dashed")+
  geom_boxplot(size=0.2,position = position_dodge(0.7),
               width=0.4)+
  # stat_summary(fun= "mean", geom = "point",
  #              position = position_dodge(0.7), size = 6)
  stat_summary(fun.y = mean,
               geom = "pointrange",
               fun.ymax = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.ymin = function(x) mean(x) - sd(x) / sqrt(length(x)),
               position = position_dodge(0.7), size = 0.3)+
  scale_x_discrete(name = "Group", labels =c("LPA", "Control"))+
  scale_color_discrete(name = "Session", labels =c("Pre-test", "Post-test"))+
  labs(y="PSA (degree)")+
  theme_bw(base_size = 20)

png("PSA_SIMU_WI.png", units="in", width=10, height=10, res=200)
p3
dev.off()