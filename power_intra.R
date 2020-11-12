library(ggplot2)
library(tidyverse)
library(MASS)

#power_intra---------------
# mu1=0
# mu2=0.5
# sigma1=1
# sigma2=1
# cor<-c(0.5,0.2)
# es=mu2
power_intra<-function(N = 20, mu1, sigma1, mu2, sigma2,es=0.1,
                      cor = c(0, 0.2, 0.4, 0.6)){
mu_pre_post <- c(pre = mu1, post = mu2)
correlations <- cor
set.seed(1)
n_sims <- 100 # we want 100 simulations
p_vals <- c()
# this vector will contain the power for each sample-size (it needs the initial 0 for the while-loop to work)
cohens_ds <- c()

powers_at_cor <- list()
cohens_ds_at_cor <- list()

for(icor in 1:length(correlations)){ # do a power-simulation for each specified simulation
  n <- 2 # sample-size 
  i <- 2 # index of the while loop for saving things into the right place in the lists
  power_at_n <- c(0) 
  cohens_ds_at_n <- c() 
  sigma <- matrix(c(sigma1^2, sigma1*sigma2*correlations[icor], 
                    sigma1*sigma2*correlations[icor], sigma2^2), 
                  ncol = 2) #var-covar matrix
  while(power_at_n[i-1] < .95){
    for(sim in 1:n_sims){
      bivnorm <- data.frame(mvrnorm(n, mu_pre_post, sigma)) # simulate the bivariate normal
      p_vals[sim] <- t.test(bivnorm$pre, bivnorm$post, paired = TRUE, var.equal = TRUE, conf.level = 0.9)$p.value # run t-test and extract the p-value
      cohens_ds[sim] <- abs((mean(bivnorm$pre)-mean(bivnorm$post))/(sqrt(sd(bivnorm$pre)^2+sd(bivnorm$post)^2-2*cor(bivnorm$pre, bivnorm$post)*sd(bivnorm$pre)*sd(bivnorm$post)))) # we also save the cohens ds that we observed in each simulation
    }
    power_at_n[i] <- mean(p_vals < .10) # check power (i.e. proportion of p-values that are smaller than alpha-level of .10)
    names(power_at_n)[i] <- n
    cohens_ds_at_n[i] <- mean(cohens_ds) # calculate means of cohens ds for each sample-size
    names(cohens_ds_at_n)[i] <- n
    n <- n+1 # increase sample-size by 1
    i <- i+1 # increase index of the while-loop by 1 to save power and cohens d to vector
  }
  power_at_n <- power_at_n[-1] # delete first 0 from the vector
  cohens_ds_at_n <- cohens_ds_at_n[-1] # delete first NA from the vector
  powers_at_cor[[icor]] <- power_at_n # store the entire power curve for this correlation in a list
  cohens_ds_at_cor[[icor]] <- cohens_ds_at_n # do the same for cohens d
  names(powers_at_cor)[[icor]] <- correlations[icor] # name the power-curve in the list according to the tested correlation
  names(cohens_ds_at_cor)[[icor]] <- correlations[icor] # same for cohens d
  list1 <- mapply(cbind, powers_at_cor, "es"=es, SIMPLIFY=F)
  list2 <- mapply(cbind, list1, 
                  "cor"= as.numeric(names(powers_at_cor)), 
                  SIMPLIFY=F)
  list3 <- mapply(cbind, list2, "sd"=sigma1, SIMPLIFY=F)
}
return(list(power=list3))
}

power_intra_multiple<-function(sd){
  y<-c(-0.2, -0.4, -0.6, -0.8)  # different values of after-effects
  ls<-lapply(y, function(x) power_intra(N= 20, mu1 = 0, mu2 = x, 
                                        sigma1 = sd, sigma2 =sd, es = x))
  df<-bind_rows(ls)
  return(df)
}
s<-c(0.8,1,1.2, 1.4) # different values of sd
ls<-lapply(s, function(x) power_intra_multiple(sd=x))

fun <-function(y, x){
  y<-as.data.frame(ls[[y]]$power[x])
  y$id<-seq(1:nrow(y))
  colnames(y)<-c("power_at_n", "effect_size", "cor", "sd", "id")
  y$cor<-paste("corr = ",y$cor,sep="")
  return(y)
}
ls
ls1<-lapply(seq(1:9),function(x) fun(x=x,y=1))
ls2<-lapply(seq(1:length(ls[[2]]$power)),function(x) fun(x=x,y=2))
ls3<-lapply(seq(1:length(ls[[3]]$power)),function(x) fun(x=x,y=3))
ls4<-lapply(seq(1:length(ls[[4]]$power)),function(x) fun(x=x,y=4))

ls11<-bind_rows(ls1) %>% mutate(sd= "SD = 0.8")
ls22<-bind_rows(ls2)%>% mutate(sd= "SD = 1")
ls33<-bind_rows(ls3)%>% mutate(sd= "SD = 1.2")
ls44<-bind_rows(ls4)%>% mutate(sd= "SD = 1.4")

df_intra<-bind_rows(ls11, ls22,ls33,ls44)
# View(df_intra)

write.table(df_intra, file = "VV_SIMU_intra", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

p<-ggplot(data= df_intra[df_intra$sd== "SD = 1",], aes(x=id, y = power_at_n, 
                            col = as.factor(effect_size)))+
  # geom_point(size=1)+
  geom_hline(yintercept = .80)+
  scale_y_continuous(breaks=seq(0,1,.1))+
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
  scale_color_discrete(name = "LPA - SHAM (°)")+
  geom_smooth()+
  labs(x = "Number of participants",
       y = expression("Power ("*alpha~"=.10)"))+
  geom_vline(xintercept = c(30, 50, 80), lty = "dashed")+
  theme_bw(base_size = 20)+
  facet_wrap(~ cor)

 
png("VV_SIMU_INTRA.png", units="in", width=14, height=10, res=200)
p
dev.off()

# fusion-------------
df_inter<-read.table("VV_SIMU_inter.txt")
df_inter$cor<-paste("Between-S")
df_inter$design <-"between"
df_intra$design <-"within"

colnames(df_inter)
colnames(df_intra)
df_all<-bind_rows(df_inter, df_intra)


write.table(df_all, file = "VV_SIMU_all.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
df_plot<-df_all[(df_all$sd== "SD = 1")&(df_all$cor != "corr = 0"),]
p<-ggplot(data= df_plot, aes(x=id, y = power_at_n,col = as.factor(effect_size)))+
  # geom_point(size=1)+
  geom_hline(yintercept = .80)+
  scale_y_continuous(breaks=seq(0,1,.1))+
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
  scale_color_discrete(name = "LPA vs. SHAM (°)")+
  geom_smooth()+
  labs(x = "Number of participants",
       y = expression("Power ("*alpha~"=.10)"))+
  geom_vline(xintercept = c(30, 50, 70), lty = "dashed")+
  theme_bw(base_size = 20)+
  facet_wrap(~ cor)


png("VV_SIMU_ALL.png", units="in", width=14, height=10, res=200)
p
dev.off()
