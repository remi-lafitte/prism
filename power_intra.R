library(ggplot2)
library(tidyverse)
library(MASS)

#power_intra---------------


power_intra<-function(N = 20, mu1, sigma1, mu2, sigma2,es=0.1,
                      cor = c(0.2, 0.4, 0.6)){
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
  list2 <- mapply(cbind, list1, "cor"=cor, SIMPLIFY=F)
  list3 <- mapply(cbind, list2, "sd"=sigma1, SIMPLIFY=F)
  
  list4<-
    lapply(seq(1:3), function(x){
  c1<-as.data.frame(list3[[x]])
  c1$id<-seq(1:nrow(c1))
  colnames(c1)<-c("power_at_n", "es", "cor", "sd", "id")
  return(c1)
  }
  )
  power_final<-bind_rows(list4)

  return(power_final)
}
}

s<-c(0.8,1, 1.2, 1.4) # different values of sd
power_intra_multiple<-function(sd){
  y<-c(-0.4, -0.5,-0.7)  # different values of after-effects
  ls<-lapply(y, function(x) power_intra(N= 15, mu1 = 0, mu2 = x, 
                                        sigma1 = sd, sigma2 =sd, es = x))
  df<-bind_rows(ls)
  return(df)
}

ls2<-lapply(s, function(x) power_intra_multiple(sd = x))
df_intra<-bind_rows(ls2)
write.table(df_intra, file = "VV_SIMU_intra", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


df_intra6<-df_intra[df_intra$cor == 0.6,]
df_intra4<-df_intra[df_intra$cor == 0.4,]
df_intra2<-df_intra[df_intra$cor == 0.2,]
 
 p6<-ggplot(data= df_intra4, aes(x=id, y = power_at_n, 
                                 col = as.factor(es)))+
   # geom_point(size=1)+
   geom_hline(yintercept = .80)+
   scale_y_continuous(breaks=seq(0,1,.1))+
   scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
   scale_color_discrete(name = "VV after-effect")+
   geom_smooth()+
   labs(x = "Number of participants",
        y = expression("Power ("*alpha~"=.10)"))+
   geom_vline(xintercept = c(30, 50, 80), lty = "dashed")+
   theme_bw(base_size = 20)+
   facet_wrap(~ sd)
 
 p4<-ggplot(data= df_intra4, aes(x=id, y = power_at_n, 
                              col = as.factor(es)))+
  # geom_point(size=1)+
  geom_hline(yintercept = .80)+
  scale_y_continuous(breaks=seq(0,1,.1))+
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
  scale_color_discrete(name = "VV after-effect")+
  geom_smooth()+
  labs(x = "Number of participants",
       y = expression("Power ("*alpha~"=.10)"))+
  geom_vline(xintercept = c(30, 50, 80), lty = "dashed")+
  theme_bw(base_size = 20)+
  facet_wrap(~ sd)
 
 p2<-ggplot(data= df_intra2, aes(x=id, y = power_at_n, 
                                 col = as.factor(es)))+
   # geom_point(size=1)+
   geom_hline(yintercept = .80)+
   scale_y_continuous(breaks=seq(0,1,.1))+
   scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
   scale_color_discrete(name = "VV after-effect")+
   geom_smooth()+
   labs(x = "Number of participants",
        y = expression("Power ("*alpha~"=.10)"))+
   geom_vline(xintercept = c(30, 50, 80), lty = "dashed")+
   theme_bw(base_size = 20)+
   facet_wrap(~ sd)
library(patchwork) 
png("VV_SIMU_INTRA.png", units="in", width=36, height=12, res=200)
p2+p4+p6
dev.off()
