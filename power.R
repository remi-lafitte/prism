library(ggplot2)
library(tidyverse)

power_inter<-function(N = 20, mu1, sigma1, mu2, sigma2,es=0.1){
set.seed(1)
n_sims <- 100 # we want 1000 simulations
p_vals <- c()
power_at_n <- c(0) # this vector will contain the power for each sample-size (it needs the initial 0 for the while-loop to work)
cohens_ds <- c()
cohens_ds_at_n <- c() 
n <- N # sample-size 
i <- 2
while(power_at_n[i-1] < .95){
  for(sim in 1:n_sims){
    group1 <- rnorm(n,mu1,sigma1) # simulate group 1
    group2 <- rnorm(n,mu2,sigma2) # simulate group 2
    p_vals[sim] <- t.test(group1, group2, paired = FALSE, var.equal = TRUE, conf.level = 0.9)$p.value # run t-test and extract the p-value
    cohens_ds[sim] <- abs((mean(group1)-mean(group2))/(sqrt((sd(group1)^2+sd(group2)^2)/2))) # we also save the cohens ds that we observed in each simulation
  }
  power_at_n[i] <- mean(p_vals < .10) # check power (i.e. proportion of p-values that are smaller than alpha-level of .10)
  cohens_ds_at_n[i] <- mean(cohens_ds) # calculate means of cohens ds for each sample-size
  n <- n+1 # increase sample-size by 1
  i <- i+1 # increase index of the while-loop by 1 to save power and cohens d to vector
}
power_at_n <- power_at_n[-1] # delete first 0 from the vector
cohens_ds_at_n <- cohens_ds_at_n[-1] # delete first NA from the vector

# plot(N:(n-1), power_at_n, xlab = "Number of participants per group", ylab = "Power", ylim = c(0,1), axes = TRUE)
# abline(h = .8, col = "red")
# abline(v = target, col = "green")

power_at_n<-as.data.frame(power_at_n)
power_at_n$id<-seq(N:(n-1))
power_at_n$effect_size<-es
return(power_at_n)
}
y<-c(-0.1,-0.2, -0.3, -0.4, -0.5,-0.7)
ls<-lapply(y, function(x) power_inter(N= 10, mu1 = 0, mu2 = x, 
              sigma1 = 0.8, sigma2 =0.8, es = x))
df<-bind_rows(ls)
ggplot(data= df, aes(x=id, y = power_at_n, 
                     col = as.factor(effect_size)))+
  geom_point(size=1)+
  geom_hline(yintercept = .80)+
  scale_y_continuous(breaks=seq(0,1,.1))+
  scale_x_continuous(limits = c(0,100), breaks=seq(0,100,10))+
  scale_color_discrete(name = "VV after-effect")+
  geom_smooth()+
  labs(x = "Number of participants",
       y = expression("Power ("*alpha~"=.10)"))+
  geom_vline(xintercept = c(10, 30, 40, 80), lty = "dashed")+
  theme_bw(base_size = 20)


