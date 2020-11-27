out <- function(model, data, cook = 4, hat=3){
  library(ggfortify)
  library(broom)
  library(magrittr)
  library(dplyr)
  library(patchwork)
  
  # Outlier detection + Assumptions 
  #=========================================
  cutoff_cook <- cook/((nrow(data)-length(model$coefficients)-2))
  df_cook<-data[which(cooks.distance(model) > cutoff_cook),] %>% 
    as_tibble() %>% 
    mutate(tool = "cookd")
  
  cutoff_hat <- mean(hatvalues(model)) * hat
  df_hat<-data[which(hatvalues(model) > cutoff_hat),] %>% 
    as_tibble() %>% 
    mutate(tool = "hat")
  
  n = nrow(data) 
  cutoff_sdr = qt(1 - 0.05 / (2*n), (n - 4)) 
  df_sdr<-data[which(abs(rstudent(model)) > cutoff_sdr),] %>% 
    as_tibble() %>% 
    mutate(tool = "sdr")
  
# plot 1
  multi_plot<-autoplot(model)
  p1<-multi_plot[2]+
    theme_bw(base_size=10)
# plot 2
  p2<-multi_plot[1]+
    theme_bw(base_size=10)
# plot 3 = d of cook
  p3<-cooks.distance(model) %>% 
  as_tibble() %>% rownames_to_column("id") %>% 
  ggplot()+aes(x=value, y = id, label=id)+
  labs(y =NULL)+guides(y = "none")+
  geom_label(size = 4, col = "blue")+
  theme_bw(base_size = 10)+
  geom_vline(xintercept = cutoff_cook, 
             lty="dashed", col="red",size=1)+
  labs(x = "Cook D", y ="Participant")
# plot 4 = hat value
  p4<-hatvalues(model) %>% 
    as_tibble() %>% rownames_to_column("id") %>% 
    ggplot()+aes(x=value, y = id, label=id)+
    geom_label(size = 4, col = "blue")+
    labs(y =NULL)+guides(y = "none")+
    theme_bw(base_size = 10)+
    geom_vline(xintercept = cutoff_hat, lty="dashed", col="red",size=1)+
    labs(x = "Leverage", y ="Participant")
  # plot 5 = SDR
 p5<- rstudent(model) %>% 
    as_tibble() %>% rownames_to_column("id") %>% 
   ggplot()+aes(x=value, y = id, label=id)+
   geom_label(size = 4, col = "blue")+
   labs(y =NULL)+guides(y = "none")+
    theme_bw(base_size = 10)+
    geom_vline(xintercept = c(-cutoff_sdr,cutoff_sdr), lty="dashed", col="red",size=1)+
    labs(x = "R Student", y ="Participant")

  
#p_out <- recordPlot()
  outlier_name <- 
    bind_rows(df_cook, df_hat, df_sdr) %>% # list of outliers
    select(id) %>% unique()
  
#New dataframe without outlier------- 
#===========================================
data_no_outlier <- data %>% 
    dplyr::filter(!id %in% outlier_name)
old_model<-model
new_model <- update(model,data=data_no_outlier)
coeff<-bind_rows(tidy(old_model, conf.int = T) %>% mutate(model = "old"),
                 tidy(new_model,conf.int = T) %>% mutate(model = "new"))
  

  p6<-ggplot(data = coeff, 
         aes(x = term, y = estimate,  
         ymin = conf.low, ymax = conf.high, col = model)) +
    geom_point(size = 2, position = position_dodge(0.2)) +
    geom_errorbar(width = 0.01,  position = position_dodge(0.2)) +
    geom_hline(yintercept = 0, lty = "dashed",size=1,col="red") +
    coord_flip()+ 
    labs(x ="" , y = "Estimate")+ 
    scale_color_discrete(name = "Model", 
                         labels = c("New", "Old"))+
    theme_bw(base_size=10)
  plot<-p1+p2+p3+p4+p5+p6

  return(list(plot=plot,outlier = outlier_name, 
              new_data = data_no_outlier, coeff=coeff))
}
