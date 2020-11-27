


p_txt<-function(p){
  
  numformat <- function(val, digits) { 
    sub("^(-?)0.", "\\1.", sprintf(paste("%.",digits,"f",sep=""), val))
    }
  # function that removes the "0" before the decimal. Digits specify rounding
  a = 'p < '
  b = 'p = '
  p <- ifelse(is.character(p), as.numeric(p), p)
  p2 <- numformat(val = p, digits = 2)
  p3 <- numformat(val = p, digits = 3)
  
  pv<-
    ifelse(p < 0.000001,
           "*p* < 10^-6^",
           ifelse(p < .00001,
                  "*p* < 10^-5^",
                  ifelse(p < .0001,
                         "*p* < 10^-4^",
                         ifelse(p < .001,
                                "*p* < .001",
                                ifelse(p < .01,
                                       "*p* < .01",
                                       ifelse(p < .05,
                                       paste("*p* < ", p3, sep = ""),
                                       paste("*p* = ", p2, sep = "")))))))
  
  return(pv)
}
#function that returns p value with APA style
# EX  = p_txt(0.852)

p_txt_word<-function(p){
  
  numformat <- function(val, digits) { 
    sub("^(-?)0.", "\\1.", sprintf(paste("%.",digits,"f",sep=""), val))
  }
  # function that removes the "0" before the decimal. Digits specify rounding
  a = 'p < '
  b = 'p = '
  p <- ifelse(is.character(p), as.numeric(p), p)
  p2 <- numformat(val = p, digits = 2)
  p3 <- numformat(val = p, digits = 3)

  pv<-
    ifelse(p < 0.000001,
           "p < 10^-6^",
           ifelse(p < .00001,
                  "p < 10^-5^",
                  ifelse(p < .0001,
                         "p < 10^-4^",
                         ifelse(p < .001,
                                "p < .001",
                                ifelse(p < .01,
                                       "p < .01",
                                       ifelse(p < .05,
                                              paste("p < ", p3, sep = ""),
                                              paste("p = ", p2, sep = "")))))))
  
  return(pv)
}
#function that returns p value with APA style
# EX  = p_txt(0.852)

aov4_txt<- function(model){
  
  pes<-
    effectsize::eta_squared(model, partial = TRUE, ci = 0.9) %>%
    data.frame %>%
    mutate_if(is.numeric, round, digits=3)
  
  table<-
    model$anova_table %>%
    data.frame %>%
    mutate(Parameter = rownames(.)) %>% 
    group_by(Parameter) %>% 
    mutate_at(vars(Pr..F.), funs(p=p_txt(Pr..F.))) %>% 
    mutate_if(is.numeric, round, digits=2) %>%
    inner_join(., pes, by = "Parameter")
  
  table$Parameter <-gsub(":", "_", table$Parameter)
  
  
  txt<-
    table %>%
    mutate(F = paste("*F*(",num.Df, ", ",den.Df,") = ",F, sep = "")) %>%
    mutate(ges = paste("$\\hat{\\eta}^2_G$ = ",round(ges,2))) %>%
    mutate(pes = paste("$\\hat{\\eta}^2_p$ = ",round(Eta_Sq_partial,2))) %>%
    mutate(pes_ci =  paste("90% CI [",round(CI_low,2),", ",round(CI_high,2),"]", sep = "")) %>%
    mutate(pes_full =  paste(pes,", ", pes_ci, sep = "")) %>%
    mutate(full= paste(F,", ",p,", ", ges,", ", pes,", ", pes_ci, sep = "")) %>%
    mutate(small= paste(F,", ",p, sep = ""))
  
  rownames(txt)<- txt$Parameter
  list <- setNames(split(txt, seq(nrow(txt))), rownames(txt))
  return(list)
}
#function that returns ANOVA with APA style

tt_txt <- function(model, beta = T, unit = ""){
  ttest<-model
  # one sample t test
  
  q<-round(ttest$statistic,2)# statistic q
  dof<-round(ttest$parameter,2)# global degree of freedom
  pv<-p_txt(ttest$p.value)# raw p value. Beyond 10^-6, the round will give 0.
  b <- round(ttest$estimate[1] - ifelse(is.na(ttest$estimate[2]),0,ttest$estimate[2]),2)
  # estimation of the beta slope.
  
  full<- paste(
    ifelse(isTRUE(beta),
           paste("*M* = ",b,unit, ", ", sep = ""), ""),
    "95% CI [",round(ttest$conf.int[1],2),unit,", ",
    round(ttest$conf.int[2],2),unit,"], *t*(",dof,") = ",q,", ",pv,
    sep ="")
  
  small<- paste("*t*(",dof,") = ",q,", ",pv,
                sep ="")
  
  M<- paste("*M* = ",b,unit, sep ="")
  CI<- paste("95% CI [",round(ttest$conf.int[1],2),unit,", ",
    round(ttest$conf.int[2],2),unit,"],", sep ="")
  CI_raw <- paste("[",round(ttest$conf.int[1],2),unit,", ",
                  round(ttest$conf.int[2],2),unit,"],", sep ="")
  
  M_CI <-paste(M, CI, sep= "")
  
  
  
  return(list(full=full, small=small, M = M, CI = CI, M_CI = M_CI,
              M_raw = b, CI_raw = CI_raw, p = pv))
}
# function that returns a t-test with APA style

cor_txt <- function(model){
  cor<-model
  b<-round(cor$estimate,2) # estimation of the beta slope.
  # q<-round(cor$statistic,2)# statistic q
  dof<-round(cor$parameter,2)# global degree of freedom
  pv<-p_txt(cor$p.value)# raw p value. Beyond 10^-6, the round will give 0.
  
  full<- paste("*r*(",dof,") = ",b,
               ", 95% CI [",
               round(cor$conf.int[1],2),
               ", ",
               round(cor$conf.int[2],2),
               "], ",
               pv,"",
               sep ="")
  
  small<- paste("*r*(",dof,") = ",b,", ",pv,sep ="")
  
  
  return(list(full=full, small=small))
  
}
# function that returns a Pearson correlation with APA style

lm_txt<- function(model, unit = ""){

  library(broom)
  par<-broom::glance(model)
  txt<-
    broom::tidy(model, conf.int = T, ddl = T) %>%
    group_by(term) %>% 
    mutate(p = p_txt(p.value)) %>%
    mutate_if(is.numeric, round, digits=2) %>%
    mutate(t = paste("*t*(",par$df.residual,") = ",statistic, sep = "")) %>%
    mutate(slope = paste("$B$ = ", estimate,unit, sep = "")) %>%
    mutate(slope_ci= paste("95% CI [",conf.low,unit,", ",conf.high,unit,"]", sep = "")) %>%
    mutate(slope_full =  paste(slope,", ", slope_ci, sep = "")) %>%
    mutate(full= paste(t,", ",p,", ", slope_full, sep = "")) %>%
    mutate(small= paste(t,", ",p, sep = ""))
  
  rownames(txt)<- txt$term
  
  list <- setNames(split(txt, seq(nrow(txt))), rownames(txt))

  return(list)
}
# function that returns a Linear Model with APA style

khi_txt_2x2<-function(table){
  # https://stats.stackexchange.com/questions/188651/statistical-reporting-of-chi-square-and-odds-ratio
  # https://www.rdocumentation.org/packages/esc/versions/0.5.1/topics/esc_chisq
 library("esc")
  or<-esc_chisq(chisq=khi$statistic,totaln=sum(khi$observed),es.type="or")
  khi<-chisq.test(table)
  stati<-round(khi$statistic,2)
  pv<-khi$p.value
  or_ci<-paste("OR = ",round(or$es, 2),", 95% CI[",
                      round(or$ci.lo,2), ", ", round(or$ci.hi,2),
                      "]", sep = "")
  ddl<-khi$parameter

  inline<-paste("$\\chi^2$","(",ddl,", N = ",sum(khi$observed),") = ",stati,", ", p_txt(pv),", ",  or_ci, sep ="")
  return(inline)
}
# function that returns a 2 x 2 Chi Square with APA style


khi_txt<-function(table){
  # https://stats.stackexchange.com/questions/188651/statistical-reporting-of-chi-square-and-odds-ratio
  # https://www.rdocumentation.org/packages/esc/versions/0.5.1/topics/esc_chisq
  khi<-chisq.test(table)
  stati<-round(khi$statistic,2)
  pv<-khi$p.value
  ddl<-khi$parameter
  inline<-paste("$\\chi^2$","(",ddl,", N = ",sum(khi$observed),") = ",stati,", ", p_txt(pv), sep ="")
  return(inline)
}
# function that returns a Chi Square with APA style

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

q1_q3 <-function(data){  
  q1_q3<-quantile(data, probs = c(0.25, 0.75), na.rm = TRUE) 
  q1_q3<-paste("[",round(q1_q3[1],2), ", ", round(q1_q3[2],2), "]", sep = "")
  return(q1_q3)
}
# function that returns "[Q1, Q3]
# EX  = q1_q3(all$ben_angle_line1)

q2_q1_q3 <-function(data){  
  q2_q1_q3<-quantile(data, probs = c(0.25,0.5, 0.75), na.rm = TRUE) 
  q2_q1_q3<-paste(round(q2_q1_q3[2],2)," [",round(q2_q1_q3[1],2), ", ", round(q2_q1_q3[3],2), "]", sep = "")
  return(q2_q1_q3)
}
# function that returns "Q2 [Q1, Q3]
# EX  = q2_q1_q3(all$ben_angle_line1)

percent <-function(number){ 
  percentage<-paste(round(number*100,0),"%",sep = "")
  return(percentage)
}
# function that returns a percentage from a proportion
# EX  = percent(7/10)

table_wilcox <- function(iv, dv, data){
  
  data<-data[,c(iv,dv)]
  iv_name<-iv
  dv_name<-dv
  iv<- data[,iv]
  dv<- data[,dv]
  
  formula<-as.formula(paste(dv_name, "~", iv_name))
  
  W<-wilcox.test(formula = formula,data= data, conf.int = T, conf.level = 0.95)
  P<-p_txt(W$p.value)
  
 W_CI<-paste("95% CI [", W$conf.int[1] %>% round(.,2),
            ", ",
            W$conf.int[2] %>% round(.,2),
            "]", sep = "")
 # W_CI_raw<-paste("[", W$conf.int[1] %>% round(.,2),
 #           ", ",
 #           W$conf.int[2] %>% round(.,2),
 #           "]", sep = "")
  
  # Mdiff <- round(W$estimate,2)
  # tapply(data$vd, data$vi, median, na.rm=T)
 
   Z_P<-
    paste("Z = ",round(qnorm(W$p.value/2),2),", ", 
          P, sep="")
   Z<-round(qnorm(W$p.value/2),2)
   
  # R<-round(Z/sqrt(nrow(data)),2)
  library(rcompanion)
  
  ES<-wilcoxonR(x = dv,
            g = iv,ci=T)
  
  # RCI<- paste("r = ", round(ES[1],2),
  #   ", 95% CI [", round(ES[2],2), ", ",
  #   round(ES[3],2), "]", sep = "")
  
  R<-paste("r = ", round(ES[1],2),sep ="")
  # R_raw<-  round(ES[1],2)
  CI<-paste("95% CI [", round(ES[2],2),", ",round(ES[3],2), "]", sep = "")
  # CI_raw<-paste("[", round(ES[2],2),", ",
  #               round(ES[3],2), "]", sep = "")
  
  Levels<-levels(as.factor(iv))
  Level1 <- q2_q1_q3(data[data[,iv_name] == Levels[1],dv_name])
  Level2 <- q2_q1_q3(data[data[,iv_name] == Levels[2],dv_name])


  table<-data.frame(DV = dv_name,group_1 = Level1, group_2 = Level2, 
                 Wilcox = Z_P,
                 ES = R, CI = CI ,Z)
  return(table)
}

table_chi <- function(iv, dv, data){
  
  data<-data[,c(iv,dv)]
  iv_name<-iv
  dv_name<-dv
  iv<- data[,iv]
  dv<- data[,dv]
  table<-table(iv, dv)
  
  K<-chisq.test(table)
  P<-p_txt(K$p.value)

  Result <- khi_txt(table)
  
  library(rcompanion)
  ES<-round(cramerV(table,conf=0.95, ci = T, digits = 2),2)
  V<-paste("V = ",ES[1], sep ="")
  V_CI<-paste("V = ",ES[1], ", 95% CI [",ES[2],", ",ES[3],"]", 
             sep ="")
  CI<-paste("95% CI [",ES[2],", ",ES[3],"]", sep ="")

  table<-data.frame(DV = dv_name,group_1 = "", group_2 = "", 
                    KHI = Result,
                    ES = V, CI = CI,  K = K$statistic)

  return(table)
}


qqplot <- function(x, data, ...) {
  c(qqnorm(data[[x]], main = names(data)[x], ...),
    qqline(data[[x]])
  )
}

histplot <- function(x, data, ...) {
  hist(data[[x]], main = names(data)[x], ...)
}
out <- function(model, data, cook = 4, hat=3){
  library(ggfortify)
  library(broom)
  library(magrittr)
  library(dplyr)
  library(patchwork)

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
