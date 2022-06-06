
####################################################################
# Subfunctions for IP weighted hazards model

# Parameters for choosing different scenarios/models
# treatment: 'smoking', 'sbp', 'cholesterol'
 

####################################################################
library(data.table)
library(survival)

 

# 1. Estimating treatment weight


getTreatW <- function(data_sub,treatment){
  
  
  # 1.1: Estimation of denominator
  
  
  
  if (treatment == 'smoking'){
    
    
    p.denom <- glm(treat ~ as.factor(view_ag_sex) +  
                     view_ag_age + 
                     asian + indian + maori + pacific +
                     view_en_nzdep +
                     as.factor(imp_hx_diabetes) +
                     as.factor(imp_hx_af) +
                     as.factor(pt_familyhistory) +
                     as.factor(imp_hx_lipidlowering) +
                     as.factor(imp_hx_antithrombotics) +
                     as.factor(imp_hx_antihypertensives) ,
                   data = data_sub,family = binomial(link= "logit"))
  }
  
  if  (treatment == 'sbp'){
    p.denom <- glm(treat ~ as.factor(view_ag_sex) +  
                     view_ag_age + 
                     asian + indian + maori + pacific +
                     view_en_nzdep +
                     ex_smoke +
                     cur_smoke +
                     as.factor(imp_hx_diabetes) +
                     as.factor(imp_hx_af) +
                     as.factor(pt_familyhistory) , 
                   data = data_sub,family = binomial(link= "logit"))
  }
  
  if  (treatment == 'cholesterol'){
    p.denom <- glm(treat ~ as.factor(view_ag_sex) +  
                     view_ag_age + 
                     asian + indian + maori + pacific +
                     view_en_nzdep +
                     ex_smoke +
                     cur_smoke +
                     as.factor(imp_hx_diabetes) +
                     as.factor(imp_hx_af) +
                     as.factor(pt_familyhistory) , 
                   data = data_sub,family = binomial(link= "logit"))
  }
  
  head(data_sub)
  data_sub$pd.treat <- predict(p.denom,data_sub,type="response") 
  rm(p.denom)
   
  # Estimation of numberator
  p.num <- glm(treat ~ 1,
               data = data_sub,family = binomial(link= "logit"))
  
  data_sub$pn.treat <- predict(p.num,data_sub,type="response")
  
  rm(p.num)
  
  # Computation of estimated weights
  data_sub$sw.a <- ifelse(data_sub$treat==1,
                          data_sub$pn.treat/data_sub$pd.treat,
                          (1-data_sub$pn.treat)/(1-data_sub$pd.treat))
  
  
  return(data_sub)
}


#########################################################
# 2. weight for loss to follow up and admin. censoring;
# treat censoring as time-varying tretment
# 
# Firstly, create person-month data for censoring
# install.packages("splitstackshape")
# library("splitstackshape")
# data.ipw <- expandRows(data_sub,"survtime",drop = F)
# library(tidyr) 

getCensW <- function(data_sub,treatment){
  # Replicate each sample by 'survtime' times
  idx <- rep(seq(1,dim(data_sub)[1]),data_sub$survtime)
  length(idx) == sum(data_sub$survtime)
  data.c.ipw <- data_sub[idx,]
  
  data.c.ipw$time <- sequence(rle(as.character(data.c.ipw$VSIMPLE_INDEX_MASTER))$lengths)
  head(data.c.ipw,61)
  
  data.c.ipw$event.c <- ifelse((data.c.ipw$time == data.c.ipw$survtime) 
                               & (data.c.ipw$censored==1),
                               1, 0)
  
  
  # Estimation of denominator
  
  if (treatment=='smoking'){
    p.c.denom <- glm(event.c==0 ~ as.factor(treat) + 
                       time  + 
                       as.factor(view_ag_sex) +  
                       view_ag_age + 
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       as.factor(imp_hx_diabetes) +
                       as.factor(imp_hx_af) +
                       as.factor(pt_familyhistory) +
                       as.factor(imp_hx_lipidlowering) +
                       as.factor(imp_hx_antithrombotics) +
                       sbp + 
                       imp_index2y_tchdl_ratio +
                       view_ag_age*as.factor(imp_hx_diabetes) + 
                       view_ag_age*sbp+
                       sbp*as.factor(imp_hx_antihypertensives) ,
                     data = data.c.ipw,family = binomial(link= "logit"))
    
  }
  
  
  
  if( (treatment == 'sbp'  ) |  (  treatment =='cholesterol' )){
    p.c.denom <- glm(event.c==0 ~ as.factor(treat) + 
                       time  + 
                       as.factor(view_ag_sex) +  
                       view_ag_age + 
                       asian + indian + maori + pacific +
                       view_en_nzdep+
                       ex_smoke +
                       cur_smoke  +
                       as.factor(imp_hx_diabetes) +
                       as.factor(imp_hx_af) +
                       as.factor(pt_familyhistory) +
                       as.factor(imp_hx_lipidlowering) +
                       as.factor(imp_hx_antithrombotics) +
                       as.factor(imp_hx_antihypertensives)  +
                       sbp + 
                       imp_index2y_tchdl_ratio +
                       view_ag_age*as.factor(imp_hx_diabetes) + 
                       view_ag_age*sbp+
                       sbp*as.factor(imp_hx_antihypertensives), 
                     data = data.c.ipw,family = binomial(link= "logit"))  
    
    
  }
  
  
  data.c.ipw$pd.csrk <- predict(p.c.denom,data.c.ipw,type="response")
  
  rm(p.c.denom) 
  
  
  
  ## Estimation of numerator for censoring weights
  
  p.c.num <- glm(event.c==0 ~ as.factor(treat) + time ,
                 data = data.c.ipw,family = binomial(link= "logit"))
  
  data.c.ipw$pn.csrk <- predict(p.c.num,data.c.ipw,type="response")
  
  rm(p.c.num)
  
  # Computation of estimated weights 
  
  data.c.ipw[event.c==0, sw.csrk:= pn.csrk/pd.csrk]
  data.c.ipw[event.c==1, sw.csrk:= (1-pn.csrk)/(1-pd.csrk)]
  data.c.ipw[, sw.c := prod(sw.csrk), by ='VSIMPLE_INDEX_MASTER']
  sw.c = data.c.ipw[time==1,.(VSIMPLE_INDEX_MASTER ,sw.c)]
  
  rm(data.c.ipw)
  dim(data_sub)
  dim(sw.c)
  sum(is.na(sw.c))
  
  
  data_sub <- merge(data_sub,sw.c,by='VSIMPLE_INDEX_MASTER')
  
  data_sub$sw <- data_sub$sw.a * data_sub$sw.c
  # sum(is.na(data_sub$sw))
  
  return(data_sub)
}



########################################
# 3. Compute Relative Risk (RR) from hazard model


getRR <- function(data,indices){

  # 3.1 Creation of person-month data for effect model; note diff from data.c.ipw
  
  data_sub <- data[indices,]
  idx <- rep(seq(1,dim(data_sub)[1]),data_sub$survtime)
  data.ipw <- data_sub[idx,]
  data.ipw$time <- sequence(rle(as.character(data.ipw$VSIMPLE_INDEX_MASTER))$lengths)-1
  data.ipw$event <- ifelse((data.ipw$time == data.ipw$survtime-1) & data.ipw$cvd==1,
                           1, 0)
  data.ipw$timesq <- data.ipw$time^2
  
  
  ## fit of weighted hazards model
  ipw.model.linint <- glm(event==1 ~ treat +  time + I(treat*time),
                          family = binomial(),
                          weight = sw, 
                          data = data.ipw)
  
 
  
  maxSur = max(data.ipw$survtime)
  
  # Creation of survival curves
  ipw.treat0 <- data.frame(cbind(seq(0,maxSur),0,(seq(0,maxSur))^2))
  ipw.treat1 <- data.frame(cbind(seq(0,maxSur),1,(seq(0,maxSur))^2))
  
  colnames(ipw.treat0) <- c("time","treat","timesq")
  colnames(ipw.treat1) <- c("time","treat","timesq")
   
  
  ipw.treat0$p.noevent0 <- 1-predict(ipw.model.linint,ipw.treat0,type = "response")
  ipw.treat1$p.noevent1 <- 1- predict(ipw.model.linint,ipw.treat1,type = "response")
  
  
  # Computation of survival for each person-month 
  ipw.treat0$surv0 <- cumprod(ipw.treat0$p.noevent0)
  ipw.treat1$surv1 <- cumprod(ipw.treat1$p.noevent1)
  
  
  ipw.graph <- merge(ipw.treat0, ipw.treat1, by =c('time','timesq'))
  ipw.graph$RR = (1-ipw.graph$surv1)/(1-ipw.graph$surv0)
  
  RR60 = ipw.graph[ipw.graph$time==60,'RR']
  
  return(RR60)
  
  
  
  
  # write.csv(ipw.graph[order(as.numeric(ipw.graph$time)),],RRoutputfname, row.names = FALSE)
  
  # ipw.graph$survdiffnw <- ipw.graph$surv1nw-ipw.graph$surv0nw
}





      