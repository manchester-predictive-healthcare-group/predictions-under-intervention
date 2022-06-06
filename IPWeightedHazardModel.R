
####################################################################
# Subfunctions for IP weighted hazards model

# Parameters for choosing different scenarios/models
# treatment: 'smoking', 'sbp', 'cholesterol'
# treatType: 'level', 'medication'
# ifFullModel: TRUE - using full model in the Post-baseline model
#              FALSE -using intermediate model in the Post-baseline model 

# Corr. to scenarios in the paper:
# I: treatment = 'smoking', treatType = None
# IIa: treatment = 'sbp', treatType = 'level'
# IIb: treatment = 'sbp', treatType = 'medication'
# IIIa: treatment = 'cholesterol', 'level'
# IIIb: treatment = 'cholesterol', treatType = 'medication'

####################################################################

# 1. Estimating treatment weight

getTreatW <- function(data_sub, treatment, ifFullModel, treatType){

  treatment, ifFullModel, treatType  
  
  # 1.1: Estimation of denominator
  
  # smoking 
  if (treatment == 'smoking'){ # Scenario I
    if (ifFullModel == F){  # if using 'intermediate' model
      p.denom <- glm(treat ~ as.factor(view_ag_sex) +  
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       view_ag_age_BL + 
                       as.factor(imp_hx_diabetes_BL) +
                       as.factor(imp_hx_af_BL) +
                       as.factor(pt_familyhistory) ,
                     data = data_sub,family = binomial(link= "logit"))
    }else{  # if using 'full' model
      p.denom <- glm(treat ~ as.factor(view_ag_sex) +  
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       view_ag_age_BL + 
                       sbp_BL+
                       imp_index2y_tchdl_ratio_BL+
                       as.factor(pt_familyhistory) + 
                       as.factor(imp_hx_diabetes_BL)  +
                       as.factor(imp_hx_af_BL) +
                       as.factor(imp_hx_lipidlowering_BL) +
                       as.factor(imp_hx_antithrombotics_BL) ),
                     data = data_sub,family = binomial(link= "logit"))
      
    }
  }  # smoking 
  
  
  # SBP
  if  (treatment == 'sbp'){
    if (ifFullModel == F){  # if using 'intermediate' model 
      p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       view_ag_age_BL + 
                       ex_smoke_BL +
                       cur_smoke_BL +
                       as.factor(imp_hx_diabetes_BL) +
                       as.factor(imp_hx_af_BL) +
                       as.factor(pt_familyhistory)+
                       imp_index2y_tchdl_ratio_BL,
                     data = data_sub,
                     family = binomial(link= "logit"))
      
    }else {   # if using 'FULL' model
      # Scenario IIa
      if (treatType == 'level'){ 
        p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                         asian + indian + maori + pacific +
                         view_en_nzdep +
                         view_ag_age_BL + 
                         ex_smoke_BL +
                         cur_smoke_BL +  
                         imp_index2y_tchdl_ratio_BL+
                         as.factor(pt_familyhistory) + 
                         as.factor(imp_hx_diabetes_BL)  +
                         as.factor(imp_hx_af_BL) +
                         as.factor(imp_hx_lipidlowering_BL) +
                         as.factor(imp_hx_antithrombotics_BL) +
                         as.factor(imp_hx_antihypertensives_BL) +  
                         ),
                       data = data_sub, 
                       family = binomial(link= "logit"))
      }else if (treatType == 'medication'){  
        # Scenario IIb
        p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                         asian + indian + maori + pacific +
                         view_en_nzdep +
                         view_ag_age_BL + 
                         ex_smoke_BL +
                         cur_smoke_BL +  
                         imp_index2y_tchdl_ratio_BL+
                         as.factor(pt_familyhistory) + 
                         as.factor(imp_hx_diabetes_BL)  +
                         as.factor(imp_hx_af_BL) +
                         as.factor(imp_hx_lipidlowering_BL) +
                         as.factor(imp_hx_antithrombotics_BL)  ,
                       data = data_sub, 
                       family = binomial(link= "logit"))
        
        
        
    }  
  } # sbp
  }
  
  # cholesterol
  if  (treatment == 'cholesterol'){
    if (ifFullModel == F){
      p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       view_ag_age_BL + 
                       ex_smoke_BL +
                       cur_smoke_BL +
                       as.factor(imp_hx_diabetes_BL) +
                       as.factor(imp_hx_af_BL) +
                       as.factor(pt_familyhistory)+
                       sbp_BL,
                     data = data_sub,
                     family = binomial(link= "logit"))
    }else {   # if using 'FULL' model
      # Scenario IIIa
      if (treatType == 'level'){ 
        p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                       view_ag_age_BL + 
                       asian + indian + maori + pacific +
                       view_en_nzdep +
                       ex_smoke_BL +
                       cur_smoke_BL +
                       sbp_BL+
                       imp_index2y_tchdl_ratio_BL+
                       as.factor(pt_familyhistory) + 
                       as.factor(imp_hx_diabetes_BL)  +
                       as.factor(imp_hx_af_BL) +
                       as.factor(imp_hx_lipidlowering_BL) +
                       as.factor(imp_hx_antithrombotics_BL) +
                       as.factor(imp_hx_antihypertensives_BL) ,
                     data = data_sub,
                     family = binomial(link= "logit"))
      }else if (treatType == 'medication'){  
        # Scenario IIIb
        p.denom <- glm(treat==1 ~ as.factor(view_ag_sex) +  
                         view_ag_age_BL + 
                         asian + indian + maori + pacific +
                         view_en_nzdep +
                         ex_smoke_BL +
                         cur_smoke_BL +
                         sbp_BL+
                         imp_index2y_tchdl_ratio_BL+
                         as.factor(pt_familyhistory) + 
                         as.factor(imp_hx_diabetes_BL)  +
                         as.factor(imp_hx_af_BL) +
                         as.factor(imp_hx_antithrombotics_BL) +
                         as.factor(imp_hx_antihypertensives_BL) ,
                       data = data_sub,
                       family = binomial(link= "logit"))
    }  
  }
  
  }# cholesterol
  
  data_sub$pd.treat <- predict(p.denom,data_sub,type="response")
  
  rm(p.denom)
  
  # 1.2: Estimation of numberator
  
  p.num <- glm(treat==1 ~ 1,
               data = data_sub,
               #weights =  sw.select,
               family = binomial(link= "logit"))
  
  data_sub$pn.treat <- predict(p.num,data_sub,type="response")
  
  rm(p.num)
  
  
  # 1.3 Computation of estimated weights
  data_sub$sw.a <- ifelse(data_sub$treat==1,
                          data_sub$pn.treat/data_sub$pd.treat,
                          (1-data_sub$pn.treat)/(1-data_sub$pd.treat))
  
  
  return(data_sub)
}
# End of treatment weight computaion

#########################################################
# 2. Weight for loss to follow up and administrative censoring;
# Here we treated censoring as time-varying tretment

getCensW <- function(data_sub, treatment){
  
  
  # Firstly, create person-month data for censoring
  # Replicate each sample by 'survtime' times
  idx <- rep(seq(1,dim(data_sub)[1]),data_sub$survtime)
  length(idx) == sum(data_sub$survtime)
  data.c.ipw <- data_sub[idx,]
  
  data.c.ipw$time <- sequence(rle(as.character(data.c.ipw$VSIMPLE_INDEX_MASTER))$lengths)
  
  data.c.ipw$event.c <- ifelse((data.c.ipw$time == data.c.ipw$survtime) 
                               & (data.c.ipw$censored==1),
                               1, 0)
  # 2.1 Estimation of denominator
  
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
  
  # 2.2 Estimation of numerator for censoring weights
  
  p.c.num <- glm(event.c==0 ~ as.factor(treat) + time ,
                 data = data.c.ipw,family = binomial(link= "logit"))
  
  data.c.ipw$pn.csrk <- predict(p.c.num,data.c.ipw,type="response")
  rm(p.c.num)
  
  # 2.3 Computation of estimated weights 
  
  data.c.ipw[event.c==0, sw.csrk:= pn.csrk/pd.csrk]
  data.c.ipw[event.c==1, sw.csrk:= (1-pn.csrk)/(1-pd.csrk)]
  data.c.ipw[, sw.c := prod(sw.csrk), by ='VSIMPLE_INDEX_MASTER']
  sw.c = data.c.ipw[time==1,.(VSIMPLE_INDEX_MASTER ,sw.c)]
  
  rm(data.c.ipw) 
  
  
  data_sub <- merge(data_sub,sw.c,by='VSIMPLE_INDEX_MASTER')
  
  data_sub[,sw:=  sw.a * sw.c * sw.select]
  
  return(data_sub)
}

########################################
# 3. Compute Relative Risk (RR) from hazard model

getRR <- function(data, indices){
  data_sub <- data[indices,]
  
  # 3.1 creation of person-month data for effect model;  
  # Note: this is diff from data.c.ipw
  
  idx <- rep(seq(1,dim(data_sub)[1]),data_sub$survtime)
  data.ipw <- data_sub[idx,]
  data.ipw$time <- sequence(rle(as.character(data.ipw$VSIMPLE_INDEX_MASTER))$lengths)-1
  data.ipw$event <- ifelse((data.ipw$time == data.ipw$survtime-1) & data.ipw$cvd==1,
                           1, 0)
  data.ipw$timesq <- data.ipw$time^2
  head(data.ipw[survtime==1,],2)
  
  
  
  # 3.2 fit the weighted hazards model
  
  ipw.model.linint <- glm(event==0 ~ treat +  time + I(treat*time),
                          family = binomial(),
                          weights = sw, 
                          data = data.ipw)
  
  summary(ipw.model.linint)
  maxSur = max(data.ipw$survtime)
  
  # 3.3 creation of survival curves
  ipw.treat0 <- data.frame(cbind(seq(0,maxSur),0,(seq(0,maxSur))^2))
  ipw.treat1 <- data.frame(cbind(seq(0,maxSur),1,(seq(0,maxSur))^2))
  
  colnames(ipw.treat0) <- c("time","treat","timesq")
  colnames(ipw.treat1) <- c("time","treat","timesq")
  
  ipw.treat0$p.noevent0 <- predict(ipw.model.linint,ipw.treat0,type = "response")
  ipw.treat1$p.noevent1 <- predict(ipw.model.linint,ipw.treat1,type = "response")
  
  
  # 3.4 computation of survival for each person-month 
  ipw.treat0$surv0 <- cumprod(ipw.treat0$p.noevent0)
  ipw.treat1$surv1 <- cumprod(ipw.treat1$p.noevent1)
  
  
  ipw.graph <- merge(ipw.treat0, ipw.treat1, by =c('time','timesq'),
                     suffixes = c('_treat0','_treat1'))
  ipw.graph$survdiff <- ipw.graph$surv1-ipw.graph$surv0
  ipw.graph$RR = (1-ipw.graph$surv1)/(1-ipw.graph$surv0)
  
  RR60 = ipw.graph[ipw.graph$time==60,'RR']
  
  # write.csv(ipw.graph[order(as.numeric(ipw.graph$time)),],RRoutputfname, row.names = FALSE)
  
  return(RR60)
  
}



#######################################

getSubData <- function(data,treatment,treatType){
  
  strict = 'no'  # stricter def of treatment; effect reinforced
  
  
  # table(DATA$censorNoRec)
  # head(DATA,2)
  DATA <- copy(data)
  setDT(DATA)
  useSectionWeight = TRUE
  
  if( useSectionWeight == FALSE){
    DATA[,sw.select := 1]
  }
  
  
  
  if (treatment == 'smoking'){
    # keep only current smokesr or ex-smokers 
    DATA[, keep := cur_smoke_BL ==1]
  }else if (treatment== 'sbp'){
    #  keep <-  rep(TRUE,dim(DATA)[1])
    #  DATA[, keep := (sbp_BL >= 140) & ( obplm_BL ==0)]
    DATA[, keep := (sbp_BL >= 140) ]
  }else if (treatment== 'cholesterol'){
    DATA[, keep := (  imp_index2y_tchdl_ratio_BL>=5)]
  }
  
  
  data_sub = DATA[keep==TRUE,c('VSIMPLE_INDEX_MASTER',
                               'survtime', 
                               'cvd', 'censored',
                               'FEM', 'MAL',
                               'view_ag_sex',
                               'view_ag_age_BL',
                               'view_ag_age' ,
                               'asian','indian','maori','pacific','european',
                               'view_en_nzdep', 
                               'ex_smoke_BL' ,
                               'cur_smoke_BL' ,
                               'pt_diabetes_BL' ,
                               'pt_atrial_fibrillation_BL' ,
                               'obplm_BL','ollm_BL','oatm_BL',
                               'sbp_BL' , 
                               'imp_index2y_tchdl_ratio_BL' ,
                               'pt_familyhistory' ,
                               'imp_hx_diabetes_BL' ,
                               'imp_hx_af_BL' ,
                               'imp_hx_lipidlowering_BL' ,
                               'imp_hx_antithrombotics_BL' ,
                               'imp_hx_antihypertensives_BL',
                               'ex_smoke' ,
                               'cur_smoke' ,
                               'pt_diabetes' ,
                               'pt_atrial_fibrillation' ,
                               'obplm','ollm','oatm',
                               'sbp' , 
                               'imp_index2y_tchdl_ratio' ,
                               'imp_hx_diabetes' ,
                               'imp_hx_af' ,
                               'imp_hx_lipidlowering' ,
                               'imp_hx_antithrombotics' ,
                               'imp_hx_antihypertensives',
                               'sw.select') ]
  
  
  
  
  if (treatment == 'smoking'){
    data_sub$treat <- data_sub$ex_smoke
  }
  
  if (treatment== 'sbp'){
    
    if( (strict == 'yes') &  (treatType == 'medication')){
      data_sub$treat <- 1*(data_sub$obplm &( data_sub$sbp<=130))}
    
    if( (strict == 'no') &  (treatType == 'medication')){
      data_sub$treat  <- data_sub$obplm}
    
    if( (strict == 'no') &  (treatType == 'level')){
      data_sub$treat <- 1* (data_sub$sbp <=130 ) }         
    
    if( (strict == 'yes') &  (treatType == 'level')){
      print(1)
      data_sub$treat <- 1*( data_sub$obplm & (data_sub$sbp <=130 )) }         
    
  }
  
  
  
  if (treatment== 'cholesterol'){
    if( (strict == 'yes') &  (treatType == 'medication'))
    {data_sub$treat <- 1*( data_sub$ollm & (data_sub$imp_index2y_tchdl_ratio <=3.5) )
    }
    if ((strict == 'no') &  (treatType == 'medication')){
      data_sub$treat <- 1*( data_sub$ollm  )
    }
    if ((strict == 'no') &  (treatType == 'level')){
      data_sub$treat <- 1* ((data_sub$imp_index2y_tchdl_ratio_BL - 
                               data_sub$imp_index2y_tchdl_ratio) >=1.5 )
    }
    if ((strict == 'yes') &  (treatType == 'level')){
      data_sub$treat <- 1*  (data_sub$ollm & ((data_sub$imp_index2y_tchdl_ratio_BL - 
                                                 data_sub$imp_index2y_tchdl_ratio) >=1.5 ))
    }
  }
  
  return(data_sub) 
}

