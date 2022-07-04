library(survival)
library(stringr)
library(reshape2)
library(data.table)
library(plyr)

source("refitPREDICT.R")
 

############################################################################
# Compute RR for treated vs untreated
 
# sbp/td:hdl reduction (regardless of medications, scenario IIa or IIIa)  
treatType = 'level' 

# sbp/td:hdl medication initiation (regardless of effecgts, scenario IIb or IIIb)
# treatType =  'medication' 


dim(DATA)
dim(dataM)
dim(dataF)

RR_sub = NULL
ARRsum = NULL
AR_ARR_all = NULL
RR_all = NULL

 

for (treatment in c("smoking","sbp","cholesterol")){  #
  dataInt <- copy(DATA)
  setDT(dataInt)
  
  dataInt[, imp_index2y_tchdl_ratio_org:=imp_index2y_tchdl_ratio]
  dataInt[, imp_hx_antihypertensives_org:=imp_hx_antihypertensives]
  dataInt[, sbp_org:=sbp]
  dataInt[, imp_hx_lipidlowering_org:=imp_hx_lipidlowering]
  dataInt[, cur_smoke_org:=cur_smoke]
  dataInt[, ex_smoke_org:=ex_smoke]
  
  
  if (treatment == 'cholesterol'){
    # cholesterol intervention
    if( treatType == 'level'){
    
      IndInt = dataInt$imp_index2y_tchdl_ratio_org >=5
      # tchdl reduce to 3.5
      dataInt[IndInt,imp_index2y_tchdl_ratio_org:=3.5]
    
    } else{
      
      IndInt = (dataInt$imp_index2y_tchdl_ratio_org >=5) & (dataInt$imp_hx_lipidlowering_org==0)
      dataInt[IndInt,imp_index2y_tchdl_ratio_org:=3.5]
      dataInt[IndInt,imp_hx_lipidlowering:=1]
    }
    
    
    
  }
  
  if (treatment == 'sbp'){
    if (treatType == 'level'){
      # sbp int
      IndInt = dataInt$sbp_org>=140 
      # sbp reduce to 130
      dataInt[IndInt,sbp := 130]
    } else{
      IndInt = (dataInt$sbp_org>=140) & (dataInt$imp_hx_antihypertensives_org==0)
      
      dataInt[IndInt,sbp := 130]
      dataInt[IndInt,imp_hx_antihypertensives:=1]
   }
     
  }
  
  if (treatment == 'smoking'){
    # # smoking
    IndInt = dataInt$cur_smoke_org==1
    dataInt[IndInt,cur_smoke:=0]
    dataInt[IndInt, ex_smoke :=1]
  }
  
  setDF(dataInt)
  
  #colnames(dataInt)
  
  #### Baseline survival
  
  RR = NULL
  AR_ARR = NULL
  
  
  for (sex in c("F","M")){
    if (sex == "F"){
      model <-  copy(survF )
      d <- DATA[(DATA$FEM == 1) & IndInt ,]
      d_int <- dataInt[(dataInt$FEM==1) & IndInt,]
      
    }else{
      model <-  copy(survM )
      d <- DATA[(DATA$MAL == 1) & IndInt,]
      d_int <- dataInt[(dataInt$MAL==1) & IndInt,]
    }
    
    y5_cox <- summary(survfit(model),time=5)$surv
    #print(y5_cox)
    
    
    predLP <- predict(model,d,type='lp' ,reference = 'sample')
    
    prob <- y5_cox^exp(predLP)
    #print(summary(prob))
    #d[prob<0.4082,]
    #d[prob>0.9995,]
    
    predLPInt <- predict(model,d_int,type='lp')
    probInt <- y5_cox^exp(predLPInt) # survival
    
    #print(summary(probInt))
    
    SR = probInt/prob 
    length(probInt)
    length(prob)
    
    
    
    rr_modinp = (1-probInt)/(1-prob)
    summary(rr_modinp)
    d_int$RR_modinp = rr_modinp
    d_int$AR_BL = 1-prob
    d_int$AR_modinp = 1- probInt 
    d_int$ARR_modinp = probInt - prob  
    

    if (treatType == 'level') 
      {rr_others = rbind(c( .63,.72,.70,.53),
                     c(.58,.90,.93,.80),
                     c(.62,.92, .95,.70))}
    else{
      rr_others = rbind(c( .63,.72,.70,.53),
                        c(1.11,1.24,1.17,.71),
                        c(.86,1.17, 1.10,.70))
    }
    rr_others  = data.frame(rr_others,row.names = c('smoking','sbp','cholesterol'))
    colnames(rr_others) = c('ipw_bl','Post_bl_itm', 'Post_bl_full','ext')
    
    rr = rr_others[treatment,]
    rr
    for(i in seq_along(rr)){
      rr_i = rr[[i]]
      n = names(rr[i])
      d_int[,paste0('AR_',n)] = rr_i*(1-prob)
      d_int[,paste0('ARR_',n)] = (1-rr_i)*(1-prob)
    }
    
    AR_ARR =  rbind(AR_ARR,d_int[,c(grep('AR',names(d_int),value = TRUE),
                                 'european', 'view_ag_age', 'view_ag_sex')])
    
    RR = rbind(RR, d_int[,c('RR_modinp',
                                    'european', 'view_ag_age', 'view_ag_sex')])  # stack F with M

    
  
  }
  RR$intervention = treatment
  RR_all = rbind(RR_all,RR)
  
  AR_ARR$intervention = treatment
  AR_ARR_all = rbind(AR_ARR_all,AR_ARR)
  
}


dim(RR_all)
head(RR_all)
dim(AR_ARR_all)


##### Table: RR

datalong = copy(RR_all)
datalong$intervention = factor(datalong$intervention,
                               levels = c('smoking','sbp','cholesterol'),
                               labels =  c('smoking'= 'smk cessation',
                                           'sbp' = 'BP',
                                           'cholesterol'='cholesterol' ))

mu = ddply(datalong,.(intervention),
           summarise,grp.mean=mean(RR_modinp,na.rm = TRUE),
           grp.sd=sd(RR_modinp,na.rm = TRUE))


mu


##### Table: AR




AR_ARR_all$Sex = AR_ARR_all$view_ag_sex
AR_ARR_all$Age = ifelse(AR_ARR_all$view_ag_age<=50,'<=50','>50')
AR_ARR_all$Ethnicity = ifelse(AR_ARR_all$european==1,'European','Non-European')



head(AR_ARR_all)

lab = 'AR_'
varbs = grep(lab,names(AR_ARR_all),value = TRUE)

datalong = AR_ARR_all[,c(varbs,'intervention','Sex','Age','Ethnicity')]

head(datalong)

# library(ggplot2)
# ggplot(datalong,aes(y=value, x=Method)) +
#   geom_point() +
#   geom_jitter() +
#   facet_grid(rows = vars(intervention))

varbs = str_replace_all(varbs,lab,'')

colnames(datalong) = c(varbs,'intervention','Sex','Age','Ethnicity')

datalong = melt(datalong,measure.vars= varbs, 
              id.vars = c('intervention',
                          'Sex','Age','Ethnicity'),
     variable.name = 'Method',
     value.name = 'value')


head(datalong)


datalong$Method = factor(datalong$Method,
                         levels = c('BL','modinp','ipw_bl',
                                    'Post_bl_itm', 'Post_bl_full','ext'),
                         labels =  c('BL'= 'Risk w/o INT',
                                      'modinp'='Mod. Input',
                                     'ipw_bl' = 'IPW baseline',
                                     'Post_bl_itm' ='IPW post-BL interm.', 
                                     'Post_bl_full' ='IPW post-BL full',
                                     'ext' = 'External source'))

datalong$intervention = factor(datalong$intervention,
                               levels = c('smoking','sbp','cholesterol'),
                               labels =  c('smoking'= 'smk cessation',
                                           'sbp' = 'BP',
                                           'cholesterol'='cholesterol' ))
datalongOrg = copy(datalong)

mu = ddply(datalong,.(intervention,Method),
           summarise,grp.mean=mean(value,na.rm = TRUE),
           grp.sd=sd(value,na.rm = TRUE))




mu

mu$meansd = paste0(format(mu$grp.mean*100,digits = 2),
                   '% (',format(mu$grp.sd*100,digits = 2),'%)')

mu =  dcast(mu, intervention~Method,value.var = c('meansd') ,drop = TRUE)

mu
 

head(datalong)

#### print different methods

datalong =  copy(datalongOrg)

datalong = datalong[datalong$Method %in% c('Mod. Input'), ]
xlab = 'Modifying input: Absolute risk reduction'
mu_all = data.frame(intervention = unique(datalongOrg$intervention))

head(datalong)

for (g in c('Sex','Age','Ethnicity' )){
  data_g = copy(datalong)
  data_g$group = data_g[,g]
  mu = ddply(data_g,.(intervention,Method,group),
             summarise,grp.mean=mean(value,na.rm = TRUE),
             grp.sd=sd(value,na.rm = TRUE))
  
  
  
  
  
  
  mu$meansd = paste0(format(mu$grp.mean,digits = 3),' (',format(mu$grp.sd,digits = 3),')')
  mu =  dcast(mu, intervention~group,value.var = c('meansd') ,drop = TRUE)
  mu_all = merge(mu_all, mu, all.y=TRUE, by  = 'intervention') 
                 
}


mu_all





# plot AR


ARvarbs = grep('AR_',names(AR_ARR_all),value = TRUE)


datalong = AR_ARR_all[,c(ARvarbs,'intervention','Sex','Age','Ethnicity')]
 
head(datalong)



datalong = melt(datalong,measure.vars= ARvarbs[ARvarbs!='AR_BL'], 
                id.vars = c('AR_BL','intervention',
                            'Sex','Age','Ethnicity'),
                variable.name = 'Method',
                value.name = 'AR')

head(datalong)

datalong = melt(datalong, measure.vars= c('AR_BL','AR'), 
              id.vars = c('Method','intervention'),
              variable.name = 'BLvsINT',
              value.name = 'AR')
head(datalong)


datalong$BLvsINT = factor(datalong$BLvsINT,
                        levels = c('AR_BL','AR'),
                        labels =  c('AR_BL'= 'w/o INT',
                                    'AR'='wt INT'))

datalong$Method = factor(datalong$Method,
                        levels = c('AR_modinp','AR_ipw_bl',
                                   'AR_Post_bl_itm', 'AR_Post_bl_full','AR_ext'),
                        labels =  c('AR_modinp'='Mod. Input',
                                    'AR_ipw_bl' = 'IPW baseline',
                                    'AR_Post_bl_itm' ='IPW post-BL interm.', 
                                    'AR_Post_bl_full' ='IPW post-BL full',
                                    'AR_ext' = 'External source'))

datalong$intervention = factor(datalong$intervention,
                        levels = c('smoking','sbp','cholesterol'),
                        labels =  c('smoking'= 'smk cessation',
                                    'sbp' = 'BP',
                                    'cholesterol'='cholesterol' ))




mu = ddply(datalong,.(intervention,Method,BLvsINT),
           summarise,grp.mean=mean(AR,na.rm = TRUE),
           grp.len = length(AR),
           grp.sd=sd(AR,na.rm = TRUE))

mu

p <- ggplot(datalong,
            aes(x=AR,fill=BLvsINT,color=BLvsINT )) +
  geom_histogram(  alpha=.5,binwidth = 0.20/25,position = 'dodge') +
  xlim(0,.2)  + xlab('Absolute Risk') +
  #ylim(0,2000) + 
  geom_vline(data=mu,
             aes(xintercept = grp.mean,color=BLvsINT), 
             linetype='dashed') + 
  theme(text = element_text(size=14)) + 
  scale_color_discrete(name = 'wt or w/o \n intervention') +
  scale_fill_discrete(name = 'wt or w/o \n intervention') +
  facet_grid(col = vars(intervention),rows =vars( Method))


p


#### plot ARR

varbs = grep('ARR_',names(AR_ARR_all),value = TRUE)


datalong = AR_ARR_all[,c(varbs,'intervention')]


datalong = melt(datalong,measure.vars= varbs, 
              id.vars = c( 'intervention'),
              variable.name = 'Method',
              value.name = 'value')


head(datalong)


datalong$Method = factor(datalong$Method,
                       levels = c('ARR_modinp','ARR_ipw_bl',
                                  'ARR_Post_bl_itm', 'ARR_Post_bl_full','ARR_ext'),
                       labels =  c('ARR_modinp'='Mod. Input',
                                   'ARR_ipw_bl' = 'IPW baseline',
                                   'ARR_Post_bl_itm' ='IPW post-BL interm.', 
                                   'ARR_Post_bl_full' ='IPW post-BL full',
                                   'ARR_ext' = 'External source'))

datalong$intervention = factor(datalong$intervention,
                             levels = c('smoking','sbp','cholesterol'),
                             labels =  c('smoking'= 'smk cessation',
                                         'sbp' = 'BP',
                                         'cholesterol'='cholesterol' ))


mu = ddply(datalong,.(intervention,Method),
           summarise,grp.mean=mean(value,na.rm = TRUE),
           grp.len = length(value),
           grp.sd=sd(value,na.rm = TRUE))

mu

xlab = 'Absolute Risk Reduction'
library(ggplot2)
p <- ggplot(datalong,
            aes(x=value)) +
  geom_histogram( fill='red', alpha=.3,
                  binwidth = 0.03/60,position = 'identity') +
  xlim(0,0.03)  + xlab(xlab) +
  #ylim(0,2000) + 
  geom_vline(data=mu,
             aes(xintercept = grp.mean,color='red'), 
             linetype='dashed') + 
  theme(legend.position = 'none',text = element_text(size=14)) + 
  facet_grid(cols = vars(intervention),rows =vars( Method))


 
##### 



# 
# if (treatType == 'level'){
# write.table(RR_all,paste0(outDir,'refitPREDICT_RR_treatbyLevels.csv'),sep=',')
# } else {
#   write.table(RR_all,paste0(outDir,'refitPREDICT_RR_treatbyMedInit.csv'),sep=',')
#   
# }
 


# standard error of effect estimated

model = copy(survM)
data = copy(dataM)

summary(model)
model$means

L = dim(data)[2]

d1 = data[1,]
d2 = copy(d1)
d2$cur_smoke =0
d2$ex_smoke =1



d3 = copy(d2)
d3[,2:L] = d2[,2:L] - d1[,2:L]

d4 = copy(d3)


for (c in colnames(d3)[2:L]){
  d4[,c] = d4[,c] + model$means[c]
  
}

###
d1
d5 = copy(d1)
d2$cur_smoke =0
d2$ex_smoke =1



d3 = copy(d2)
d3[,2:L] = d2[,2:L] - d1[,2:L]

d4 = copy(d3)


for (c in colnames(d3)[2:L]){
  d4[,c] = d4[,c] + model$means[c]
  
}




d = rbind(d1,d2,d3,d4)
d

Predd = predict(model,d,type='risk' ,reference = 'sample',se.fit = T)
Predd

MM = data.frame(model$means)
tMM = transpose(MM)
colnames(tMM) = rownames(MM)
rownames(tMM) = colnames(MM)

PredMM = predict(model,tMM,type='risk' ,reference = 'sample',se.fit = T)
Predd$fit
PredMM$fit*Predd$fit

Predd$fit[2]/Predd$fit[1]


preddLP <- predict(model,d,type='lp' ,reference = 'sample',se.fit = T)
preddLP 

predexpLP<- predict(model,data,type='risk' ,reference = 'sample',se.fit = T)

predexpected<- predict(model,data,type='expected' ,reference = 'sample',se.fit = T)

predexpected$fit[1:10]
predexpected$se.fit[1:10]


predLP$fit[1:10]
exp(predLP$fit[1:10])
predexpLP$fit[1:10]

predLP$se.fit[1:10]

predexpLP$se.fit[1:10]










