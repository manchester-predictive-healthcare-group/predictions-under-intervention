library(data.table)
library(plyr)
library(survival)



DATA <- read.csv('data_for_refitPredict.csv')
 

dataF <- DATA[FEM==1,c('VSIMPLE_INDEX_MASTER',
                       'survtime',
                       'cvd','view_ag_age',
                       'FEM', 'MAL',
                       'asian','indian','maori','pacific','european',
                       'view_en_nzdep',
                       'ex_smoke' ,
                       'cur_smoke' ,
                       'imp_hx_diabetes' ,
                       'imp_hx_af' ,
                       'pt_familyhistory' ,
                       'sbp',
                       'imp_index2y_tchdl_ratio',
                       'imp_hx_lipidlowering' ,
                       'imp_hx_antithrombotics' ,
                       'imp_hx_antihypertensives')] 

 


dataM <- DATA[MAL==1,c('VSIMPLE_INDEX_MASTER',
                       'survtime',
                       'cvd',
                       'view_ag_age',
                       'FEM', 'MAL',
                       'asian','indian','maori','pacific','european',
                       'view_en_nzdep',
                       'ex_smoke' ,
                       'cur_smoke' ,
                       'imp_hx_diabetes' ,
                       'imp_hx_af' ,
                       'pt_familyhistory'  ,
                       'sbp.c', 'sbp',
                        'imp_index2y_tchdl_ratio',
                       'imp_hx_lipidlowering' ,
                       'imp_hx_antithrombotics' ,
                       'imp_hx_antihypertensives')] 


head(dataM)
# sanity check
sum(dataF$MAL)
sum(dataM$FEM)


colnames((DATA))

setDF(dataM)
setDF(dataF)
###################################################
###################################################

survF <- coxph(Surv(survtime,cvd==1) ~ view_ag_age+
                 maori + pacific +indian + asian + 
                 view_en_nzdep +
                 ex_smoke + cur_smoke + 
                 as.factor(pt_familyhistory) +
                 as.factor(imp_hx_af) +
                 as.factor(imp_hx_diabetes)  +
                 sbp+
                 imp_index2y_tchdl_ratio +
                 as.factor(imp_hx_antihypertensives)+
                 as.factor(imp_hx_lipidlowering) +
                 as.factor(imp_hx_antithrombotics) +
                 view_ag_age*as.factor(imp_hx_diabetes) + 
                 view_ag_age*sbp+
                 sbp*as.factor(imp_hx_antihypertensives),
                 data = dataF)



survM <- coxph(Surv(survtime,cvd==1) ~ view_ag_age+
                 maori + pacific +indian + asian + 
                 view_en_nzdep +
                 ex_smoke + cur_smoke + 
                 as.factor(pt_familyhistory) +
                 as.factor(imp_hx_af) +
                 as.factor(imp_hx_diabetes)  +
                 sbp+
                 imp_index2y_tchdl_ratio +
                 as.factor(imp_hx_antihypertensives)+
                 as.factor(imp_hx_lipidlowering) +
                 as.factor(imp_hx_antithrombotics) +
                 view_ag_age*as.factor(imp_hx_diabetes) + 
                 view_ag_age*sbp+
                 sbp*as.factor(imp_hx_antihypertensives),
               data = dataM)


cF <- as.data.frame(survF$coefficients)

colnames(cF) = 'Female'
cM <- as.data.frame(survM$coefficients)
colnames(cM) = 'Male'


rNames <- c(
  'Age',
  'Maori',
  'Pacific',
  'Indian',
  'Asian',
  'NZDep',  'Ex-smoker',
  'Current-smoker',
  'FamilyCVD',
  'HistAF',
  'HistDiab',
  'SBP',
  'TCHDL',
  "OBPLM", 
  "OLLM", 
  "OATM",
  'Age_HistDiab',
  'Age_SBP',
  'OBPLM_SBP')

rownames(cF) <- rNames
rownames(cM) <- rNames

cFM <- merge(cF,cM,by=0)
# write.table(cFM,'refitPREDICT.csv',sep=',')






############################################################################
  
model <- copy(survM) 
data <- copy(dataM)

#### Baseline survival
y5_cox <- summary(survfit(model),time=5)$surv
y5_cox


predLP <- predict(model,data)
prob <- y5_cox^exp(predLP)

head(1-prob)

head(predLP)

summary(predLP)

# C-statistics
summary(model)$concordance[1]
summary(model)$concordance[1] - 1.96*summary(model)$concordance[2]
summary(model)$concordance[1] + 1.96*summary(model)$concordance[2]
 
# calibration slop
model2 <- coxph(Surv(survtime,cvd==1) ~ predLP,data=data)
model2$coef

# 
predLP.cent <- cut(predLP,breaks=quantile(predLP,probs = c(0,0.16,0.5,0.84,1)))
#predLP.cent

plot(survfit(Surv(survtime,cvd==1)~predLP.cent,data=data),
     col=c(1:4),main="Kaplan-Meier survival estimate",xlab="Analysis time")
legend(5,.2,c('q1','q2','q3','q4'),col=c(1:4),lty=1,bty="n")


#### Baseline survival
y5_cox <- summary(survfit(model),time=5)$surv
y5_cox

# Heuristic shrinkage
null_model <- coxph(Surv(survtime,cvd==1)~1, data=data)
chi2 <- anova(null_model,model)$Chisq[2]
chi2
# which is the same as in:
df_surv <- summary(model)$logtest[2]
df_surv
# Heuristic shrinkage of van Houwelingen: (chi2-df)/chi2
vanH <- (chi2-df_surv)/chi2
vanH
heuristic_LP  <- vanH*predLP
summary(predLP)
summary(heuristic_LP)

# recalculate calibration slope using the shrunken lp
coxph(Surv(survtime,cvd==1) ~ heuristic_LP,data=data)

#### Baseline survival for shrunk model
model_shrunk <- coxph(Surv(survtime,cvd==1) ~ offset(heuristic_LP),
                      data=data)
model_shrunk
y5_cox_shrunk <- summary(survfit(model_shrunk),time=5)$surv
y5_cox_shrunk
y5_cox

prob <- y5_cox^exp(predLP)
prob_shrunk <- y5_cox_shrunk^exp(heuristic_LP)


#plot(prob,prob_shrunk)
summary(prob-prob_shrunk)

 

