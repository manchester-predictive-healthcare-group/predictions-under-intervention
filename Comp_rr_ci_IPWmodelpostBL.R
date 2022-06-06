print(getwd())
outDir <- paste0(getwd(), '/Output/')
print(outDir)

library(data.table)
library(survival)
library(boot)

###################################################################
# Compute RR & CI from the IPW model: using post-baseline data
###################################################################

source("EffectEstimate_IPW_postBL.R")



DATA = read.csv('data_for_postBLTreatment_2years_addSelectionWeight.csv')
fname <- paste0(outDir, 'IP_postBL_2years_boots_RR_LL_UL.csv')

# DATA = read.csv('data_for_postBLTreatment_2ndVisit_addSelectionWeight.csv')
# fname <- paste0(outDir, 'IP_postBL_2ndVisit_boots_RR_LL_UL.csv')

boots <- data.frame(i = 1:10, row.name = NA, RR = NA, SE =NA, LL = NA, UL = NA)


for (treatment in c('smoking','sbp','cholesterol')){
  
    for (treatType in c('level', 'medication')){
        
        if (treatment == 'smoking' & treatType == 'medication'){
          next # we don't have anything special for smoking medication
        }
        for (ifFullModel in c(FALSE,TRUE)){


            data_sub = getSubData(DATA,treatment,treatType)
            
            
            table(data_sub$cvd)
            dim(data_sub)
            
            
            ## bootstraping for getting CI for RR
            N = nrow(data_sub)

            d <- copy(data_sub)
            d <- getTreatW(d,treatment,ifFullModel)
            d <- getCensW(d,treatment)
            
            results <- boot(data=d,statistic = getRR,R=500)
            m = results$t0
            se <- sd(results$t[,1])
            ll <- m - qnorm(0.975)*se
            ul <- m + qnorm(0.975)*se
            boots[i,'RR'] = m
            boots[i,'SE'] = se
            boots[i,'LL'] = ll
            boots[i,'UL'] = ul
 
           
            boots[i,'row.names'] = paste0(treatment,'_',treatType,
                           ifelse(ifFullModel,'_full','_interm')) 

        }
    }
}



write.csv(boots,fname,row.names = TRUE)

