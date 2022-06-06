print(getwd())
outDir <- paste0(getwd(), '/Output/')
print(outDir)

library(data.table)
library(survival)
library(boot)


###################################################################
# Compute RR & CI from the IPW model
###################################################################



# DATA = read.csv('data_for_postBLTreatment_2years_addSelectionWeight.csv')
# fname <- paste0(outDir, 'IP_2years_boots_RR_LL_UL.csv')

# DATA = read.csv('data_for_postBLTreatment_2ndVisit_addSelectionWeight.csv')
# fname <- paste0(outDir, 'IP_2ndVisit_boots_RR_LL_UL.csv')


RR_CI = data.frame(row.names = c('RR','SE','LL','UL'))

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
             
             
            
            boots <- data.frame(i = 1:N, RR = NA, SE =NA, LL = NA, UL = NA)
            
            fname <- paste0(outDir, 'IPW_RR_boots_',treatment,'_',treatType,
                            ifelse(ifFullModel,'_full','_interm'), '.csv')
            cat(fname,'\n')
            
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
 
            
            write.csv(boots, fname,row.names = FALSE)
            
            res  = as.data.frame(apply(boots[,c(-1)],2, function(x) exp(mean(log(x)))))
            
            colnames(res) = paste0(treatment,'_',treatType,
                           ifelse(ifFullModel,'_full','_interm')) 
            
            
            RR_CI = merge(RR_CI,res,by = 'row.names')  
            rownames(RR_CI) = RR_CI$Row.names
            RR_CI$Row.names = NULL
            

        }
    }
}



write.csv(RR_CI[c('SE','RR','LL','UL'),],fname,row.names = TRUE)

