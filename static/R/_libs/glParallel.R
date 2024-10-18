# a function for optimism-correct AUC
# by resampling the same size of the original sample with B times
# http://cainarchaeology.weebly.com/r-function-for-optimism-adjusted-auc.html
auc.adjust <- function(data, fit, B){
    # get overfitted AUC using the orignal data and the fitted model (probably glm)
    fit.model <- fit
    data$pred.prob <- fitted(fit.model)
    #data$pred.prob <- predict(fit.model,type=c("response")) # same as the above
    auc.app <- pROC::roc(data[,1], data$pred.prob, data=data)$auc # require 'pROC'. 
                                                                  # NB, "y" (outcome) is the *first* column) 

    # set the following three vectors with the length of B
    auc.boot <- vector (mode = "numeric", length = B)
    auc.orig <- vector (mode = "numeric", length = B)
    o <- vector (mode = "numeric", length = B)

    # bootstrap for B times
    for(i in 1:B){
        boot.sample <- kimisc::sample.rows(data, nrow(data), replace=TRUE) # require 'kimisc'
        fit.boot <- glm(formula(fit.model), data = boot.sample, family = "binomial")
        boot.sample$pred.prob <- fitted(fit.boot)
        # get bootstrapped AUC
        auc.boot[i] <- pROC::roc(boot.sample[,1], boot.sample$pred.prob, data=boot.sample)$auc
        # get original data boot AUC
        data$pred.prob.back <- predict.glm(fit.boot, newdata=data, type="response")
        auc.orig[i] <- pROC::roc(data[,1], data$pred.prob.back, data=data)$auc
        # calculated optimism corrected version
        o[i] <- auc.boot[i] - auc.orig[i]
    }
    auc.adj <- auc.app - (sum(o)/B)
    return(auc.adj)
}

# optimism-adjusted c-stat 
# using simple bootstrap resampling such as k-fold repeated CV
get_KFCV<-function(x,mc.cores=20){
        my.data=x # isa data.frame # data.frame(x[,c(my.subsets,"y")])
        my.data$y<-factor(ifelse(my.data$y==1,'case','non_case'),levels=c("non_case","case"))

        cl <- makePSOCKcluster(mc.cores) # No. of cores to use
        registerDoParallel(cl)
        set.seed(333)
        model<-caret::train(
                as.formula( paste( 'y', '~', '.' ) ), 
                data=my.data,
                method="glm",
                family="binomial",
                trControl = trainControl(method = "repeatedcv",
                                            number=10,repeats=100,
                                            summaryFunction = twoClassSummary,
                                            classProbs = TRUE,
                                            savePredictions = T
                                        ),
                metric = "ROC"
        )
        stopCluster(cl)
        model$results$ROC * 100
}

# optimism-adjusted c-stat (LPOCV)
# leave pair out cv
get_LPOCV<-function(x=my.mat){
        my.data=x  # isa data.frame
        my.data$y<-factor(ifelse(my.data$y==1,'case','non_case'),levels=c("non_case","case")) # the outcome

        # case-control grid
        myGrid<-expand.grid(
            rownames(my.data[my.data$y=="case",]),
            rownames(my.data[my.data$y=="non_case",])
        )

        out<-apply(myGrid, 1, function(i){
            #as.numeric(as.vector(t(i))) # a pair of case-control
            #my.data[as.numeric(as.vector(t(i))),] # a pair of case-control

            #my.data1<-my.data[as.numeric(as.vector(t(i))),] # a pair of case-control
            #my.data2<-my.data[-as.numeric(as.vector(t(i))),] # remaining (excluding the pair)

            ##
            my.data1<-my.data[as.vector(t(i)),] # a pair of case-control
            my.data2<-my.data[!rownames(my.data) %in% as.vector(t(i)),] # the remaining (excluding the pair)

            fit<-glm(y~. , data = my.data2, family = "binomial") # fit the model based on the remaining
            predict.glm(fit, newdata=my.data1, type="response") # predict the outcome of the pair using the model above
            #roc(response=my.data2$y, predictor=fitted(fit))$auc
        })

        # the proportion of all pairwise combinations in which the predicted probability was greater for the case than for the control 
        c.stat<-apply(out,2, function(i) i[1] > i[2]) %>% table 
        c.stat["TRUE"] / sum(c.stat) *100 
}

# the above using the `boot` package
# see https://www.statmethods.net/advstats/bootstrapping.html
get_LPOCV_boot<-function(x,index){
        my.data=x[index,]  # isa data.frame
        my.data$y<-factor(ifelse(my.data$y==1,'case','non_case'),levels=c("non_case","case")) # the outcome

        # case-control grid
        myGrid<-expand.grid(
            rownames(my.data[my.data$y=="case",]),
            rownames(my.data[my.data$y=="non_case",])
        )

        out<-apply(myGrid, 1, function(i){
            #as.numeric(as.vector(t(i))) # a pair of case-control
            #my.data[as.numeric(as.vector(t(i))),] # a pair of case-control

            #my.data1<-my.data[as.numeric(as.vector(t(i))),] # a pair of case-control
            #my.data2<-my.data[-as.numeric(as.vector(t(i))),] # remaining (excluding the pair)

            ##
            my.data1<-my.data[as.vector(t(i)),] # a pair of case-control
            my.data2<-my.data[!rownames(my.data) %in% as.vector(t(i)),] # the remaining (excluding the pair)

            fit<-glm(y~. , data = my.data2, family = "binomial") # fit the model based on the remaining
            predict.glm(fit, newdata=my.data1, type="response") # predict the outcome of the pair using the model above
            #roc(response=my.data2$y, predictor=fitted(fit))$auc
        })

        # the proportion of all pairwise combinations in which the predicted probability was greater for the case than for the control 
        c.stat<-apply(out,2, function(i) i[1] > i[2]) %>% table 
        return(c.stat["TRUE"] / sum(c.stat) *100)
}

