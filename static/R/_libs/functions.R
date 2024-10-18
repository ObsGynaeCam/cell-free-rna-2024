
# x: a data.table (or data.frame) which contain 'y' as a outcome along with features in the columns
# value: a data.table having the following columns: feature, AUC and its 95% CI
find_DEP<-function(x,mc.cores=6,pvalue=0.05){
  stopifnot(is.data.frame(x))
  stopifnot(any(grepl("y",colnames(x)))) # should have the column 'y'
  stopifnot(all(sort(unique(data.frame(x)[,"y"]))==c(0,1))) # should be '0' and '1'

  all.proteins<-names(x)[which(names(x)!="y" & names(x)!="ID")] # remove "y" and "ID" (if any)

  # 1. t-test
  dt.ttest<-parallel::mclapply(all.proteins, function(my.protein){
    #message(my.protein)
    my.formula<-formula(paste(my.protein,"~ y",collapse= " "))
    my.test<-t.test(my.formula, data=x)
    dt.ttest <- 
    data.table(
    `feature`=my.protein,
    `case`=my.test$estimate["mean in group 1"], # case
    `ctrl`=my.test$estimate["mean in group 0"], # control
    `se`=my.test$stderr,
    `pval`=my.test$p.value
    ) # return this dt
  },mc.cores=mc.cores) %>% rbindlist

  # 2. logistic regression
  dt.lreg<-parallel::mclapply(all.proteins, function(my.protein){
    #message(my.protein)
    my.formula2<-formula(paste("y ~", my.protein))
    my.model<-glm(my.formula2, data=x, family="binomial")

    foo1<-as.data.frame(coef(summary(my.model))) # Estimate(log_odds=coef), Std. Err, z value, Pr(>[z])
    foo2<-as.data.frame(exp(cbind(coef(my.model), confint(my.model)))) # Odds Ratio & CI (95%)
    foo3<-cbind(my.protein,cbind(foo1,foo2)[my.protein,]) %>% data.table 
    colnames(foo3)<-c("feature","log_odds","se","zval","pval","odds","odds_lo","odds_hi")

    #auc.ci<-pROC::ci.auc(response=x$y,predictor=fitted(my.model), quiet=T) # same as the above
    auc.ci<-pROC::ci.auc(response=x$y[!is.na(x[,..my.protein]) %>% as.vector],
                         predictor=fitted(my.model), quiet=T) # same as the above

    cbind(foo3,auc=auc.ci[2],auc_lo=auc.ci[1],auc_hi=auc.ci[3]) # return this dt
  },mc.cores=mc.cores) %>% rbindlist

  dt.ttest[,BH:=p.adjust(pval,"BH")]
  dt.lreg[,BH:=p.adjust(pval,"BH")]

  DEPs<-c(dt.ttest[BH<=pvalue]$feature, dt.lreg[BH<=pvalue]$feature) %>% unique

  #dt.lreg[feature %in% DEPs,.(feature,auc,auc_lo,auc_hi)][order(-auc)] # return this dt
  merge(dt.ttest, dt.lreg, by="feature")[feature %in% DEPs][order(-auc)] # return this dt
}# end of find_DEP

#' Leave Pair Out Cross Validation
#'
#' This function returns an optimism-adjusted c-stat. Read more by [Gordon Am J Epi 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4108045/).
#' @param x data.frame
#' @return a numeric value of LPOCV 
#' @examples
#' input=read.csv(system.file("extdata","demo_input.csv",package="TSM")) # read the example input 
#' get_LPOCV(x=input[,c("F1","y")]) # get LPOCV of "F1" as a sole predictor
#'
#' get_LPOCV(x=input[,c("F1","F2","y")]) # get LPOCV of "F1" and "F2" as two predictors
#' @export
get_LPOCV<-function(x,mc.cores=1){
  stopifnot(is.data.frame(x))
  stopifnot(any(grepl("y",colnames(x)))) # should have the column 'y'
  stopifnot(all(sort(unique(x[,"y"]))==c(0,1))) # should be '0' and '1'
  my.data=x  # isa data.frame
  my.data$y<-factor(ifelse(my.data$y==1,'case','non_case'),levels=c("non_case","case")) # the outcome

  # case-control grid
  myGrid<-expand.grid(
      rownames(my.data[my.data$y=="case",]),
      rownames(my.data[my.data$y=="non_case",])
  )

  if(T){
    out<-parallel::mclapply(1:nrow(myGrid), function(i){ # for each row
        my.data1<-my.data[as.vector(t(myGrid[i,])),] # a pair of case-control
        my.data2<-my.data[!rownames(my.data) %in% as.vector(t(myGrid[i,])),] # the remaining (excluding the pair)

        fit<-glm(y~. , data = my.data2, family = "binomial") # fit the model based on the remaining
        predict.glm(fit, newdata=my.data1, type="response") # predict the outcome of the pair using the model above
        #roc(response=my.data2$y, predictor=fitted(fit))$auc
    },mc.cores=mc.cores)
    out<-do.call(rbind,out) %>% t
  }else{
    system.time(
    out<-apply(myGrid, 1, function(i){ # for each row
        my.data1<-my.data[as.vector(t(i)),] # a pair of case-control
        my.data2<-my.data[!rownames(my.data) %in% as.vector(t(i)),] # the remaining (excluding the pair)

        fit<-glm(y~. , data = my.data2, family = "binomial") # fit the model based on the remaining
        predict.glm(fit, newdata=my.data1, type="response") # predict the outcome of the pair using the model above
        #roc(response=my.data2$y, predictor=fitted(fit))$auc
    })
    )
  }

  # the proportion of all pairwise combinations in which the predicted probability was greater for the case than for the control 
  c.stat<-table(apply(out,2, function(i) i[1] > i[2]))
  return(c.stat["TRUE"] / sum(c.stat))
} # end of LPOCV

#' Leave Pair Out Cross Validation with CI via the `boot` package
#'
#' This function returns an optimism-adjusted c-stat. Read more by [Gordon Am J Epi 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4108045/).
#' @param x data.frame
#' @return a numeric value of LPOCV 
#' @examples
#' input=read.csv(system.file("extdata","demo_input.csv",package="TSM")) # read the example input 
#' get_LPOCV(x=input[,c("F1","y")]) # get LPOCV of "F1" as a sole predictor
#'
#' get_LPOCV(x=input[,c("F1","F2","y")]) # get LPOCV of "F1" and "F2" as two predictors
#' @export
get_LPOCV_boot<-function(x,index){
        my.data=x[index,]  # isa data.frame
        my.data$y<-factor(ifelse(my.data$y==1,'case','non_case'),levels=c("non_case","case")) # the outcome

        # case-control grid
        myGrid<-expand.grid(
            rownames(my.data[my.data$y=="case",]),
            rownames(my.data[my.data$y=="non_case",])
        )

        out<-apply(myGrid, 1, function(i){
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

#set.seed(333)
#LPOCV.boot<-boot(data=test_data, statistic=get_LPOCV_boot,R=100,parallel="multicore",ncpus=20) # R: the No. of resampling
#LPOCV.ci<-boot.ci(LPOCV.boot, type = "bca")# use "bca" method for bias-corrected and accelerated interval
#c(LPOCV.ci$t0[[1]], LPOCV.ci$bca[4:5]) #LPOCV


#' The Smith Method
#' 
#' This function TSM (aka. The Smith Method) selects a desired number of features (4 by default) by purposefully dropping highly correlated features, *i.e,* picking up a set of representative features that can best explain the binary outcomes. In a plain language, it works like the follwoing: The first representative feature is the one that shows the highest AUC (Area Under the ROC Curve) out of all the features. The next representative feature is the one that shows the highest AUC out of the remaing features after dropping highly correlated features with the first representative feature. The third, the fourth, and so on, represenative feature will be picked up as the same way the 2nd is picked up.
#' @param x Path to the input file 
#' @param method A Character either `pearson` (default) or `spearman`, which is the same paramter `method` for `cor()`.
#' @param corr A numeric vector for the thresholds of correlation coefficients.  
#' @param k The number of desicred features (default:4)
#' @param verbose Boolean
#' @return a data.table (default) or a list of data.table (verbose=T)
#' @examples
#' input=read.csv(system.file("extdata","demo_input.csv",package="TSM")) # read the example input 
#' TSM(x=input) # run TSM with default parameters
#'
#' TSM(x=input, corr=c(0.4, 0.5)) # two correlation coefficients only 
#'
#' TSM(x=input, method="pearson") # pearson method for cor()
#' @import data.table 
#' @import magrittr
#' @export 
TSM<-function(x,method="spearman",corr=seq(0.1,0.7,by=0.1),k=4,verbose=F){
  stopifnot(is.data.frame(x))
  stopifnot(any(grepl("y",colnames(x)))) # should have the column 'y'
  stopifnot(all(sort(unique(x[,"y"]))==c(0,1))) # should be '0' and '1'

  IDs=colnames(x)[!grepl("y",colnames(x))]
  cases<-x[,"y"]==1
  li.top.rank<-list()

  message("calculating AUC for each features...")
  foo<-list()
  for(my.ID in IDs){
      my.mat<-x[,c("y",my.ID)] 
      my.model<-glm(y ~., data=my.mat, family="binomial")

      # ROC & AUC 
      #predict(my.model,type=c("response"))  # probability
      #fitted(my.model) # same as above
      my.prob<-fitted(my.model)
      # in case of NA in the input
      if(nrow(my.mat)!=length(my.prob)){
          dt.prob<-rbind(
          data.table(`index`=as.numeric(names(my.prob)),`prob`=my.prob),
          data.table(`index`=as.numeric(my.mat[,my.ID] %>% is.na %>% which), `prob`=NA)
          )[order(index)]
          my.prob<-dt.prob$prob
      }
      my.mat$prob<-my.prob
      my.roc <- pROC::roc(y ~ prob, data = my.mat, quiet=T)
      foo[[my.ID]] <- data.table(ID=my.ID, auc=my.roc$auc)
  } # end of for   
  dt.auc<-rbindlist(foo)[order(-auc)]

  # for each level of correlation
  for(my.cor in corr){
    cor.index<-paste0("cor",my.cor)
    message(cor.index)

    # round 1
    top.rank<-dt.auc[1]$ID
    mat.cor<-cor(x[cases,IDs],method=method)[top.rank,]
    hi.cor<-abs(mat.cor) > my.cor
    #mat.cor[hi.cor]
    hi.cor.feature<-names(mat.cor[hi.cor])
    li.top.rank[[cor.index]][["top.rank"]]<-top.rank
    li.top.rank[[cor.index]][["cor"]]<-hi.cor.feature
    li.top.rank[[cor.index]][["num.cor"]]<-length(hi.cor.feature)

    # round >=2
    while(length(li.top.rank[[cor.index]][["cor"]]) < length(IDs) -1){
        #print(paste0("round:",round))
        features<-IDs[!IDs %in% li.top.rank[[cor.index]][["cor"]]] # drop highly correlated features (i.e. non-highly correlated features)
        top.rank<-dt.auc[ID %in% features]$ID[1]
        mat.cor<-cor(x[cases,features],method=method)[top.rank,]
        hi.cor<-abs(mat.cor) > my.cor
        #mat.cor[hi.cor]
        hi.cor.feature<-names(mat.cor[hi.cor])
        li.top.rank[[cor.index]][["top.rank"]]<-c(li.top.rank[[cor.index]][["top.rank"]],top.rank)
        li.top.rank[[cor.index]][["cor"]]<-c(li.top.rank[[cor.index]][["cor"]],hi.cor.feature)
        li.top.rank[[cor.index]][["num.cor"]]<-c(li.top.rank[[cor.index]][["num.cor"]],length(hi.cor.feature))
    } # end of while

    ##########################################
    # Performance of k features in the model #
    ##########################################
    num.features<-ifelse(length(li.top.rank[[cor.index]][["top.rank"]])>=k,k,length(li.top.rank[[cor.index]][["top.rank"]]))
    my.features<-li.top.rank[[cor.index]][["top.rank"]][1:num.features]
    my.data=x[,c("y",my.features)] # olink NPX of the features
    my.fit<-glm(y~. , data = my.data, family = "binomial") # fit the model based on the selected features
    li.top.rank[[cor.index]][["fit"]]<-my.fit

    # get the model performance
    li.top.rank[["performance"]][[cor.index]]<-data.table(
                                                        Cor=my.cor,
                                                        `Num features`=length(li.top.rank[[cor.index]][["top.rank"]]),
                                                        Features=paste(li.top.rank[[cor.index]][["top.rank"]],collapse=","),
                                                        `Best features`=paste(my.features,collapse=","),
                                                        AIC=my.fit$aic,
                                                        BIC=BIC(my.fit),
                                                        AUC=pROC::roc(response=x$y, predictor=fitted(my.fit), quiet=T)$auc,
                                                        `AUC(LPOCV)`=get_LPOCV(my.data)
                                                        )
  } # end of for(my.cor)

  if(verbose){
    return(li.top.rank)
  }else{
    return(rbindlist(li.top.rank[["performance"]])[order(-`AUC(LPOCV)`)])
  }
} # end of TSM


##
get_AUC<-function(dt.X,dt.Y){
  df.train <- dt.X %>% as.matrix(rownames="ID") %>% data.frame
  df.test <- dt.Y %>% as.matrix(rownames="ID") %>% data.frame

  my.fit<-glm(y~. , data = df.train, family = "binomial") 
  df.test$prob<- predict(my.fit, newdata = df.test, type = "response")

  auc.ci.tr<-pROC::ci.auc(response=df.train$y,predictor=fitted(my.fit), quiet=T) # same as the above
  auc.ci<-pROC::ci.auc(y ~ prob, data = df.test, quiet=T) #

  my.lpocv<-get_LPOCV(df.train)

  data.table(
             AUC_tr=round(auc.ci.tr[2],4),AUC_lo_tr=round(auc.ci.tr[1],4),AUC_hi_tr=round(auc.ci.tr[3],4),
             LPOCV=round(my.lpocv,4),
             AUC=round(auc.ci[2],4),AUC_lo=round(auc.ci[1],4),AUC_hi=round(auc.ci[3],4)
  )
}

##
get_AUC2<-function(dt.X,dt.Y){
  df.mat1<- dt.X %>% as.matrix(rownames="ID") %>% data.frame
  df.mat2<- dt.Y %>% as.matrix(rownames="ID") %>% data.frame

  my.fit1<-glm(y~. , data = df.mat1, family = "binomial") 
  auc.ci1<-pROC::ci.auc(response=df.mat1$y,predictor=fitted(my.fit1), quiet=T) # same as the above

  my.fit2<-glm(y~. , data = df.mat2, family = "binomial") 
  auc.ci2<-pROC::ci.auc(response=df.mat2$y,predictor=fitted(my.fit2), quiet=T) # same as the above

  #my.lpocv<-get_LPOCV(df.train)
  data.table(
             AUC_tr=round(auc.ci1[2],4),AUC_lo_tr=round(auc.ci1[1],4),AUC_hi_tr=round(auc.ci1[3],4),
             LPOCV=NA,
             AUC=round(auc.ci2[2],4),AUC_lo=round(auc.ci2[1],4),AUC_hi=round(auc.ci2[3],4)
  )
}


#################
# Lasso-pathway #
#################
run_lassoP<-function(x,N=ncol(x)-1){
  stopifnot(is.matrix(x))
  stopifnot(any(grepl("y",colnames(x)))) # should have the column 'y'
  stopifnot(all(sort(unique(data.frame(x)[,"y"]))==c(0,1))) # should be '0' and '1'

  # run Lasso #
  set.seed(333) # set a random seed for a reproducibility
  cv.fit<-glmnet::cv.glmnet(
                  x=x[,!colnames(x) %in% "y"],
                  y= x[,"y"],
                  family="binomial",
                  alpha=1, # default (i.e. lasso)
                  keep=T, # FALSE by default
                  type.measure = "auc" #type.measure="class" # default for 'binomial'
  )
  # Method 1: get the first lambda sequence having all the rows
  mat.beta <- cv.fit$glmnet.fit$beta %>% as.matrix
  my.lambdas<-apply(mat.beta, 2, function(i){table(i!=0)["TRUE"]})>=N
  # the first index having desured number of features 
  this.index<-my.lambdas[my.lambdas & !is.na(my.lambdas)][1] %>% names 
  stopifnot(is.na(this.index)==F)
  #this.index.num <- (strsplit(this.index,"s")[[1]][2] %>% as.integer) +1
  #nZero<-sum(mat.beta[,this.index]!=0, na.rm=T) # number of non-zero coefficient
  dt.bar<-mat.beta[mat.beta[,this.index]!=0,this.index,drop=F] %>% as.data.table(keep.rownames=T)
  dt.bar[,lambda:=colnames(dt.bar)[2]]
  setnames(dt.bar,c("feature","beta","lambda"))
  dt.lasso <- dt.bar[order(-abs(beta))][,rank:=1:.N]

  # Method 2: incrementally find featues from the lambda sequence by favouring firstly found feature
  dt.beta<-mat.beta %>% as.data.table(keep.rownames=T) %>% melt.data.table(id.vars="rn",variable.name="lambda",value.name="beta")
  #dt.beta[,lambda_N:=tstrsplit(lambda,"s",keep=2)[[1]] %>% as.integer]
  dt.foo<-dt.beta[beta!=0,.SD[1],by=rn]
  dt.lasso2<-dt.foo[order(lambda,-abs(beta))][,rank:=1:.N]
  setnames(dt.lasso2,c("feature","lambda","beta","rank"))
  #dt.lasso2[1:N]

  ## Below from Paul
  currentNames       <- NULL
  numberOfSelections <- 1
  maxNumberOfSelections<-N
  allSelections<- vector(mode = "list", length = maxNumberOfSelections)
  for(j in 1:ncol(mat.beta))
  {
    currentBetas <- sort(abs(mat.beta[mat.beta[,j]!=0,j]), decreasing = T)
    newNames <- (names( currentBetas))
    if(length(newNames) > length(currentNames)){
      differenceInNumberOfSelections <- length(newNames) - length(currentNames)
      for(j in 1:differenceInNumberOfSelections)
      {
        allSelections[[numberOfSelections]] <- sort(newNames[1:(j+length(currentNames))])
        numberOfSelections <- numberOfSelections + 1
        if(numberOfSelections > maxNumberOfSelections){break}
      }
      currentNames <- newNames
    }

    if(numberOfSelections > maxNumberOfSelections){break}
  }

  list(`M1`=dt.lasso[1:N],`M2`=dt.lasso2[1:N],`M3`=allSelections) %>% return
}

clean_olinkID<-function(){
  stop("Done already")
  dt.olinkID<-fread("data/dt.olinkID.csv") # this is the old  file
  multi.genes<-c("CXCL8.Oncology","IDO1.Cardiometabolic_II",
                  "IL6.Oncology","LMOD1.Cardiometabolic_II",
                  "SCRIB.Cardiometabolic_II","TNF.Oncology") # n=6
  dt.olinkID[,Assay2:=tstrsplit(Assay,"\\.",keep=1)[[1]]]
  dt.olinkID[,Drop:=F]
  dt.olinkID[Assay2 %in% tstrsplit(multi.genes, "\\.")[[1]],Drop:=ifelse(Assay %in% multi.genes,F,T)]
  dt.olinkID[Assay2=="FLT1",Drop:=T] # OID21301 Oncology   FLT1  P17948  3216   FLT1    Yes
  df.olinkID<- as.matrix(dt.olinkID, rownames="OlinkID") %>% data.frame
  fwrite(dt.olinkID, file="data/dt.olinkID.csv")
}
