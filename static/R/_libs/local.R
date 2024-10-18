# set global options
#options(knitr.table.format="latex") # should pick this up automatically
library(knitr)
library(kableExtra)
library(ggplot2) # graph
library(ggrepel)
library(ggthemes)
library(ggsci)
#library(ggbeeswarm)
#library(ggvenn) # for Venn diagram
library(data.table) # iloveit
library(caret) # train/test workflow 
library(glmnet) # lasso & elastic net 
#library(glmulti) # subset selections via (g)lm
library(magrittr) # for piping (%>%)
library(readxl) # to read excel file
library(doParallel) # for parallel jobs
library(pROC) # ROC calculation
library(ComplexHeatmap) # for heatmap
library(pheatmap) # yet another heatmap
#library(DT)
#library(gprofiler2)
#library(clusterProfiler)
#library(GO.db) # for GOBPPARENTS
#library(GOstats) # for GOGraph() function
#library(R.utils) # to read .gz file via fread
library(GenomicRanges)
library(tximeta)

##
## Plasma samples 
##
# Below from Suzanne Rohrback <srohrback@illumina.com> of illumina
# Re-processed samples. 3 samples were processed through enrichment twice, I kept the data separate from each sequencing run.
#   - GS-B-374-UW is data from the same plasma sample as GS-B-374-b. Both failed QC
#   - GS-59-CX is data from the same plasma sample as GS-59-CX-b. Both failed QC
#   - GS-179-HQ is data from the same plasma sample as GS-179-HQ-b. Only GS-179-HQ-b passed QC
# Re-sequenced samples. 66 samples were sequenced in 2 sequencing runs to obtain sufficient depth.
#   - Their data was processed through DRAGEN RNA twice, once with the first sequencing run, once with data from both sequencing runs combined
#   - They are indicated in the excel file by the “ProcessingBatch” column, if 2 processing batches are listed (eg A0004;B0011), it was sequenced twice
#   - You should only look at bams from the project that corresponds to the second ProcessingBatch, that will contain the full depth data
if(TRUE){
  #dt.meta.info<-readxl::read_excel("data/AS SENT POPS078B with ILMN sample ID.xlsx") %>% as.data.table # from pilot data Steve
  # meta-data from Steve
  dt.foo<-readxl::read_excel("data/preterm_case_control_and_term_fetal_sex_info.xlsx",sheet=1) %>% data.table # preterm (phase1)
  #dt.bar<-readxl::read_excel("data/preterm_case_control_and_term_fetal_sex_info.xlsx",sheet=2) %>% data.table # term (phase2)
  dt.bar<-readxl::read_excel("data/term_grouping_status_n907.xlsx") %>% data.table # term with case info (phase2)
  dt.meta.info<-rbind(
    dt.foo[,.(Type="preterm",setnumber,Batch=Box,Sample_number,IlluminaID,Location=Box,Sublocation,GA=Sampletakenat,PI=PatientInitials,pn_female,case)],
    dt.bar[,.(Type="term",setnumber,Batch=Location,Sample_number,IlluminaID,Location,Sublocation,GA=Sampletakenat,PI=PatientInitials,pn_female,case)]
  )

  dt.fastq<-fread("data/FASTQ.list.r1.txt",col.names=c("SampleID","FastQ")) # All fastq files (n=3421)
  dt.fq.samples<-dt.fastq[!grepl("PPC",SampleID),.N,SampleID] # PPC: exclude the placenta samples sequenced together
  dt.fq.samples[,Type:=ifelse(grepl("^GS-B",SampleID),"term","preterm")] # GS-B-xx: term 
  dt.fq.samples[Type=="preterm",IlluminaID:=tstrsplit(SampleID,"-",keep=3L)]
  dt.fq.samples[Type=="term",IlluminaID:=tstrsplit(SampleID,"-",keep=4L)] # n=755 samples
  dt.fq.samples[,.N,.(Type,IlluminaID)][order(-N)] # three samples sequenced >1 (therefore, 755-3 unique samples)

  (dt.samples<-merge(dt.fq.samples, # n=755 samples
                     dt.meta.info[,.N,.(Type,IlluminaID,PI,setnumber,Batch,pn_female,GA,case)] # n=752 samples
                     ,by=c("Type","IlluminaID")
  )) # n=755
  dt.samples[,Sex:=ifelse(pn_female,"Female","Male")]
  dt.samples[,GA:=gsub(" weeks","wk",GA)]
  dt.samples[,Condition:=ifelse(case==1,"Case","Control")]
  dt.samples[,Batch:=gsub(" ","_",Batch)]
  dt.samples[,Batch:=gsub("&","",Batch)]

  # make as `factor` to be used in deseq design
  dt.samples$Batch<-factor(dt.samples$Batch)
  dt.samples$Condition<-factor(dt.samples$Condition,levels=c("Control","Case"))
  dt.samples$Sex<-factor(dt.samples$Sex,levels=c("Male","Female")) # F-M
  dt.samples$GA<-factor(dt.samples$GA)

  #save(dt.samples, file="RData/dt.samples.RData")
  mat.samples<-data.frame(dt.samples)
  rownames(mat.samples)<-dt.samples$SampleID
}
message("dt.samples loaded")


##
## Sample QC
##
if(F){
  xtabs(~Type+Condition, data=dt.samples) %>% addmargins

  ###########
  # samples #
  ###########
  dt.meta.info[,.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)][order(-N)] # n=752 samples
  xtabs(~N+Type, data=dt.meta.info[,.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)]) %>% addmargins # N for the number of tubes

  xtabs(~Type+GA, data=dt.meta.info[,.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)]) %>% addmargins # n=752 unique samples
  xtabs(~case+GA, data=dt.meta.info[Type=="preterm",.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)]) %>% addmargins # n=277 unique samples in total
  xtabs(~case+GA, data=dt.meta.info[Type=="term",.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)]) %>% addmargins # n=475 unique samples in total

  dt.meta.info[,.N,.(Type,IlluminaID)][order(IlluminaID)] # n=752 samples (Type & IlluminaID are unique for a sample)
  dt.meta.info[,.N,.(Type,IlluminaID)][,.N,Type]

  ##################
  # No. of patient #
  ##################
  dt.meta.info[,.N,.(Type,PI,setnumber)][order(N)] #n= 195 patients #N=8 (two tubes per GA per patient); N=4 (one tube per GA per patient)
                                                  # Type, PI, setnumber are unique for a patient
  dt.meta.info[,.N,.(Type,PI,setnumber,pn_female,case)][order(N)] #n= 195 patients #N=8 (two tubes per GA per patient); N=4 (one tube per GA per patient)
  dt.meta.info[,.N,.(Type,PI,setnumber,pn_female,case)][,.N,pn_female]
  dt.meta.info[,.N,.(Type,PI,setnumber,pn_female,case)][,.N,.(Type,pn_female)]
  dt.meta.info[,.N,.(Type,PI,setnumber,pn_female,case)][,.N,.(Type,pn_female,case)]
  dt.meta.info[,.N,.(Type,PI,setnumber,case)][,.N,case]
  dt.meta.info[,.N,.(Type,PI,setnumber,case)][,.N,.(Type,case)]

  # find samples seuqneced more than once
  dt.samples[,.N,.(Type,IlluminaID,PI,setnumber,pn_female,GA,case)][order(-N)] # n=752 samples
  dt.samples[,.N,.(Type,IlluminaID)][N>1] # 3 samples sequenced >1
  dt.samples[grepl("-b$",SampleID)] # 3 samples sequecned >1
  dt.samples[(Type=="preterm" & IlluminaID %in% c("CX","HQ")) | (Type=="term" & IlluminaID=="UW")] # sequenced >1

  dt.samples[,.N,.(Type,PI,setnumber)][order(N)] #n= 195 patients
  dt.samples[!grepl("-b$",SampleID),.N,.(Type,PI,setnumber)][N<4][order(N)] # 23 patients with less than 4 plasma samples 
  dt.samples[Type=="preterm" & PI=="MP"]

##
## QC-failed
##
#   preterm:
#   - GS-59-CX is data from the same plasma sample as GS-59-CX-b. Both failed QC
#   - GS-179-HQ is data from the same plasma sample as GS-179-HQ-b. Only GS-179-HQ-b passed QC
  dt.colDataAll[grepl("-b$",names)] # CX (CX-b): both failed; HQ (HQ-b): only HQ-b passed QC according to illumina
  dt.colData<-dt.colData[!names %in% c("GS-59-CX-b","GS-179-HQ")] # remove these two samples (n=277)
#   term:
#   - GS-B-374-UW is data from the same plasma sample as GS-B-374-b. Both failed QC
  dt.colDataTerm<-dt.colDataAll[Type=="term" & !names %in% c("GS-B-374-UW","GS-B-374-UW-b")] # remove these two samples as both of them flagged as "failed" by illumina (it is a 28wk control sample) n=474

}

##
## samples.csv for EGA
##
if(F){
  dt.pops<-  rbind(
  readxl::read_excel("data/cfRNA_samples_preterm_POPSID.xlsx") %>% data.table,
  readxl::read_excel("data/cfRNA_samples_term_POPSID.xlsx") %>% data.table,
  fill=T
  )
  dt.pops.samples<-merge(dt.meta.info, dt.pops[,.(Sample_number,POPSID)])
  dt.pops.samples$POPSID<-as.integer(dt.pops.samples$POPSID)
  dt.pops.samples[order(POPSID)]
  dt.pops.samples[,.N,.(Type,POPSID)][order(-N)] # n=195 patients
  dt.pops.samples[POPSID==682]

  merge(dt.pops.samples[POPSID==682], dt.samples[,.(Type,IlluminaID,SampleID)], by=c("Type","IlluminaID"))
  merge(dt.pops.samples[POPSID==682], dt.samples[,.(Type,IlluminaID,SampleID,Sex,Condition)], by=c("Type","IlluminaID"))[,.N,.(SampleID,GA,POPSID,Type,FetalSex=Sex,Condition)]

  # alias,title,description,biological_sex,subject_id,phenotype,biosample_id,case_control,organism_part,cell_line
  (dt.foo<-merge(dt.pops.samples, dt.samples, by=c("Type","IlluminaID"))[,.N,.(SampleID,POPSID,Condition,gestational_age=GA.x,Type,fetal_sex=Sex)][,.(alias=SampleID,title="Maternal plasma",description="Maternal plasma",biological_sex="female",subject_id=POPSID,phenotype=ifelse(Condition=="Case","Pre-eclampsia with fetal growth restriction","control"),biosample_id=NA,case_control=ifelse(Condition=="Case","case","control"),organism_part="Human",cell_line=NA,gestational_age,fetal_sex)])

  fwrite(dt.foo,file="data/samples.csv")
}


##########################################
# DRAGEN stat (from Suzanne of illumina) #
##########################################
#dt.dragen<-readxl::read_excel("data/POPS_SEQUENCING_QC_METRICS.xlsx") %>% data.table

##############
## Lib Size ##
##############
if(F){
dt.reads<-fread("data/SLX-ILM-Plasma2021.fq.read.base.cnt.txt")
dt.reads[,SampleID:=tstrsplit(V1,".",fixed=T,keep=2L)]
dt.reads[grepl("^GS-PPC",SampleID),.(NumRead=sum(V2),NumFile=.N),SampleID][order(NumRead)]
dt.reads[!grepl("^GS-PPC",SampleID),.(NumRead=sum(V2),NumFile=.N),SampleID][order(NumRead)]
}

#############################################################
## Tx & gene meta info from FASTA used for salmon indexing ##
#############################################################
my.salmon.index="POPS-2022.GRCh38.88"
if(my.salmon.index=="POPS-2022.GRCh38.88"){
  dt.tx2gene<-fread("~/results/RNA-Seq/Placentome/gffcompare/POPS-2022/POPS-Placenta-Transcriptome/POPS-2022.GRCh38.88.Novel.Known.Freq.0.1.TPM.0.1.tr.reconstruction.tx2gene.txt", header=F,col.names=c("transcript_id","gene_id"))
}else if(my.salmon.index=="CHESS-POPS-2022.GRCh38.88"){
  dt.tx2gene<-fread("~/data/Annotation/CHESS_POPS/CHESS2.2.POPS.GRCh38.88.Novel.Freq.0.1.TPM.0.1.tx2gene.txt", header=F,col.names=c("transcript_id","gene_id"))
}else{
  my.RData<-"RData/dt.tx2gene.Homo_sapiens.GRCh38.88.cdna.all.ncrna.fa.gz.RData"
  if(!file.exists(my.RData)){
  my.cmd<-"zcat data/Homo_sapiens.GRCh38.88.cdna.all.ncrna.fa.gz | grep -P '^>' | cut -d ' ' -f1,3,4,5,6,7 | awk '{split($2,foo,\":\"); print substr($1,2), foo[3], substr($3,6), substr($4,14), substr($5,20), substr($6,13)}'"
  dt.tx2gene<-(tstrsplit(system(my.cmd, intern=T), " ") %>% simplify2array %>% data.table)
  my.colNames<-c("transcript_id","chromosome_name","gene_id","transcript_biotype","gene_biotype","hgnc_symbol")
  setnames(dt.tx2gene,my.colNames)
    #dt.tx2gene[,.N,chromosome_name] %>% nrow # 380
    #dt.tx2gene[,.N,gene_id] %>% nrow # 62803 gene
    #dt.tx2gene[,.N,transcript_id] %>% nrow # 217082 tr 
    #dt.tx2gene[!chromosome_name %in% c(seq(1:22),"X","Y","MT"),.N,chromosome_name][order(-N)]
    message("saving dt.tx2gene...")
    save(dt.tx2gene,file="RData/dt.tx2gene.Homo_sapiens.GRCh38.88.cdna.all.ncrna.fa.gz.RData")
  }else{
    load("RData/dt.tx2gene.Homo_sapiens.GRCh38.88.cdna.all.ncrna.fa.gz.RData")
    message("dt.tx2gene loaded")
  }
}

######################
## Downsample reads ##
######################
prep_salmon<-function(my.colData,my.salmon="Salmon"){
  li.salmon<-list()
  stopifnot(my.salmon %in% c("Salmon","Salmon_aln","Salmon_dedup"))

  li.salmon[[my.salmon]][["coldata"]]<-my.colData
  li.salmon[[my.salmon]][["se"]]<-tximeta(li.salmon[[my.salmon]][["coldata"]]) # tx level [todo: parallel?]

  if(my.salmon=="Salmon"){
    li.salmon[["Salmon"]][["gse"]]<-summarizeToGene(li.salmon[["Salmon"]][["se"]])
    # gene level
    #li.salmon[[my.salmon]][["counts"]]<-assays(li.salmon[[my.salmon]][["se.g"]])[["counts"]]
    #li.salmon[[my.salmon]][["TPM"]]<-assays(li.salmon[[my.salmon]][["se.g"]])[["abundance"]] # TPM
  }
  if(my.salmon=="Salmon_aln"){
    li.salmon[["Salmon_aln"]][["gse"]]<-tximeta(li.salmon[["Salmon_aln"]][["coldata"]],skipMeta=T,tx2gene=dt.tx2gene[,.(transcript_id,gene_id)],txOut=F) 
    # remove the trailing version number (e.g. ENSG.1 => ENSG)
    #rownames(li.salmon[["Salmon_aln"]][["counts"]])<-li.salmon[["Salmon_aln"]][["counts"]] %>% rownames %>% substr(1,15)
    #rownames(li.salmon[["Salmon_aln"]][["TPM"]])<-li.salmon[["Salmon_aln"]][["TPM"]] %>% rownames %>% substr(1,15)
  }

  return(li.salmon)
} # end of function

# downsampled to 1M, 5M and 10M reads
# only 100 samples
# to see which ones better: 1) `Salmon` (i.e. Salmon SA mode) or 2) `Salmon_aln` (HiSat2+Salmon)
prep_ds_salmon<-function(my.slx,my.salmon.runs=c("Salmon","Salmon_aln")){
  li.TPM<-list()

  # Salmon & Salmon_aln
  for(my.salmon in my.salmon.runs){
    message(my.salmon)
    dt.foo<-data.table(`files`=system(paste0("ls ", "~/results/",my.slx,"/",my.salmon,"/GRCh38.88/*/quant.genes.sf"), intern=T))
    dt.foo[,names:=tstrsplit(files,"/",keep=8)]

    li.quant.genes<-lapply(dt.foo$files, function(i) fread(i) %>% as.matrix(rownames="Name"))

    li.TPM[[my.salmon]][["TPM"]]<-sapply(li.quant.genes, function(i) i[,"TPM",drop=F])
    li.TPM[[my.salmon]][["counts"]]<-sapply(li.quant.genes, function(i) i[,"NumReads",drop=F])

    colnames(li.TPM[[my.salmon]][["TPM"]])<-dt.foo$names
    colnames(li.TPM[[my.salmon]][["counts"]])<-dt.foo$names

    # remove the trailing version number (e.g. ENSG.1 => ENSG)
    rownames(li.TPM[[my.salmon]][["TPM"]])<-li.quant.genes[[1]] %>% rownames %>% substr(1,15)
    rownames(li.TPM[[my.salmon]][["counts"]])<-li.quant.genes[[1]] %>% rownames %>% substr(1,15)
  }

  return(li.TPM)
} # end of function


##
## Parse the count data # 
##

init_li_tpm<-function(){
  my.slx<-"SLX-ILM-Plasma2021.Homo_sapiens.v1"
  my.RData<-paste0("RData/li.TPM.RData")
  if(!file.exists(my.RData)){
    li.TPM<-prep_ds_salmon(my.slx,my.salmon.runs="Salmon") # NB, this run is only for `Salmon (SA)`
    save(li.TPM, file=my.RData)
    message("li.TPM saved")
  }else{
    system.time(attach(my.RData))
    message("li.TPM loaded")
  }
}

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

# the above using the `boot` package
# see https://www.statmethods.net/advstats/bootstrapping.html
get_LPOCV_boot2<-function(x,index,fit){
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

            #fit<-glm(y~. , data = my.data2, family = "binomial") # fit the model based on the remaining
            predict.glm(fit, newdata=my.data1, type="response") # predict the outcome of the pair using the model above
            #roc(response=my.data2$y, predictor=fitted(fit))$auc
        })

        # the proportion of all pairwise combinations in which the predicted probability was greater for the case than for the control 
        c.stat<-apply(out,2, function(i) i[1] > i[2]) %>% table 
        return(c.stat["TRUE"] / sum(c.stat) *100)
}

# 
# x; dataset
# my.fold: an arbitrary fold name
# my.num: number of desired features
# my.func: a functioin (isa `list`) to be used by rfeControl 
# replaced with `runRFE2`
get_rfe_rank<-function(x,my.fold,my.num=4, my.func=rfFuncs){
  #x isa `matrix` and should contain 'y' column
  all.features<-colnames(x)[colnames(x)!="y"]

  ##########
  # run RF #
  ##########
  cl <- makePSOCKcluster(8) # No. of cores to use
  registerDoParallel(cl)
  set.seed(333) # set a random seed for a reproducibility
  cv.rfe<-caret::rfe(
             x=as.data.frame(x[,all.features]),
             y=factor(ifelse(x[,"y"]==1,'case','non_case'),levels=c("non_case","case")),
             rfeControl=rfeControl(functions=my.funcs, # lmFuncs, rfFuncs, treebagFuncs, nbFuncs
                                   method="cv", 
                                   number=10, # 10 by default for cv
                                   saveDetails=T,
                                   verbose=T),
  )
  stopCluster(cl)

  cv.rfe
  predictors(cv.rfe)
  cv.rfe$optVariables
  cv.rfe$optsize
  cv.rfe$bestSubset
  varImp(cv.rfe) # not the same order of the above
  varImp(cv.rfe, useModel=F, nonpara=F, scale=T)

  cv.rfe$fit  # isa `randomForest`

  dt.foo<-data.table(cv.rfe$variable)
  dt.foo[,.N,.(Resample)]
  dt.foo[,.N,.(Resample,Variables)]
  dt.foo[Resample=="Fold01"]

  data.table(`fold`=my.fold,`feature`=cv.rfe$optVariables[1:my.num],`rank`=1:my.num)
} # end of get_rfe_rank

