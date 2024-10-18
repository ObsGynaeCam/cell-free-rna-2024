dt.foo<-fread("~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/SLX-ILM-Plasma2021.fq.read.base.cnt.txt", col.names=c("FQ","Read","Base"))
#dt.foo[,Type:=ifelse(grepl("PPC",FQ),"Placenta",ifelse(grepl("GS-B",FQ),"term","preterm"))]
dt.foo[,SampleID:=tstrsplit(FQ,".",fixed=T,keep=2L)]
dt.foo[,FR:=tstrsplit(FQ,".",fixed=T,keep=5L)]


dt.fq<-merge(dt.foo,dt.samples[,.(SampleID,Type)])

## term
dt.fq[Type=="term",.N,SampleID][,.N,N]
dt.fq[Type=="term",.N,SampleID][order(-N)]

dt.fq[Type=="term",.N,SampleID]
dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-Sum)]
dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)]
dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][,.N,N]

N8_sample<-dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==8]$SampleID
N10_sample<-dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==10]$SampleID
N16_sample<-dt.fq[Type=="term",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==16]$SampleID

# 4 files per sample (4 FW + 4 RV; 8 files total)
dt.fq[SampleID %in% N8_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/4 # 5M
dt.fq[SampleID %in% N8_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/4]

# 5 files per sample (5 FR + 5 RV; 10 files total)
dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/5 # 4M
dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/5]
bad.fq<-dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/5]$FQ

# 8 files per sample (8 FR + 8 RV; 16 files total)
dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/8 # 2.5M
dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/8]
bad.fq2<-dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/8]$FQ

dt.fq[SampleID %in% N16_sample,.N,SampleID]
dt.fq[SampleID %in% N16_sample & !FQ %in% bad.fq2,.N,SampleID]

##
dt.fq[Type=="term",.(Sum=sum(Read),.N),SampleID][order(Sum)]
dt.bar<-dt.fq[Type=="term" & !FQ %in% c(bad.fq, bad.fq2),.(Sum=sum(Read),.N),SampleID][,N_sample:=20*1e+6/N*2]

merge(dt.fq,dt.bar)
merge(dt.fq,dt.bar)[,.(sum(N_sample)),.(SampleID,FR)]

fwrite(merge(dt.fq,dt.bar), file="~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/term.downsample.FQ.20M.csv")

#for i in `grep fq.gz ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/downsample.FQ.20M.csv`; do FQFILE=`echo $i | cut -d, -f2`; DS=`echo $i | cut -d, -f9`; FNAME=SLX-ILM-Plasma2021-ds20M.`echo $FQFILE | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo  $FNAME; time seqtk sample -s11 $FQFILE $DS | gzip -9 > ~/data/fastq/SLX-ILM-Plasma2021-ds20M/$FNAME; done


## preterm
dt.fq[Type=="preterm",.N,SampleID][,.N,N]
dt.fq[Type=="preterm",.N,SampleID][order(-N)]

dt.fq[Type=="preterm",.N,SampleID]
dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-Sum)]
dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)]
dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][,.N,N]

N8_sample<-dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==8]$SampleID
N10_sample<-dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==10]$SampleID
N16_sample<-dt.fq[Type=="preterm",.(.N,Sum=sum(Read)/1e+6), SampleID][order(-N)][N==16]$SampleID

# 4 files per sample (4 FW + 4 RV; 8 files total)
dt.fq[SampleID %in% N8_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/4 # 5M
dt.fq[SampleID %in% N8_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/4]

# 5 files per sample (5 FR + 5 RV; 10 files total)
dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/5 # 4M
dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/5]
bad.fq<-dt.fq[SampleID %in% N10_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/5]$FQ

# 8 files per sample (8 FR + 8 RV; 16 files total)
dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)]
20*1e+6/8 # 2.5M
dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/8]
bad.fq2<-dt.fq[SampleID %in% N16_sample, .(SampleID,FQ,Read/1e+6)][order(V3)][V3<20/8]$FQ

dt.fq[SampleID %in% N16_sample,.N,SampleID]
dt.fq[SampleID %in% N16_sample & !FQ %in% bad.fq2,.N,SampleID]

##
dt.fq[Type=="preterm",.(Sum=sum(Read),.N),SampleID][order(Sum)]
dt.bar<-dt.fq[Type=="preterm" & !FQ %in% c(bad.fq, bad.fq2),.(Sum=sum(Read),.N),SampleID][,N_sample:=20*1e+6/N*2]

merge(dt.fq,dt.bar)
merge(dt.fq,dt.bar)[,.(sum(N_sample)),.(SampleID)]
merge(dt.fq,dt.bar)[,.(sum(N_sample)),.(SampleID,FR)]

fwrite(merge(dt.fq,dt.bar), file="~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/preterm.downsample.FQ.20M.csv")
