###########################
# Random sampling for 100 #
###########################
for i in `ls ~/data/fastq/SLX-ILM-Plasma2021/*r_1.fq.gz`; do BARCODE=`echo $i | cut -d/ -f7 | cut -d. -f2`; echo $BARCODE;done | grep -v PPC | sort | uniq -c | awk '$1==4{print $0}' | sort -R | head -n 100 | awk '{print $2}' > ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.txt
join -t "," <(sort ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.txt) <(sort -t "," -k1,1 ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/FASTQ.list.r1.txt) > ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r1.txt

######################################
# Downsampling to 1M for 100 samples #
######################################
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r1.txt`; do FNAME=SLX-ILM-Plasma2021-ds1M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 250000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds1M/$FNAME; done
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r2.txt`; do FNAME=SLX-ILM-Plasma2021-ds1M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 250000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds1M/$FNAME; done

######################################
# Downsampling to 5M for 100 samples #
######################################
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r1.txt`; do FNAME=SLX-ILM-Plasma2021-ds5M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 1250000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds5M/$FNAME; done
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r2.txt`; do FNAME=SLX-ILM-Plasma2021-ds5M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 1250000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds5M/$FNAME; done

#######################################
# Downsampling to 10M for 100 samples #
###################]###################
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r1.txt`; do FNAME=SLX-ILM-Plasma2021-ds10M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 2500000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds10M/$FNAME; done
 for i in `awk -F"," '{print $2}' ~/results/SLX-ILM-Plasma2021.Homo_sapiens.v1/random.sample.100.FASTQ.list.r2.txt`; do FNAME=SLX-ILM-Plasma2021-ds10M.`echo $i | cut -d/ -f7 | cut -d. -f2,3,4,5,6,7`; echo $i; time seqtk sample -s11 $i 2500000 | gzip -9 >  ~/data/fastq/SLX-ILM-Plasma2021-ds10M/$FNAME; done
