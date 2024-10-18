# set parameters
my.disease="PE"
my.type="preterm"
my.salmon="Salmon"
my.salmon.index="POPS-2022.GRCh38.88"
my.slx<-"SLX-ILM-Plasma2021.Homo_sapiens.v1"

# set filter
minCPM=0.1; minRead=10; minFreq=0.1; minFC=1.2

# which p.adjust methods?
adjust.methods=c(`Benjamini & Yekutieli`="BY",
            `Benjamini & Hochberg`="BH",
            `Bonferroni`="bonferroni")
