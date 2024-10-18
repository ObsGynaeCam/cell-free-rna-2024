$HOME/Install/SalmonTools/scripts/generateDecoyTranscriptome.sh \
  -j 32 \
  -m $HOME/Install/mashmap/mashmap \
  -g $HOME/data/genome/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
  -a $HOME/results/RNA-Seq/Placentome/gffcompare/POPS-2022/POPS-Placenta-Transcriptome/POPS-2022.GRCh38.88.Novel.Known.Freq.0.1.TPM.0.1.tr.reconstruction.gffread.gtf \
  -t $HOME/results/RNA-Seq/Placentome/gffcompare/POPS-2022/POPS-Placenta-Transcriptome/POPS-2022.GRCh38.88.Novel.Known.Freq.0.1.TPM.0.1.tr.reconstruction.gffread.gtf.fa \
  -o $HOME/data/Salmon/decoy/Homo_sapiens/Ensembl/GRCh38/POPS-2022.GRCh38.88.Novel.Known.Freq.0.1.TPM.0.1.tr.reconstruction
