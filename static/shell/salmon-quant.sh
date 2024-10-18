salmon quant -p 32 \
             -i POPS_TR_INDEX  \
             -l A \
             -1 S1_FQ_1.fq \
             -2 S1_FQ_2.fq \
             --seqBias \
             --gcBias \
             --posBias \
             --discardOrphanQuasi \
             --writeUnmappedNames \
             --writeMapping \
             -o OUTPUT_DIR | samtools view -bS > OUTPUT_DIR/S1_salmon.bam
