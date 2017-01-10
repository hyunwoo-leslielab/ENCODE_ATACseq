# Convert .bam to .bed and cut 9-bp long artifact for Tn5 transposase
bamToBed -i ENCFF053CGD.bam | adjustBedTn5.sh | gzip > ENCFF053CGD.Tn5.bed.gz
# Save .bed to GRanges Bioconductor object
Rscript convertbedtogr.R ENCFF053CGD.Tn5.bed.gz ENCFF053CGD.bed.gr.rds ATAC
# Pool replicates and transform .bed files to ATAC-specific macs2 running protocol
gunzip -c ENCFF053CGD.Tn5.bed.gz ENCFF557JSA.Tn5.bed.gz | sortBed -i stdin | slopBed -i stdin -g mm10.chrom.sizes -l 75 -r -75 -s | gzip > Sample_rTr.Tn5.sorted.75shifted.bed.gz
# ATAC-specific macs2 running code
macs2 callpeak --tempdir . -t Sample_rTr.Tn5.sorted.75shifted.bed.gz -f BED -n Sample_rTr -g mm -p 1e-2 --nomodel --shift 75 --extsize 200 --keep-dup all --call-summits
# The output is uploaded as a .narrowPeak file
