#!/usr/bin/env Rscript
library(data.table)
library(GenomicRanges)
peaksRaw = lapply(dir(pattern="*_peaks.narrowPeak"), function(a) fread(a))
peaksRaw = lapply(peaksRaw, function(a) a[, list(chr=V1, start=V2, stop=V3, p.value=V8, signal.value=V7, q.value=-1)])
peaksSummits = lapply(dir(pattern="*_summits.bed"), function(a) fread(a))
peaksRaw = lapply(1:length(peaksRaw), function(a) data.table(peaksRaw[[a]], summit=peaksSummits[[a]]$V2 + 1))
peaks.GR = lapply(peaksRaw, function(a) GRanges(seqnames=a$chr, ranges=IRanges(start=a$start + 1, end=a$stop), strand="*"))
peaks.GR = lapply(peaks.GR, function(a) {names(a) = 1:length(a); a})
fnrep = dir(pattern="*.bed.gr.rds", recursive=T)
fn = gsub("_peaks.narrowPeak", "", dir(pattern="*_peaks.narrowPeak"))
readsRep.GR = list()
for (j in 1:length(fn)) {
  readsRep.GR[[fn[j]]] = list("rep1"=readRDS(fnrep[2 * j - 1]), "rep2"=readRDS(fnrep[2 * j]))
}
names(readsRep.GR) = fn
mm10chromSizes = fread("mm10.chrom.sizes")
mm10chromSizes = mm10chromSizes[, list(chr = V1, size = V2)]
mm10chromSizes = mm10chromSizes[chr %in% c(paste0("chr", 1:19), "chrX", "chrY", "chrM")]

source("/ifs/e63data/leslielab/hyc2004/programs/ATACseqcode/prepareforIDR.R")

library(gtools)
# Functions for IDR

library(foreach)
library(doMC)
registerDoMC(length(fnrep))

wrapper1 = function(myReadsGR, myPeaksGR, myChromSizes, outer.window=10000, gsize=2.7e9) {
  return(determine.peaks.significance(readsGR=myReadsGR, peaksGR=myPeaksGR, tag.size=min(width(myReadsGR)), outer.window=outer.window, gsize=2.7e9, chromSizes=myChromSizes))
}

wrapper2 = function(a, b, myChromSizes, half.width=NULL, overlap.ratio=0, is.broadpeak=F, sig.value="p.value") {
  return(consistency.analysis(
    peaktable1=as.data.frame(data.table(
      chr=as.character(seqnames(a)), start=start(a), stop=end(a),
      p.value=a$pvalue, summit=a$summit, signal.value=a$score, q.value=-1)),
    peaktable2=as.data.frame(data.table(
      chr=as.character(seqnames(b)), start=start(b), stop=end(b),
      p.value=b$pvalue, summit=b$summit, signal.value=b$score, q.value=-1)),
    half.width=half.width, overlap.ratio=overlap.ratio, is.broadpeak=is.broadpeak, sig.value=sig.value, chr.size=as.data.frame(myChromSizes))$out.table)
}

peaksWithPvals = foreach (j = 1:length(peaksRaw)) %:%
  foreach (k = 1:(length(fnrep) / length(peaksRaw))) %dopar% {
    wrapper1(readsRep.GR[[j]][[k]], peaks.GR[[j]], mm10chromSizes)
  }

rDF = foreach (j = 1:length(peaksRaw)) %:%
  foreach (pair = combn(1:(length(fnrep) / length(peaksRaw)), 2, simplify=F)) %dopar% {
    wrapper2(peaksWithPvals[[j]][[pair[1]]], peaksWithPvals[[j]][[pair[2]]], mm10chromSizes)
  }

rDT = list(as.data.table(rDF[[1]][[1]]), as.data.table(rDF[[2]][[1]]), as.data.table(rDF[[3]][[1]]))
names(rDT) = fn
saveRDS(rDT, file="rDT.rds")

idr.cutoff = 5e-3
idr.results = lapply(rDT, function(a) GRanges(a$chr1, IRanges(a$start1, a$stop1), IDR=a$IDR))
Pks = lapply(idr.results, function(a) subset(a, IDR <= idr.cutoff))

library(rtracklayer)
# wget -c https://www.encodeproject.org/files/ENCFF790DJT/@@download/ENCFF790DJT.bed.gz; gunzip ENCFF790DJT.bed.gz
blacklist.mm10 = import.bed("ENCFF790DJT.bed")
indices.in.DJT = lapply(Pks, function(a) unique(findOverlaps(a, blacklist.mm10)@queryHits))

Pks.noblacklist = lapply(1:length(peaksRaw), function(a) Pks[[a]][-indices.in.DJT[[a]]])
names(Pks.noblacklist) = fn
for (j in 1:length(fn)) {
  write.table(data.frame(chrom=seqnames(Pks.noblacklist[[j]]), chromStart=start(Pks.noblacklist[[j]]) - 1, chromEnd=end(Pks.noblacklist[[j]]), name=paste0(fn[j], "_", 1:length(Pks.noblacklist[[j]])), score=Pks.noblacklist[[j]]$IDR, strand="."), file=paste0(fn[j], "_IDR0.005.bed"), sep="\t", row.names=F, quote=F, col.names=F)
}

mergeCTpeaks(pkLst=Pks.noblacklist, GRfile="ATACpeakatlasGR.rds", DTfile="ATACpeakatlasDT.rds")
combined.peaks = readRDS("ATACpeakatlasGR.rds")
write.table(data.frame(chrom=seqnames(combined.peaks), chromStart=start(combined.peaks) - 1, chromEnd=end(combined.peaks), name=paste0("Atlas_", 1:length(combined.peaks)), score=1, strand="."), file="ATAC-seq-atlas.bed", sep="\t", row.names=F, quote=F, col.names=F)
#write.table(data.frame(chrom=seqnames(combined.peaks), chromStart=start(combined.peaks) - 1, chromEnd=end(combined.peaks), name=paste0("Atlas_", 1:length(combined.peaks)), score=1, strand="."), file=gzfile("ATAC-seq-atlas.bed.gz"), sep="\t", row.names=F, quote=F, col.names=F)
