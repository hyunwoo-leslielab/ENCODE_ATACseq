#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
input1 = args[1]
output1 = args[2]
exptype = args[3]

if (exptype == "ChIP" & !is.null(args[4])) {
  shiftsize = as.numeric(args[4])
  shiftnum = ceiling(c("pos"=shiftsize / 2, "neg"=-shiftsize / 2))
}

library(data.table)
library(GenomicRanges)

readsRep = fread(sprintf("zcat %s", input1))
# validate the width of .bed file
print(dim(readsRep))
readsRep = readsRep[readsRep$V3 > readsRep$V2, ]
if (exptype == "ChIP") {
  readsRep.GR = GRanges(seqnames=readsRep$V1, ranges=IRanges(
    start=readsRep$V2 + 1 + (readsRep$V6 == "+") * shiftnum["pos"] + (readsRep$V6 == "-") * shiftnum["neg"],
    end=readsRep$V3 + (readsRep$V6 == "+") * shiftnum["pos"] + (readsRep$V6 == "-") * shiftnum["neg"]),
    strand=readsRep$V6)
} else if (exptype == "ATAC") {
  readsRep.GR = GRanges(seqnames=readsRep$V1, ranges=IRanges(start=readsRep$V2 + 1, end=readsRep$V3), strand=readsRep$V6)
}

saveRDS(readsRep.GR, file=output1)

