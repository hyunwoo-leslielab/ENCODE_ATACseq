#!/usr/bin/env Rscript

# prepare for IDR
source("functions-all-clayton-12-13.r")
source("granges-functions.r")

consistency.analysis <- function (peaktable1, peaktable2, half.width=NULL,
                                  overlap.ratio=0, is.broadpeak=FALSE, sig.value, chr.size) {
   
  ## ########### process the data
  ## process data, summit: the representation of the location of summit
  rep1 <- process.narrowpeak(peaktable1, chr.size, half.width=half.width, summit="offset", broadpeak=is.broadpeak)
  rep2 <- process.narrowpeak(peaktable2, chr.size, half.width=half.width, summit="offset", broadpeak=is.broadpeak)

  ## Log information
  show (sprintf ("Peaks from table 1: %d. Peaks after cleaning: %d", nrow (rep1$data.ori), nrow (rep1$data.cleaned)))
  show (sprintf ("Peaks from table 2: %d. Peaks after cleaning: %d", nrow (rep2$data.ori), nrow (rep2$data.cleaned)))

  ## compute correspondence profile (URI)
  show ('Running URI (Step 1 of 3)...')
  # time.start <- get.time ()
  uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, sig.value1=sig.value, sig.value2=sig.value, overlap.ratio=overlap.ratio)
  # time.end <- get.time ()
  # show (sprintf ("Time of determining URI: %.2f", (time.end - time.start)/60))

  ## EM procedure for inference
  show ('Running EM (Step 2 of 3)....')
  # time.start <- get.time ()
  em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)
  # time.end <- get.time ()
  # show (sprintf ("Time of running EM: %.2f", (time.end - time.start)/60))

  ## add on 3-29-10
  ## output both local idr and IDR
  show ('Determining IDR (Step 3 of 3)...')
  # time.start <- get.time ()
  idr.local <- 1-em.output$em.fit$e.z
  IDR <- c()
  o <- order(idr.local)
  IDR[o] <- cumsum(idr.local[o])/c(1:length(o))
  # time.end <- get.time ()
  # show (sprintf ("Time of determining IDR: %.2f", (time.end - time.start)/60))
 
  dp <- em.output$data.pruned
  write.out.data <- data.frame (chr1 = as.character (seqnames (dp$sample1)),
                                start1 = start (dp$sample1),
                                stop1 = end (dp$sample1),
                                sig.value1 = values (dp$sample1)$sig.value,
                                chr2 = as.character (seqnames (dp$sample2)),
                                start2 = start (dp$sample2),
                                stop2 = end (dp$sample2),
                                sig.value2 = values (dp$sample2)$sig.value,
                                idr.local = 1-em.output$em.fit$e.z, IDR=IDR)
                                 
 
  ## Compute marginal means
  mar.mean <- get.mar.mean(em.output$em.fit)

  ## Set up results object
  show ('Packaging results....')
  results <- list ()
  results$URI <- uri.output
  results$EM <- em.output
  results$out.table <- write.out.data
  results$marginal.mean <- mar.mean

  return (results)
}

# Function for calculating macs2 p-values given peak locations and number of reads per peak.
# original function written by Manu in sequencingUtils.
# gsize 2.7e9 corresponds to the mouse genome, roughly.
determine.peaks.significance <- function (readsGR, peaksGR, tag.size = 50,
                                          outer.window = 10000, gsize = 2.7e9,
                                          countDistroPlot, outlierCutoff = NA,
                                          chromSizes) {

     # keep only peaks and reads in desired chromosomes.

     readsGR <- readsGR[seqnames(readsGR) %in% chromSizes$chr]
     peaksGR <- peaksGR[seqnames(peaksGR) %in% chromSizes$chr]

     ## Determine counts at peaks and extended peaks.

     peaksGR$score <- countOverlaps (peaksGR, readsGR, minoverlap = floor(tag.size/2))
     extended.peaks <- resize (peaksGR, outer.window, fix = "center")
     peaksGR$score.outer <- countOverlaps (extended.peaks, readsGR, minoverlap = floor(tag.size/2))
  peaksGR$height = -1
  peaksGR$summit = -1

     # get summits and their heights.

     for (chrom in chromSizes$chr) {
          cov <- coverage (readsGR[seqnames(readsGR) == chrom])
          command <- paste0("cov <- cov$", chrom)
          eval(parse(text = command))
          pksInChr <- peaksGR[seqnames(peaksGR) == chrom]
          v <- Views (cov, start = start(pksInChr), end = end(pksInChr))
           heights <- max(v)
           summits <- which.max(v)
           summits <- summits - start(pksInChr)
           peaksGR[names(pksInChr)]$height <- heights
           peaksGR[names(pksInChr)]$summit <- summits
     }

     ## Local lambda.
    
     total.tag.count <- length(readsGR)
     effective.widths <- pmax (tag.size * 2, width (peaksGR))
     lambda.local <- effective.widths * total.tag.count / gsize
     # => seems to me like the number of reads the peak would get given its size
     # if the reads were spread uniformly across the genome.

     ## Background lambda.

     lambda.outer <- peaksGR$score.outer/outer.window * effective.widths
     # => hmm ... this one looks to me more like the local lambda.

     lambda.selected <- pmax (lambda.local, lambda.outer)

     ## Calculate pvalues.
    
     peaksGR$pvalue <- round (-10 * log10 (apply (cbind (lambda.selected, peaksGR$score * effective.widths / width (peaksGR)), 1,
                                                           function (x) { ppois (x[2], x[1], lower.tail = FALSE) })), 2)

     cat(paste0(length(peaksGR[is.infinite(peaksGR$pvalue)]), " peaks with inf pvalue.\n"))

     ## Remove infinities.

     peaksGR$pvalue[is.infinite (peaksGR$pvalue)] <- max(peaksGR[is.finite(peaksGR$pvalue)]$pvalue)

     return (peaksGR)
}

# create atlas.

# Inspired in EncodeEpigenomics/merged_peaks/merge_cell_type_peaks.R.

collapse.pairwise.celltype.peaks <- function (peaks1, peaks2, overlap.ratio.cutoff=0.75) {

  ## Function to find overlapping ratios
  find.overlap.ratio <- function (gr1, gr2) {
   
    ## Find Overlaps
    overlaps <- as.matrix (findOverlaps (gr1, gr2))

    ## Build ranges
    ranges <- cbind (start (gr1)[overlaps[,1]], end (gr1)[overlaps[,1]],
                     start (gr2)[overlaps[,2]], end (gr2)[overlaps[,2]])
    ranges <- t(apply (ranges, 1, sort))

    ## Min widths
    widths <- pmin (width (gr1)[overlaps[,1]], width (gr2)[overlaps[,2]])

    ## Overlap ratios
    overlap.ratio <- (ranges[,3] - ranges[,2])/widths

    ## Find best in both 
    best.gr1 <- tapply (1:nrow (overlaps), overlaps[,1], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
    best.gr2 <- tapply (1:nrow (overlaps), overlaps[,2], function (x) { x[overlap.ratio[x]==max(overlap.ratio[x])][1] } )
    common <- intersect (best.gr1, best.gr2)
   
    return (list (overlap.ratio = overlap.ratio,
                  ranges = ranges,
                  overlaps = overlaps,
                  best.common = common))
  }

  ## Reset metadata
  values (peaks1) <- values (peaks2) <- NULL

  ## Overlapping peaks which exceed overlap.ratio.cutoff
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common[or$overlap.ratio[or$best.common] > overlap.ratio.cutoff]
  ## Create GRanges object with common regions
  union <- GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                    IRanges (or$ranges[common,2], or$ranges[common,3]))

  ## Determine disjoint peaks
  peaks1 <- peaks1[countOverlaps (peaks1, union) == 0]
  peaks2 <- peaks2[countOverlaps (peaks2, union) == 0]
 
  ## Other overlapping peaks
  or <- find.overlap.ratio (peaks1, peaks2)
  common <- or$best.common
  ## Create GRanges with union of the regions
  union <- c(union,GRanges (seqnames (peaks1)[or$overlaps[common,1]],
                            IRanges (or$ranges[common,1], or$ranges[common,4])))

  ## Non overlapping peaks
  union <- c(union, peaks1[countOverlaps (peaks1, union) == 0])
  union <- c(union, peaks2[countOverlaps (peaks2, union) == 0])

  return (union)
}

# assumes function "collapse.pairwise.celltype.peaks".
mergeCTpeaks <- function(pkLst, GRfile, DTfile) {

    combined.peaks <- collapse.pairwise.celltype.peaks (pkLst[[1]], pkLst[[2]])
    if (length(pkLst) > 2) {
      for (pkCtr in 3:length(pkLst)) {
        combined.peaks <- collapse.pairwise.celltype.peaks (combined.peaks, pkLst[[pkCtr]])
      }
    }

    ## annotate cell type where peak was called.

    for (cellType in names(pkLst)) {

        overlaps <- countOverlaps (combined.peaks, pkLst[[cellType]])
        mcols (combined.peaks)[,paste0(cellType, ".peak")] <- 0
        mcols (combined.peaks)[,paste0(cellType, ".peak")][overlaps > 0] <- 1

    }

    names(combined.peaks) <- 1:length(combined.peaks)
    combined.peaks$pattern <- apply(as.matrix(mcols(combined.peaks)), 1, function(row) {
                                                                                argsToPaste0 <- as.list(row)
                                                                                do.call(paste0, argsToPaste0)
                                                                            })

    combined.peaks.DT <- data.table(pk = names(combined.peaks),
                                    chr = as.character(seqnames(combined.peaks)),
                                    start = as.integer(start(combined.peaks)),
                                    end = as.integer(end(combined.peaks)),
                                    pattern = combined.peaks$pattern)

    saveRDS(combined.peaks, file = GRfile)
    saveRDS(combined.peaks.DT, file = DTfile)

}
