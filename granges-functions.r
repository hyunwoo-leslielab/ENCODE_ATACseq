reduce.granges <- function (gr) {
  vals <- values (gr)
  reduced.gr <- GRanges (seqnames (gr)[1], IRanges (min (start (gr)), max(end (gr))),
                         start.t = min (vals$start.t), stop.t = max (vals$stop.t),
                         sig.value = mean (vals$sig.value), signal.value = mean (vals$signal.value),
                         p.value = mean (vals$p.value), q.value = mean (vals$q.value))
  return (reduced.gr)
}

convert.to.gr <- function (data, sig.value) {
  gr <- GRanges (data[,'chr'], IRanges (data[,'start.ori'], data[,'stop.ori']),
                  start.t=data[,'start'], stop.t=data[,'stop'],
                  sig.value=data[,sig.value], signal.value=data[,'signal.value'],
                  p.value=data[,'p.value'], q.value=data[,'q.value'])
  return (gr)
}


convert.to.bed <- function (gr) {
  vals <- values (gr)
  bed <- data.frame (id=1:length (gr), sig.value=vals$sig.value,
                     start = vals$start.t, stop=vals$stop.t,
                     signal.value = vals$signal.value,
                     p.value = vals$p.value, q.value = vals$q.value,
                     chr = as.character (seqnames (gr)),
                     start.ori = start (gr), stop.ori = end (gr))
  return (bed)
}


merge.peaks.gr <- function (gr1, gr2, p.value.impute) {

  ## Global union
  global.union <- reduce (union (gr1, gr2))

  ## Overlaps and counts
  gr1.overlaps <- as.matrix (findOverlaps (global.union, gr1)); gr1.counts <- countOverlaps (global.union, gr1)
  gr2.overlaps <- as.matrix (findOverlaps (global.union, gr2)); gr2.counts <- countOverlaps (global.union, gr2)
  rownames (gr1.overlaps) <- sprintf ("S%d", gr1.overlaps[,1])
  rownames (gr2.overlaps) <- sprintf ("S%d", gr2.overlaps[,1])

  ## Unique in both
  inds <- sprintf ("S%d", which (gr1.counts == 1 & gr2.counts == 1))
  merge1 <- gr1[gr1.overlaps[inds,2]]; merge2 <- gr2[gr2.overlaps[inds,2]]

  ## Reduce overlaps
  for (ind in c(which (gr1.counts > 1 | gr2.counts > 1))) {
    m1 <- gr1[gr1.overlaps[gr1.overlaps[,1] == ind,2]]
    if (length (m1) > 1)
      m1 <- reduce.granges (m1)

    m2 <- gr2[gr2.overlaps[gr2.overlaps[,1] == ind,2]]
    if (length (m2) > 1)
      m2 <- reduce.granges (m2)
    
    merge1 <- c(merge1, m1)
    merge2 <- c(merge2, m2)
  }

  ## Non overlapping elements
  gr1.non <- gr1[countOverlaps (gr1, gr2) == 0]
  if (length (gr1.non) > 0 ) {
    merge1 <- c(merge1, gr1.non)
    values (gr1.non)$p.value <- values (gr1.non)$q.value <- p.value.impute
    values (gr1.non)$sig.value <- values (gr1.non)$signal.value <- p.value.impute
    merge2 <- c(merge2, gr1.non)
  }
  
  gr2.non <- gr2[countOverlaps (gr2, gr1) == 0]
  if (length (gr2.non) > 0 ){
    merge2 <- c(merge2, gr2.non)
    values (gr2.non)$p.value <- values (gr2.non)$q.value <- p.value.impute
    values (gr2.non)$sig.value <- values (gr2.non)$signal.value <- p.value.impute
    merge1 <- c(merge1, gr2.non)
  }

  ## Sort according to starts
  inds <- sort (values (merge1)$start.t, index.return=TRUE)$ix
  merge1 <- merge1[inds]; merge2 <- merge2[inds]

  return (list (merge1=merge1, merge2=merge2))
  
  return (list (merge1=convert.to.bed (merge1), merge2=convert.to.bed (merge2)))
}


rm.unmatch.gr <- function(sample1, sample2, p.value.impute=0){

  inds <- values (sample1)$sig.value > p.value.impute & values (sample2)$sig.value > p.value.impute
  sample1.prune <- sample1[inds]
  sample2.prune <- sample2[inds]
 
  invisible(list(sample1=sample1.prune, sample2=sample2.prune))
}
  
  
