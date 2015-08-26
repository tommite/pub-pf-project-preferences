library(plyr)
library(hitandrun)

opt.weights <- function(nd, samples, proj.perfs) {
  nd.perfs <- nd %*% as.matrix(proj.perfs)  
  w.sums <- aaply(samples, 1, '%*%', t(nd.perfs))
  aaply(w.sums, 1, which.max)
}

comp.acc <- function(nd, samples, proj.perfs) {
  w.opts <- opt.weights(nd, samples, proj.perfs)

  share.w.ok <- function(proj) {
    proj.pfs <- which(nd[,proj] > 0)
    w.ok <- w.opts %in% proj.pfs
    sum(w.ok) / nrow(samples)
  }
  
  laply(1:ncol(nd), share.w.ok)
}

estimate.ci <- function(nd, samples, proj.perfs) {
  w.opts <- opt.weights(nd, samples, proj.perfs)
  
  popt.pfs <- which(1:nrow(nd) %in% w.opts)
  
  ci <- function(proj) {
    proj.pfs <- which(nd[,proj] > 0)
    sum(proj.pfs %in% popt.pfs) / length(popt.pfs)
  }
  
  list(ci=laply(1:ncol(nd), ci), popt.pfs=popt.pfs)
}

sample.weights.with.pp <- function(proj.perfs, nd, ok.popts, w.constr, n.samples=1E4) {
  state <- har.init(w.constr)

  final.samples <- matrix(nrow=0, ncol=n)

  nd.perfs <- nd %*% as.matrix(proj.perfs)

  while (nrow(final.samples) < n.samples) {
    result <- har.run(state, n.samples-nrow(final.samples))
    samples <- result$samples

    w.sums <- aaply(samples, 1, '%*%', t(nd.perfs), .drop=FALSE)
    w.opts <- aaply(w.sums, 1, which.max)
    w.ok <- aaply(w.opts, 1, '%in%', ok.popts)

    final.samples <- rbind(final.samples, samples[w.ok,])

    state <- result$state
  }
  final.samples
}

analysis.pp <- function(proj.perfs, res.prev, x.in, x.out, w.constr) {
  ok.nds.in <- as.numeric(which(rowSums(res.prev$nd[,x.in, drop=FALSE]) == length(x.in)))
  ok.nds.out <- as.numeric(which(rowSums(res.prev$nd[,x.out, drop=FALSE]) == 0))
  ok.nds <- intersect(ok.nds.in, ok.nds.out)
  ok.popts <- intersect(res.prev$popt.pfs, ok.nds)

  if(length(ok.popts) == 0) {
    error('No portfolios fulfilling project portfolio preferences')
  }

  w <- sample.weights.with.pp(proj.perfs, res.prev$nd, ok.popts, w.constr)

  acc <- comp.acc(res.prev$nd, w, proj.perfs)
  res <- estimate.ci(res.prev$nd, w, proj.perfs)

  list(acc=acc, ci=res$ci, popt.pfs=res$popt.pfs)
}
