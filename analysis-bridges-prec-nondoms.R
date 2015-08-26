options(java.parameters = "-Xmx4g")
library(xlsx)
library(rpm)
library(hitandrun)
library(plyr)
library(xtable)
library(ror)
library(ggplot2)
library(gridExtra)
source('weights.R')
source('funcs.R')

## DPS = VPS in the original publication (Finnish abbreviation for Damage Point Sum)
cnames <- c('DPS', 'Traffic', 'Carry', 'Width', 'Salt', 'Visual', 'Cost', 'DPS Red')
fname <- 'siltadata-precomp.xlsx'
data <- read.xlsx(fname, sheetName='projects')
projects <- data[,3:10]
colnames(projects) <- cnames
n.samples <- 1E4
n <- 6 # 6 criteria
set.seed(1911)
alt.perm <- sample(1:nrow(projects), nrow(projects), replace=FALSE)
projects <- projects[alt.perm,] # random permutation of rows because original data is sorted
rownames(projects) <- 1:nrow(projects)
proj.perfs <- projects[,1:n]

make.analysis <- function(type) { # 'prefs' or 'noprefs'
  nd <- read.xlsx2(fname, sheetName=paste0('nd-', type), header=TRUE)[,-1]
  nd <- t(aaply(as.matrix(nd), 2, as.numeric))
  nd <- nd[,alt.perm]

  w.constr <- if (type == 'prefs') silta.constraints() else simplexConstraints(n)
  
  verts <- make.vertices(n, w.constr)
  
  ## Compute non-doms
  nd.perfs <- nd %*% as.matrix(proj.perfs)
  cind.nd <- colSums(nd) / nrow(nd)
  
  ## Compute pot. opts
  pot.opts <- potopt.indices(nd.perfs, w.constr)
  po.projs <- which(colSums(nd[pot.opts,]) > 0)
  cind.popt <- colSums(nd[pot.opts,]) / length(pot.opts)
  
  ## Compute weight shares
  samples <- hitandrun(w.constr, n.samples=n.samples)

  w.ok.opts <- comp.acc(nd, samples, proj.perfs)
  
  proj.info <- cbind(projects, w.ok.opts, cind.popt, cind.nd)
  
  list(proj=proj.info, nd=nd, popt.pfs=pot.opts)
}

out.projs <- c(17, 1, 3, 6, 8, 18, 20, 28)

res.noprefs <- make.analysis('noprefs')
res.prefs <- make.analysis('prefs')
res.xin1 <- analysis.pp(proj.perfs, res.prefs, 45, NULL, silta.constraints())
res.xout1 <- analysis.pp(proj.perfs, res.prefs, 45, out.projs, silta.constraints())

filter.res <- function(res) {
  popt.projs <- res$nd[res$popt.pfs,]
  cs <- colSums(popt.projs)
  exterior <- as.logical(cs == 0)
  core <- as.logical(cs == nrow(popt.projs))
  to.ret <- res$proj[,(ncol(res$proj)-2):(ncol(res$proj)-1)]
  to.ret
}

res.all <- cbind(projects[,1:n], filter.res(res.noprefs), filter.res(res.prefs))
colnames(res.all)[(ncol(res.all)-3):ncol(res.all)] <- c('AI.no', 'CI.no', 'AI.ord', 'CI.ord')

res.all <- cbind(res.all, res.xin1$acc, res.xin1$ci)
colnames(res.all)[(ncol(res.all)-1):ncol(res.all)] <- c('AI.xi1', 'CI.xi1')

res.all <- cbind(res.all, res.xout1$acc, res.xout1$ci)
colnames(res.all)[(ncol(res.all)-1):ncol(res.all)] <- c('AI.xo1', 'CI.xo1')

res.all <- cbind(res.all, t(res.prefs$nd[res.xout1$popt.pfs,]))

## Order re: acceptabilities
res.all <- res.all[order(res.all[,7], decreasing=TRUE),]
res.cols <- c(cnames[1:n], rep(c('AI', 'CI'), 4), '$p^1$', '$p^2$', '$p^3$')
colnames(res.all) <- res.cols
acc.ci <- res.all[,-(1:6)]

acc.ci.flt <- do.call(cbind, llply(1:4, function(idx) {
  rng <- (((idx-1)*2)+1):(idx*2)
  d <- acc.ci[,rng]
  nas <- acc.ci[,idx*2] == 0
  d[nas,] <- c(NA, NA)
  d
}))
res.all.flt <- cbind(res.all[,1:6], acc.ci.flt, res.all[,(ncol(res.all)-2):ncol(res.all)])
colnames(res.all.flt) <- res.cols

to.plot <- which(rowSums(acc.ci) != 0)
acc.to.plot <- acc.ci[to.plot,c(1,3,5,7)]
ci.to.plot <- acc.ci[to.plot,c(2,4,6,8)]

df <- as.data.frame(cbind(rownames(acc.to.plot), acc.to.plot, ci.to.plot))
colnames(df) <- c('x', 'a1', 'a2', 'a3', 'a4', 'c1', 'c2', 'c3', 'c4')
df$x2 <- factor(df$x, rev(as.character(df$x)))

theme_set(theme_grey(base_size = 8))

pdf('caseres.pdf')
pa1 <- ggplot(data=df, aes(x=x2, y=a1, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(AI(x^k, W^0))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pa2 <- ggplot(data=df, aes(x=x2, y=a2, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(AI(x^k, W^1))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pa3 <- ggplot(data=df, aes(x=x2, y=a3, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(AI(x^k, W^2))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pa4 <- ggplot(data=df, aes(x=x2, y=a4, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(AI(x^k, W^3))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pc1 <- ggplot(data=df, aes(x=x2, y=c1, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(CI(x^k, W^0))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pc2 <- ggplot(data=df, aes(x=x2, y=c2, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(CI(x^k, W^1))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pc3 <- ggplot(data=df, aes(x=x2, y=c3, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(CI(x^k, W^2))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
pc4 <- ggplot(data=df, aes(x=x2, y=c4, fill=NULL)) + geom_bar(stat='identity') + coord_flip() + ylab(bquote(CI(x^k, W^3))) + xlab(bquote(x^k)) + ylim(0, 1) + theme(axis.text.y=element_text(size=4, face='plain'))
grid.arrange(pa1, pc1, pa2, pc2, pa3, pc3, pa4, pc4, ncol=2)
dev.off()

## Format
data.digits <- c(0, 3, 1, 0, 0, 2, 2, 0, 0)
res.digits <- c(0, 3, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0)
data.align <- c('rcccccccc')
res.align <- c('rcccccc|cc|cc|cc|cc|ccc')

## Write to latex
sink('results.tex')
print(xtable(res.all.flt, label='results', digits=res.digits, align=res.align,
             caption=paste0("Project criteria measurements, acceptability (AI) and core indices (CI) for core- and borderline projects in 4 iterations of the analysis: without weight information, with linear weight constraints, and with addition of project preference statements In($x^{45}$) and Out($",
               paste0('x^{', out.projs, '}', collapse=', '),
               "$). The projects are sorted according to AI in the preference-free analysis")),
      latex.environments=c('tiny', 'center'),
      sanitize.text.function=identity,
      caption.placement='top')
sink(NULL)

sink('bridges.tex')
print(xtable(projects, digits=data.digits, align=data.align,
             caption="Bridge criteria measurements, costs, and DPS reductions"),
      latex.environments=c('tiny', 'center'),
      caption.placement='top')
sink(NULL)
