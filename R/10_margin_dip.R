# 10_margin_dip.R — Hartigan dip tests on the trait margins of every analysis subsample.
# Raw summed-Likert scores are granular (ties), which the dip test reads as departures from
# unimodality at these sample sizes; a small tie-breaking jitter (uniform within half the minimum
# score gap) removes that artifact. Reported in the Results power section: after tie-breaking,
# 64 of the 65 margins retain unimodality. The copula null inherits each margin exactly, so
# marginal granularity cannot separate real data from null in any case.
suppressMessages(library(diptest))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
tot<-0; rej<-0
for(nm in names(slices)){
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); x<-d[sample(nrow(d),min(slices[[nm]],nrow(d))),,drop=FALSE]
  gap<-min(diff(sort(unique(x[,1]))))
  set.seed(20260620); xj<-x+matrix(runif(length(x),-gap/2,gap/2),nrow(x))
  p_raw<-apply(x,2,function(v) dip.test(v)$p.value)
  p_jit<-apply(xj,2,function(v) dip.test(v)$p.value)
  tot<-tot+ncol(x); rej<-rej+sum(p_jit<.05)
  cat(sprintf("%-14s margins=%2d | raw p<.05: %2d | jittered p<.05: %2d | max dip D=%.4f\n",
    nm,ncol(x),sum(p_raw<.05),sum(p_jit<.05),max(apply(x,2,function(v) unname(dip.test(v)$statistic)))))
}
cat(sprintf("TOTAL: %d of %d margins reject unimodality after tie-breaking\n",rej,tot))
