suppressPackageStartupMessages(library(RTaxometrics))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"/Users/menghao/Downloads/study1_repro/phase2"
neo<-as.matrix(read.csv(file.path(P2,"data","conf_NEO120_facet.csv"),header=FALSE))   # 30 facets
hex<-as.matrix(read.csv(file.path(P2,"data","conf_HEXACO_facet.csv"),header=FALSE))    # 24 facets
neo_b<-list(NEO_Neuroticism=1:6,NEO_Extraversion=7:12,NEO_Openness=13:18,NEO_Agreeableness=19:24,NEO_Conscientiousness=25:30)
hex_b<-list(HEX_HonestyHumility=1:4,HEX_Emotionality=5:8,HEX_eXtraversion=9:12,HEX_Agreeableness=13:16,HEX_Conscientiousness=17:20,HEX_Openness=21:24)
band<-function(c) if(is.na(c))"NA" else if(c<.45)"dimensional" else if(c>.55)"TAXONIC" else "ambiguous"
runtr<-function(data,cols,label){
  x<-data[,cols]; set.seed(20260620); x<-x[sample(nrow(x),min(1300,nrow(x))),]
  res<-tryCatch(RunCCFIProfile(as.data.frame(x),num.p=11,n.pop=5000,n.samples=18,graph=0),
                error=function(e){cat("ERR",label,":",conditionMessage(e),"\n");NULL})
  saveRDS(res,file.path(P2,"out",paste0("ccfi_",label,".rds")))
  ccfi<-if(is.null(res)) NA_real_ else {v<-suppressWarnings(as.numeric(unlist(res))); mean(v[v>0&v<1],na.rm=TRUE)}
  cat(sprintf("%-26s CCFI=%.3f -> %s\n",label,ccfi,band(ccfi))); flush.console(); ccfi
}
A<-c()
for(nm in names(neo_b)) A[nm]<-runtr(neo,neo_b[[nm]],nm)
for(nm in names(hex_b)) A[nm]<-runtr(hex,hex_b[[nm]],nm)
saveRDS(A,file.path(P2,"out","ccfi_all.rds"))
cat("\n===== H2 verdict (standard band: <.45 dimensional, >.55 taxonic) =====\n")
cat("NEO mean CCFI:   ",round(mean(A[grep('^NEO',names(A))],na.rm=T),3),"\n")
cat("HEXACO mean CCFI:",round(mean(A[grep('^HEX',names(A))],na.rm=T),3),"\n")
cat("any trait >= .55 (taxonic)?",any(A>=.55,na.rm=T)," | all < .45?",all(A<.45,na.rm=T),"\n")
cat("DONE\n")
