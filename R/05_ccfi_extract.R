P2<-Sys.getenv("P2"); if(P2=="")P2<-"/Users/menghao/Downloads/study1_repro/phase2"
band<-function(c) if(is.na(c))"NA" else if(c<.45)"dimensional" else if(c>.55)"TAXONIC" else "ambiguous"
traits<-c("NEO_Neuroticism","NEO_Extraversion","NEO_Openness","NEO_Agreeableness","NEO_Conscientiousness",
 "HEX_HonestyHumility","HEX_Emotionality","HEX_eXtraversion","HEX_Agreeableness","HEX_Conscientiousness","HEX_Openness")
A<-setNames(rep(NA_real_,length(traits)),traits)
for(t in traits){ f<-file.path(P2,"out",paste0("ccfi_",t,".rds"))
 if(file.exists(f)){ r<-readRDS(f); A[t]<-r$CCFI.mean
  cat(sprintf("%-24s CCFI=%.3f  (MAMBAC %.2f  MAXEIG %.2f  LMode %.2f)  -> %s\n",
   t,r$CCFI.mean,r$CCFI.MAMBAC,r$CCFI.MAXEIG,r$CCFI.LMode,band(r$CCFI.mean)))
 } else cat(sprintf("%-24s (not done)\n",t)) }
saveRDS(A,file.path(P2,"out","ccfi_summary.rds"))
cat("\n===== H2 verdict (standard band: <.45 dimensional, .45-.55 ambiguous, >.55 taxonic) =====\n")
cat("NEO mean CCFI:   ",round(mean(A[grep('^NEO',names(A))],na.rm=T),3),"\n")
cat("HEXACO mean CCFI:",round(mean(A[grep('^HEX',names(A))],na.rm=T),3),"\n")
cat("all traits < .45?",all(A<.45,na.rm=T)," | any trait >= .55 (taxonic)?",any(A>=.55,na.rm=T),"\n")
