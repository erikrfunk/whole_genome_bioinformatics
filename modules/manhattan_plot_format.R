setwd("/data2/rosyfinches/Fst_files/")
FST<-read.table("temp.fst",header=T)
# FILE<-read.table("GC_vs_BL.weir.fst",header=T)
cat("New file is", length(FST[,1]),"snps long.\n")
snp_column<-c(1:length(FST[,1]))
FST<-cbind(SNP=as.numeric(snp_column),FST)
cat("Added SNP column.")

write.table(FST,"temp2.fst")
