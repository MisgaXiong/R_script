#!/usr/bin/Rscript

#
# NAME
#
#   Build_lysinDB.R
#
# DESCRIPTION1
#
# AUTHOR
#
#   Dongyan Xiong
#
# VERSION
#
#   0.0.1   2019-06-17
#	
#
# USAGE
#	Rscript Build_lysinDB.R gbk_file_PATH gbk_file_name out_lysin_nt_path out_lysin_aa_path
#
# DESCRIPTION2
#	Batch extracte phage lysin


arg<-commandArgs(T)
setwd(arg[1])
library(genbankr)
library(Biostrings)
phage_gb<-readGenBank(arg[2])
phage_lysin<-c("endolysin","amidase","hydrolase","lysozyme","lysin","glucosidase","endopeptidase","CHAP",
               "N-acetylmuramidase","muramidase","transglycosylase","endo-beta-N-acetylglucosaminidase",
               "autolysin","LysM","N-acetylglucosaminidase","endo-β-N-acetylglucosaminidase","peptidase","Lysozyme")
lysin_NT<-DNAStringSet()
lysin_AA<-AAStringSet()
lysin_nt<-DNAStringSet()
lysin_aa<-AAStringSet()
for (i in 1:length(phage_gb@cds)) {
  cds_product<-phage_gb@cds$product[i]
  cds_product<-strsplit(cds_product," ")[[1]]
  if(length(intersect(cds_product,phage_lysin))>0){
    lysin_aa<-AAStringSet(phage_gb@cds$translation[i])
    names(lysin_aa)<-paste(names(phage_gb@sequence),phage_gb@cds$product[i],sep = " ")
    l<-as.data.frame(phage_gb@cds@ranges[i])[1,1]
    r<-as.data.frame(phage_gb@cds@ranges[i])[1,2]
    lysin_nt<-DNAStringSet(phage_gb@sequence[[1]][l:r])
    names(lysin_nt)<-paste(names(phage_gb@sequence),phage_gb@cds$product[i],sep = " ")
    if(toString(lysin_nt[[1]][1:3])!="ATG"){
      lysin_nt<-DNAStringSet(reverseComplement(lysin_nt))
    }else{
      lysin_nt<-lysin_nt
    }
  }else{
    next
  }
  lysin_NT<-c(lysin_NT,DNAStringSet(lysin_nt))
  lysin_AA<-c(lysin_AA,AAStringSet(lysin_aa))
}
out1<-tempfile(tmpdir = arg[3])
writeXStringSet(lysin_NT,out1)
out2<-tempfile(tmpdir = arg[4])
writeXStringSet(lysin_AA,out2)
