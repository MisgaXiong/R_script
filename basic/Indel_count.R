#!/usr/bin/Rscript

#
# NAME
#
#   Indel_count2plot
#
# DESCRIPTION1
#
# AUTHOR
#
#   Dongyan Xiong
#
# VERSION
#
#   0.0.1   2019-05-31
#	0.0.2	2019-06-01
#
# USAGE
#	Rscript Indel_count.R SNP_file_PATH Ref_gbk_file_path Ref_gbk_file output_name
#
# DESCRIPTION2
#	SNP_file is the output by MUMmer


#function_1: Insert count
Indel_insert_count<-function(x){
  Indel_data<-x
  Indel_position<-NULL
  Indel_insert_info<-data.frame()
  Indel_insert_START<-NULL
  Indel_insert_LEN<-NULL
  for (i in 1:(nrow(Indel_data)-1)) {
    TRF <- Indel_data[i,1] %in% Indel_position 
    if(as.character(Indel_data[i,2])=="." && TRF=="FALSE" ){
      Indel_insert_start=Indel_data[i,1]
      Indel_insert_end<-Indel_data[i,1]
      Indel_insert_len=1
      for (j in (i+1):nrow(Indel_data)) {
        if(as.character(Indel_data[j,2])=="." && Indel_data[j,5]==0){
          Indel_insert_len=Indel_insert_len+1
          Indel_insert_end<-c(Indel_insert_end,Indel_data[j,1])
        }else{break}
      }
      Indel_position<-c(Indel_position,Indel_insert_end)
      Indel_insert_START<-c(Indel_insert_START,Indel_insert_start)
      Indel_insert_LEN<-c(Indel_insert_LEN,Indel_insert_len)
    }else{next}
  }
  Indel_insert_info<-data.frame(position=Indel_insert_START,length=Indel_insert_LEN)
}

#function_2: Deletion count
Indel_delete_count<-function(x){
  Indel_data<-x
  Indel_position<-NULL
  Indel_delete_info<-data.frame()
  Indel_delete_START<-NULL
  Indel_delete_LEN<-NULL
  for (i in 1:(nrow(Indel_data)-1)) {
    TRF <- Indel_data[i,1] %in% Indel_position 
    if(as.character(Indel_data[i,3])=="." && TRF=="FALSE" ){
      Indel_delete_start=Indel_data[i,1]
      Indel_delete_end<-Indel_data[i,1]
      Indel_delete_len=1
      for (j in (i+1):nrow(Indel_data)) {
        if(as.character(Indel_data[j,3])=="." && Indel_data[j,5]==1){
          Indel_delete_len=Indel_delete_len+1
          Indel_delete_end<-c(Indel_delete_end,Indel_data[j,1])
        }else{break}
      }
      Indel_position<-c(Indel_position,Indel_delete_end)
      Indel_delete_START<-c(Indel_delete_START,Indel_delete_start)
      Indel_delete_LEN<-c(Indel_delete_LEN,Indel_delete_len)
    }else{next}
  }
  Indel_delete_info<-data.frame(position=Indel_delete_START,length=Indel_delete_LEN)
}

#Main Program
library(stringr)
arg<-commandArgs(T)
setwd(arg[1])
dir1<-dir(arg[1])
Insert_data<-data.frame()
Delet_data<-data.frame()
for (i in 1:length(dir1)) {
  Insert_data<-rbind(Insert_data,Indel_insert_count(read.table(dir1[i],sep = "\t")))
  Delet_data<-rbind(Delet_data,Indel_delete_count(read.table(dir1[i],sep = "\t")))
}
Indel_all<-rbind(Insert_data,Delet_data)

#Coding && Non-coding region Indel count
library(genbankr)
setwd(arg[2])
gene_region<-readGenBank(arg[3])
Indel_region<-NULL
gene_region<-as.data.frame(gene_region@genes@ranges)
gene_range<-NULL
for (i in 1:nrow(gene_region)) {
  gene_range<-c(gene_range,c(gene_region[i,1]:gene_region[i,2]))
}
for (i in 1:nrow(Indel_all)) {
  if(Indel_all[i,1] %in% gene_range){
    Indel_region<-c(Indel_region,"Coding_region")
  }else{Indel_region<-c(Indel_region,"Non-coding_region")}
}
Indel_region<-data.frame(type=Indel_region)
All_Indel_data<-cbind(Indel_all,Indel_region)

#Point Plot
library(ggplot2)
fout<-arg[4]
p<-ggplot(data = All_Indel_data)+geom_point(alpha=0.2,mapping = aes(x=position,y=length, color=type))
p<-p+theme_classic() +
  theme(axis.ticks.y = element_blank(), )
ggsave(fout, p, device = "svg")
