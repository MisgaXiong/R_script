#=================================================
#Find_diffgene2genetype.R
#Usage: Find_diffgene2China.R path gene_type.tsv
#Athor: Dongyan Xiong
#=================================================
#!/usr/bin/Rscript
arg<-commandArgs(T)
setwd(arg[1])
CHA_gene_tbl<-read.table(arg[2])
rownames(CHA_gene_tbl)<-as.character(CHA_gene_tbl[,1])
CHA_gene_tbl<-CHA_gene_tbl[,2:ncol(CHA_gene_tbl)]
diff_gene<-NULL
for (i in 1:ncol(CHA_gene_tbl)) {
  if(sd(as.numeric(CHA_gene_tbl[2:nrow(CHA_gene_tbl),i]))!=0){
    diff_gene<-cbind(diff_gene,as.vector(CHA_gene_tbl[,i]))
  }
}
rownames(diff_gene)<-rownames(CHA_gene_tbl)
diff_gene<-data.frame(diff_gene)
colnames(diff_gene)<-NULL
write.table(diff_gene,"diffgene_CHA.tsv",sep = "\t")
