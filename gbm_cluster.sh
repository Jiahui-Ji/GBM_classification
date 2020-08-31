

#load packages
library(ConsensusClusterPlus)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(NMF)
library(CancerSubtypes)
library(mixtools)
library(splus2R)
library(RTCGA)


#download tcga dataset
query=GDCquery(project = "TCGA-GBM",
					  legacy = FALSE, 
					  data.category="Transcriptome Profiling",
					  data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")
GDCdownload(query,directory="./project")
            
data=GDCprepare(query,directory = "./project")
count_data=assay(data)
rna.data=as.matrix(count_data)






#download tcga fpkm dataset
query_fpkm=GDCquery(project = "TCGA-GBM",
					  legacy = FALSE, 
					  data.category="Transcriptome Profiling",
					  data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM")
GDCdownload(query_fpkm,directory="./project2")
            
data_fpkm=GDCprepare(query_fpkm,directory = "./project2")
count_data_fpkm=assay(data_fpkm)
rna.data.fpkm=as.matrix(count_data_fpkm)



#download clinical information
query=GDCquery(project = "TCGA-GBM",data.category="Clinical",
               data.type="Clinical Supplement", data.format="BCR Biotab")
GDCdownload(query)
data=GDCprepare(query)					  





library(ConsensusClusterPlus)
library(dplyr)
library(SummarizedExperiment)
library(NMF)
library(CancerSubtypes)
library(mixtools)
library(splus2R)
setwd('/rds/general/user/jj1419/ephemeral/gbm')
load('11_05_2020_gbmtcga.RData')

rna.data

rna.data.fpkm




#library size (sum), mean, SD and filter the data
rna.library.size=matrix(nrow=length(rownames(rna.data)),ncol=3)
rownames(rna.library.size)=rownames(rna.data)
colnames(rna.library.size)=c('library_size','mean','SD')
for (i in 1:length(rownames(rna.data)))
{
	rna.library.size[i,1]=sum(rna.data[i,])
	rna.library.size[i,2]=mean(rna.data[i,])
	rna.library.size[i,3]=sd(rna.data[i,])
	
}	

#rna.library.size

q=which(rna.library.size[,3]!=0)
length(q)
dim(rna.data)
rna.data.filter=rna.data[q,]
dim(rna.data.filter)

rna.library.size.filter=rna.library.size[q,] 

png("rna_expression.png")
hist(log2(rna.library.size.filter[,1])) #sum, library size
hist(log2(rna.library.size.filter[,2])) #mean
hist(log2(rna.library.size.filter[,3])) #sd
dev.off()
####### to this step, we acquire filtered gene matrix






#fpkm to tpm
fpkmToTpm=function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
rna.data.tpm=apply(rna.data.fpkm,2,fpkmToTpm)



gene.name.filter=intersect(rownames(rna.data)[q],rownames(rna.data.tpm)) 
rna.data.tpm.filter=rna.data.tpm[gene.name.filter,]

#we get rna.data.tpm.filter 53453x174








#keep samples that have at least 10 counts in at least 10 samples
keep=rowSums(rna.data.tpm.filter>3)>=2 
summary(keep)
tpm_data_filtered=rna.data.tpm.filter[keep,]
#we get tpm_data_filtered 26100x174



tpm.library.size=matrix(nrow=length(colnames(rna.data.tpm.filter)),ncol=3)
rownames(tpm.library.size)=colnames(rna.data.tpm.filter)
colnames(tpm.library.size)=c('library_size','mean','SD')
for (i in 1:length(colnames(rna.data.tpm.filter)))
{
	tpm.library.size[i,1]=sum(rna.data.tpm.filter[,i])
	tpm.library.size[i,2]=mean(rna.data.tpm.filter[,i])
	tpm.library.size[i,3]=sd(rna.data.tpm.filter[,i])
	
}

tpm.library.size

hist(log2(tpm.library.size[,1])) #sum, library size
hist(log2(tpm.library.size[,2])) #mean
hist(log2(tpm.library.size[,3])) #sd

pdf('tpm_boxplot2.pdf')
boxplot(log2(rna.data.tpm.filter[,1:50]))
dev.off()





#log2(count number) log2 transfered filter matrix 
tpm.data.lg2=log2(tpm_data_filtered+0.1)





#EM normalization use tpm matrix
exprs=log2(tpm_data_filtered+0.1)
M=apply(exprs,1,mean,na.rm=T)

set.seed(101)
yy=normalmixEM(M, mu=quantile(M,c(0.25,0.75), lambda=c(0.5,0.5)))

Mth=qnorm(0.95,min(yy$mu),yy$sigma[which.min(yy$mu)])
EM=exprs[M>Mth,]
EM[is.na(EM)] <- 0
#Head(EM) #  data.frame dimensions :	  matrix dimensions :     12183 8 
pdf('tpm_em_filter.pdf')
plot(density(M, na.rm=T))
plot(yy, which=2) # inbuilt
lines(density(M, na.rm=T), lty=2, lwd=2)
hist(EM)
dev.off()








#pca analysis 
EM

dataFit=TCGAanalyze_Filtering(tabDF=EM, method="quantile", qnt.cut=0.25)
pca_plot=TCGAvisualize_PCA(EM, 

pdf('pca_gbm.pdf')
plot(pca_plot)
dev.off()





















#change count 0 to 0.1
for (i in 1:length(colnames(rna.data.tpm)))
{
	for (j in 1:length(rownames(rna.data.tpm)))
	{
		rna.data.tpm[j,i]
	
}
}



















#normalization

#calculate the average expression of each gene to generate gene list
gene.list=matrix(NA,ncol=1,nrow=length(rownames(rna.data.tpm_lg2)))
rownames(gene.list)=rownames(rna.data.tpm_lg2)

for (i in 1:length(rownames(rna.data.tpm_lg2)))
{
	gene.list[i,1]=sum(rna.data.tpm_lg2[i,])/length(colnames(rna.data.tpm_lg2))

}



#normalmixEM
set.seed(100)
out=normalmixEM(gene.list, arbvar = FALSE, epsilon = 1e-03)
str(out)

pdf(file="gene_his.pdf")
hist(gene.list)
dev.off()

pdf(file="gene_em.pdf")
plot(out,whichplots=2)
dev.off()

















#consenesus clustering
rcc3 = ConsensusClusterPlus(rna.data.tpm,maxK=10,reps=174,pItem=0.8,pFeature=1,
title="gbm_rna_hc",distance="euclidean",clusterAlg="hc",plot="png")

rna.clusters=rcc3[[3]]$consensusClass

cluster1=c()
cluster2=c()
cluster3=c()
for (i in 1:length(rna.clusters))
{
	if (rna.clusters[[i]]==1)
	{cluster1=append(cluster1,i)}
	
	if (rna.clusters[[i]]==2)
	{cluster2=append(cluster2,i)}
	
	if (rna.clusters[[i]]==3)
	{cluster3=append(cluster3,i)}
	
}

cluster1.rna=rna.data[,cluster1]
cluster2.rna=rna.data[,cluster2]
cluster3.rna=rna.data[,cluster3]










#plot trees
hcd=hclust(dist(rna.data))
plot(hcd)
plot(as.phylo(hcd), cex = 0.2, type = "unrooted")
plot(as.phylo(hcd), cex = 0.2, type = "fan")


pdf("gbm_tree1.pdf")
plot(hcd)
dev.off()

pdf("gbm_tree2.pdf")
plot(as.phylo(hcd), cex = 0.2, type = "unrooted")
dev.off()

pdf("gbm_tree3.pdf")
plot(as.phylo(hcd), cex = 0.2, type = "fan")
dev.off()













#consensus non-negative matrix factorization
result=ExecuteCNMF(rna.data.tpm,clusterNum=3,nrun=30)


ncluster1=c()
ncluster2=c()
ncluster3=c()
for (i in 1:length(rna.clusters))
{
	if (result[[1]][i]==1)
	{ncluster1=append(ncluster1,i)}
	
	if (result[[1]][i]==2)
	{ncluster2=append(ncluster2,i)}
	
	if (result[[1]][i]==3)
	{ncluster3=append(ncluster3,i)}
	
}

ncluster1.rna=rna.data[,ncluster1]
ncluster2.rna=rna.data[,ncluster2]
ncluster3.rna=rna.data[,ncluster3]






###PCA 
library('ggpubr')

#plot pc1 vs pc2
pca.info=prcomp(t(EM))

pdf('gbm_pca.pdf')
autoplot(pca.info)
dev.off()

#paris plots 
pdf('gbm_pca_pair.pdf')
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca.info$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(pca.info$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()




#clinical data
clinical_gbm=GDCquery_clinic('TCGA-GBM','clinical')

clinical_gbm_matrix=as.matrix(clinical_gbm)
clinical_id=clinical_gbm_matrix[,1]


#deal with the names of clinical samples 
num_clinic=c()
for (i in 1:length(clinical_id))
{
	q=grep(clinical_id[i],colnames(EM),value=F)
	if (identical(q,integer(0)))		
	{}
	else
	{num_clinic=append(num_clinic,i)}
	
}


clinical_gbm_matrix_filter=clinical_gbm_matrix[num_clinic,]



#change submitter id into full id 
clinical_id_v2=clinical_gbm_matrix_filter[,1]


name=c()
for (i in 1:length(clinical_id_v2))
{
	
	q=grep(clinical_id_v2[i],colnames(EM),value=T)
	name=append(name,q[1])
}

dim(clinical_gbm_matrix_filter) #166*158


clinical_gbm_matrix_filter=cbind(clinical_gbm_matrix_filter,name)





#based on the clinical data to seperate the data into different groups

#pca on gender
EM.filter=EM[,name]
dim(EM.filter) #14717*166

type=clinical_gbm_matrix_filter[,"gender"]
type=factor(type,levels=c("male","female","NA"))
colType=c("CornflowerBlue", "ForestGreen","DarkOrange")[type]
pchType <- c(16)
pca=prcomp(t(EM.filter))

pdf("gbm_gender_pca.pdf")
plot(pca$x,
  col = colType,
  pch = pchType,
  cex = 1.0)
legend(
  "bottomright",
  bty = "n",
  c("male","female"),
  fill = c("CornflowerBlue", "ForestGreen"),
  cex = 2.0)
dev.off()



pdf("gbm_gender_pca_v2.pdf")

plot(pca$x, type="n", main="Principal components analysis bi-plot", xlab="PC1", ylab="PC2")
points(pca$x, col=colType, pch=16, cex=1)

plot(pca$x[,1], pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab="PC1", ylab="PC3")
points(pca$x[,1], pca$x[,3], col=colType, pch=16, cex=1)

plot(pca$x[,2], pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab="PC2", ylab="PC3")
points(pca$x[,2], pca$x[,3], col=colType, pch=16, cex=1)

dev.off()

#pca on tumor type
type=tumor.type
type=factor(type,levels=c('solid01','recurrent02','normal11'))
colType=c("CornflowerBlue", "ForestGreen","DarkOrange")[type]
pchType <- c(16)

pdf("gbm_tumortype_pca.pdf")
plot(pca$x,
  col = colType,
  pch = pchType,
  cex = 1.0)
legend(
  "bottomright",
  bty = "n",
  c('solid01','recurrent02','normal11'),
  fill = c("CornflowerBlue", "ForestGreen","DarkOrange"),
  cex = 1.0)
dev.off()


















#seperate data
sample.01=c(grep("-01A",name,value=T),grep("-01B",name,value=T),grep("-01C",name,value=T))
sample.02=grep("-02A",name,value=T)
sample.11=grep("-11A",name,value=T)

#pca1-3 anova test 
tumor.type=c()
for (i in 1:length(rownames(pca.1_3)))
{
	if (rownames(pca.1_3)[i] %in% sample.01)
	{tumor.type=append(tumor.type,"solid01")}
	if (rownames(pca.1_3)[i] %in% sample.02)
	{tumor.type=append(tumor.type,"recurrent02")}
	if (rownames(pca.1_3)[i] %in% sample.11)
	{tumor.type=append(tumor.type,"normal11")}

}

pca.1_3.ano=cbind(as.numeric(pca$x[,1]),as.numeric(pca$x[,2]),as.numeric(pca$x[,3]),tumor.type)
colnames(pca.1_3.ano)=c('PC1','PC2','PC3','tumor.type')
pca.1_3.ano=as.data.frame(pca.1_3.ano)

pca.1_3.ano$PC1=as.numeric(as.character(pca.1_3.ano$PC1))
pca.1_3.ano$PC2=as.numeric(as.character(pca.1_3.ano$PC2))
pca.1_3.ano$PC3=as.numeric(as.character(pca.1_3.ano$PC3))

#anova test 
aov1=aov(PC1~tumor.type,data=pca.1_3.ano)
summary(aov1)
aov2=aov(PC2~tumor.type,data=pca.1_3.ano)
summary(aov2)
aov3=aov(PC3~tumor.type,data=pca.1_3.ano)
summary(aov3)







#fit linear model 
head(pca.1_3.ano)
dim(clinical_gbm_matrix_filter)

lm1=lm(PC1~tumor.type, data=pca.1_3.ano)
summary(lm1)
lm2=lm(PC2~tumor.type, data=pca.1_3.ano)
summary(lm2)
lm3=lm(PC3~tumor.type, data=pca.1_3.ano)
summary(lm3)






pca.lm=pca.1_3.ano


clinic=as.data.frame(clinical_gbm_matrix_filter)
col.name=c()
for (i in 1:length(names(clinic)))
{
	q=length(levels(clinic[,i]))
	if ( q>1 & q<160 )
	{pca.lm=cbind(pca.lm,clinic[,i])
	col.name=append(col.name,i)}	
}

col.name=names(clinic)[col.name]
col.name.n=c('PC1','PC2','PC3','tumor.type',col.name)


names(pca.lm)=col.name.n
write.table(pca.lm,'gbm_clinic.csv',quote=F,sep=',')






p.matrix=matrix(nrow=18,ncol=3)
for (i in 1:length(names(pca.lm))-3)
{
	lm.1=lm(pca.lm[,1]~pca.lm[,i+3])
	lm.2=lm(pca.lm[,2]~pca.lm[,i+3])
	lm.3=lm(pca.lm[,3]~pca.lm[,i+3])
	
	f1=summary(lm.1)$fstatistic
	p_value.1= pf(as.numeric(f1[1]), as.numeric(f1[2]), as.numeric(f1[3]), lower.tail = FALSE)
	
	f2=summary(lm.2)$fstatistic
	p_value.2= pf(as.numeric(f2[1]), as.numeric(f2[2]), as.numeric(f2[3]), lower.tail = FALSE)
	
	f3=summary(lm.3)$fstatistic
	p_value.3= pf(as.numeric(f3[1]), as.numeric(f3[2]), as.numeric(f3[3]), lower.tail = FALSE)
	
	
	p.matrix[i,1]=p_value.1
	p.matrix[i,2]=p_value.2
	p.matrix[i,3]=p_value.3
}	
colnames(p.matrix)=c('PC1','PC2','PC3')
rownames(p.matrix)=c('tumor_type',col.name)
write.table(p.matrix,'gbm_variates_lm.csv',quote=F,sep=',')




#####test 

p.matrix=matrix(nrow=10,ncol=3)
for (i in 114:124)
{
	lm.1=lm(pca.lm[,1]~pca.lm[,i+3])
	f1=summary(lm.1)$fstatistic
	p_value.1= pf(as.numeric(f1[1]), as.numeric(f1[2]), as.numeric(f1[3]), lower.tail = FALSE)
	print(f1)
	
	lm.2=lm(pca.lm[,2]~pca.lm[,i+3])
	f2=summary(lm.2)$fstatistic
	p_value.2= pf(as.numeric(f2[1]), as.numeric(f2[2]), as.numeric(f2[3]), lower.tail = FALSE)
	print(p_value.2)
	
	lm.3=lm(pca.lm[,3]~pca.lm[,i+3])
	f3=summary(lm.3)$fstatistic
	p_value.3= pf(as.numeric(f3[1]), as.numeric(f3[2]), as.numeric(f3[3]), lower.tail = FALSE)
	print(p_value.3)
	
	#p.matrix[i,1]=p_value.1
	#p.matrix[i,2]=p_value.2
	#p.matrix[i,3]=p_value.3
}	



sample.01=c(grep("-01A",name,value=T),grep("-01B",name,value=T),grep("-01C",name,value=T))
sample.02=grep("-02A",name,value=T)
sample.11=grep("-11A",name,value=T)



t=grep("-11A",name,value=F)
test=append(c(1:160),c(164))

test=


pca.lm=pca.lm[test,]
write.table(p.matrix,'gbm_pvalue_tumor.csv',quote=F,sep=',')







# Biospecimen - technical differences
clin=GDCquery_clinic("TCGA-GBM", type = "biospecimen", save.csv = TRUE)


clin=as.matrix(clin) 

sub.name=clinical_gbm_matrix_filter[,1]
name






clin.name=c()
for (i in 1:length(clin[,25]))
{
	a=grep(clin[,25][i],name)
	
	if (identical(a,integer(0)))		
	{}
	else
	{clin.name=append(clin.name,i)}
}


clin.filter=clin[clin.name,]




biospecimen.name=c()
for (i in sub.name)
{
	biospecimen.name=append(biospecimen.name,grep(i,clin.filter[,25]))

}

clin.filter.2=clin.filter[biospecimen.name,]






pca.lm=pca.1_3.ano


clinic=as.data.frame(clin.filter.2)
col.name=c()
for (i in 1:length(names(clinic)))
{
	q=length(levels(clinic[,i]))
	if ( q>1 & q<160 )
	{pca.lm=cbind(pca.lm,clinic[,i])
	col.name=append(col.name,i)}	
}

col.name=names(clinic)[col.name]
col.name.n=c('PC1','PC2','PC3','tumor.type',col.name)


names(pca.lm)=col.name.n
write.table(pca.lm,'gbm_biospecimen.csv',quote=F,sep=',')






p.matrix=matrix(nrow=7,ncol=3)
for (i in 1:length(names(pca.lm))-3)
{
	lm.1=lm(pca.lm[,1]~pca.lm[,i+3])
	lm.2=lm(pca.lm[,2]~pca.lm[,i+3])
	lm.3=lm(pca.lm[,3]~pca.lm[,i+3])
	
	f1=summary(lm.1)$fstatistic
	p_value.1= pf(as.numeric(f1[1]), as.numeric(f1[2]), as.numeric(f1[3]), lower.tail = FALSE)
	
	f2=summary(lm.2)$fstatistic
	p_value.2= pf(as.numeric(f2[1]), as.numeric(f2[2]), as.numeric(f2[3]), lower.tail = FALSE)
	
	f3=summary(lm.3)$fstatistic
	p_value.3= pf(as.numeric(f3[1]), as.numeric(f3[2]), as.numeric(f3[3]), lower.tail = FALSE)
	
	
	p.matrix[i,1]=p_value.1
	p.matrix[i,2]=p_value.2
	p.matrix[i,3]=p_value.3
}	
colnames(p.matrix)=c('PC1','PC2','PC3')
rownames(p.matrix)=c('tumor_type',col.name)
write.table(p.matrix,'gbm_biospecimen_tumor_lm.csv',quote=F,sep=',')





#get final EM.filter and exclude the effects from covariates
EM.filter #14717*166



EM.final=EM.filter[,sample.01.n] #14717*153 
EM.recurrent=EM.filter[,sample.02]
EM.normal=EM.filter[,sample.11]



sample.01.n=c(grep("-01A",name,value=F),grep("-01B",name,value=F),grep("-01C",name,value=F))
pca.filter=pca.lm[sample.01.n,]






corrected=lm(t(EM.final)~pca.filter[,"year_of_diagnosis"]+pca.filter[,"ethnicity"])$residuals
EM.corrected=t(corrected)



#consenesus clustering
rcc3 = ConsensusClusterPlus(t(corrected),maxK=10,reps=50,pItem=0.8,pFeature=1,
title="gbm_rna_hc",distance="pearson",clusterAlg="hc",seed=1262118388.71279,plot="png")

rna.clusters=rcc3[[5]]$consensusClass


cluster1=c()
cluster2=c()
cluster3=c()
cluster4=c()
cluster5=c()
for (i in 1:length(rna.clusters))
{
	if (rna.clusters[[i]]==1)
	{cluster1=append(cluster1,i)}
	
	if (rna.clusters[[i]]==2)
	{cluster2=append(cluster2,i)}
	
	if (rna.clusters[[i]]==3)
	{cluster3=append(cluster3,i)}
	
	if (rna.clusters[[i]]==4)
	{cluster4=append(cluster4,i)}
	
	if (rna.clusters[[i]]==5)
	{cluster5=append(cluster5,i)}
		
}


cluster1.name=rownames(corrected)[cluster1]
cluster1.rna.data=rna.data[,cluster1.name]
write.table(file='gbm_cluster1.txt',cluster1.rna.data,sep='\t',quote=F)

cluster2.name=rownames(corrected)[cluster2]
cluster2.rna.data=rna.data[,cluster2.name]
write.table(file='gbm_cluster2.txt',cluster2.rna.data,sep='\t',quote=F)

cluster3.name=rownames(corrected)[cluster3]
cluster3.rna.data=rna.data[,cluster3.name]
write.table(file='gbm_cluster3.txt',cluster3.rna.data,sep='\t',quote=F)

cluster4.name=rownames(corrected)[cluster4]
cluster4.rna.data=rna.data[,cluster4.name]
write.table(file='gbm_cluster4.txt',cluster4.rna.data,sep='\t',quote=F)

cluster5.name=rownames(corrected)[cluster5]
cluster5.rna.data=rna.data[,cluster5.name]
write.table(file='gbm_cluster5.txt',cluster5.rna.data,sep='\t',quote=F)



heatmap.name=c(cluster1.name,cluster2.name,cluster3.name,cluster4.name,cluster5.name)

anno_col=c(rep('cluster1',length(cluster1)),rep('cluster2',length(cluster2)),
rep('cluster3',length(cluster3)),rep('cluster4',length(cluster4)),
rep('cluster5',length(cluster5)))

marker=c('BCL3','TGFBI','ITGB1','LOX','COL1A2','VDR','IL6','MMP7',
'HOXD3','ERBB3','CDKN1C','PDGFRA','HDAC2','EPHB1',
'PTPRA','ELOVL2','SOX9','PAX6','CDH4','MEOX2','FGFR3')
anno_row=c(rep('MES',8),rep('PN',6),rep('CL',7))





sample.name=c(cluster1.name,cluster2.name,cluster3.name,cluster4.name,cluster5.name)


sample=cbind(sample.name,anno_col)
write.table(file='cluster_label.txt',sample,sep='\t',quote=F)



#######how to change ensembl id to gene symbol

brain=data.frame(gene_id=rownames(rna.data.tpm.filter),vals=runif(length(rna.data.tpm.filter)))
brain$gene_id=as.character(brain$gene_id)
brain$gene_id=sub("[.][0-9]*","",brain$gene_id)

mart=useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes=brain$gene_id
gene_IDs=getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
               values = genes, mart= mart,useCache=F)


library("org.Hs.eg.db")
symbols=mapIds(org.Hs.eg.db, keys = rownames(rna.data.tpm.filter), keytype = "ENSEMBL", column="SYMBOL")


#tpm.final[,heatmap.name]
#write.table(file='gbm_enrichment.txt',tpm.final[,heatmap.name],sep='\t',quote=F)
#write.table(file='gbm_enrichment.gct',tpm.final[,heatmap.name],sep='\t',quote=F)
#write.table(file='gbm_enrichment_c1.gct',tpm.final[,cluster1.name],sep='\t',quote=F)

write.table(file='gbm_cluster1.txt',tpm.final[,cluster1.name],sep='\t',quote=F)
write.table(file='gbm_cluster2.txt',tpm.final[,cluster2.name],sep='\t',quote=F)
write.table(file='gbm_cluster3.txt',tpm.final[,cluster3.name],sep='\t',quote=F)
write.table(file='gbm_cluster4.txt',tpm.final[,cluster4.name],sep='\t',quote=F)
write.table(file='gbm_cluster5.txt',tpm.final[,cluster5.name],sep='\t',quote=F)

tpm.final=rna.data.tpm.filter
rownames(tpm.final)=symbols

heatmap.matrix=tpm.final[marker,heatmap.name]

heatmap.lg=log2(heatmap.matrix+1)


anno_col=as.data.frame(anno_col)
anno_row=as.data.frame(anno_row)

rownames(anno_col)=colnames(heatmap.lg)
rownames(anno_row)=rownames(heatmap.lg)



pdf('gbm_heatmap_2.pdf')
pheatmap(heatmap.lg, 
            scale="row", 
            #clustering_method="average",
            cluster_rows=FALSE,
            cluster_cols=FALSE, 
            color=colorRampPalette(c("green", "black","red"))(100), 
            annotation_col=anno_col, 
            annotation_row=anno_row, 
            annotation_names_row=TRUE, 
            annotation_names_col=FALSE, 
            show_colnames=FALSE,
            #annotation_colors=ann_colors, 
            border_color = NA)
dev.off()














heatmap.name=c(cluster1.name,cluster2.name,cluster3.name,cluster4.name,cluster5.name)

anno_col=c(rep('cluster1',length(cluster1)),rep('cluster2',length(cluster2)),
rep('cluster3',length(cluster3)),rep('cluster4',length(cluster4)),
rep('cluster5',length(cluster5)))

marker=c('BCL3','TGFBI','ITGB1','LOX','COL1A2','VDR','IL6','MMP7',
'HOXD3','ERBB3','CDKN1C','PDGFRA','HDAC2','EPHB1',
'PTPRA','ELOVL2','SOX9','PAX6','CDH4','MEOX2','FGFR3',
'RPL22L1','H4C11','SEL1L2','HOPX','CCDC182','UBD','SLPI','HMHB1',
'ATP1B4','C16orf71','MYLK3','IQCD','SLC7A3','LINC01206','CFAP53','FBXO15')

anno_row=c(rep('MES',8),rep('PN',6),rep('CL',7),rep('c3_marker',8),
rep('c5_marker',8))




library("org.Hs.eg.db")
symbols=mapIds(org.Hs.eg.db, keys = rownames(rna.data.tpm.filter), keytype = "ENSEMBL", column="SYMBOL")


tpm.final=rna.data.tpm.filter
rownames(tpm.final)=symbols

heatmap.matrix=tpm.final[marker,heatmap.name]

heatmap.lg=log2(heatmap.matrix+1)


anno_col=as.data.frame(anno_col)
anno_row=as.data.frame(anno_row)

rownames(anno_col)=colnames(heatmap.lg)
rownames(anno_row)=rownames(heatmap.lg)



pdf('gbm_heatmap_3.pdf')
pheatmap(heatmap.lg, 
            scale="row", 
            #clustering_method="average",
            cluster_rows=FALSE,
            cluster_cols=FALSE, 
            color=colorRampPalette(c("green", "black","red"))(100), 
            annotation_col=anno_col, 
            annotation_row=anno_row, 
            annotation_names_row=TRUE, 
            annotation_names_col=FALSE, 
            show_colnames=FALSE,
            #annotation_colors=ann_colors, 
            border_color = NA)
dev.off()










#plot trees
library(ape)

hcd=hclust(dist(corrected))
plot(hcd)
plot(as.phylo(hcd), cex = 0.2, type = "unrooted")
plot(as.phylo(hcd), cex = 0.2, type = "fan")


pdf("gbm_filter_tree1.pdf")
plot(hcd)
dev.off()
  
pdf("gbm_filter_tree2.pdf")
plot(as.phylo(hcd), cex = 0.2, type = "unrooted")
dev.off()

pdf("gbm_filter_tree3.pdf")
plot(as.phylo(hcd), cex = 0.2, type = "fan")
dev.off()





















#differentional expression analysis
rownames(corrected)

count.data=rna.data[,rownames(corrected)] #56499*152


symbols=mapIds(org.Hs.eg.db, keys = rownames(rna.data), keytype = "ENSEMBL", column="SYMBOL")
rownames(count.data)=symbols




#create a factor for different clusters
cluster=c()
for (i in cluster1)
{
	cluster[i]="cluster1"
}	

for (i in cluster2)
{
	cluster[i]="cluster2"
}

for (i in cluster3)
{
	cluster[i]="cluster3"
}

for (i in cluster4)
{
	cluster[i]="cluster4"
}

for (i in cluster5)
{
	cluster[i]="cluster5"
}



#covariates
#rownames(corrected)

#which(colnames(EM.filter) %in% rownames(corrected) )

q.number=c()
for (i in 1:length(rownames(corrected)))
{
	for (j in 1:length(colnames(EM.filter)))
	{
		if (colnames(EM.filter)[j]==rownames(corrected)[i])
		{q.number=append(q.number,j)}
		
	}
}

#colnames(EM.filter)[q.number]


'''
cluster.type=corrected[heatmap.name,]
anno_col
type=anno_col
type=factor(type,levels=c('cluster1','cluster2','cluster3','cluster4','cluster5'))
colType=c("CornflowerBlue", "ForestGreen", "DarkOrange", "Orchid", "DarkRed")[type]
pchType <- c(16)

pca.cluster=prcomp(cluster.type)

pdf("gbm_cluster_pca.pdf")
plot(pca.cluster$x,
  col = colType,
  pch = pchType,
  cex = 1.0)
legend(
  "bottomright",
  bty = "n",
  c("cluster1","cluster2","cluster3","cluster4","cluster5"),
  fill = c("CornflowerBlue", "ForestGreen", "DarkOrange", "Orchid", "DarkRed"),
  cex = 2.0)
dev.off()
'''






#normal
sample.de.11=grep("-11A",colnames(EM.filter),value=T)

normal.data=rna.data[,sample.de.11] #56499*5
rownames(normal.data)=symbols

normal=c("normal","normal","normal","normal","normal")






#do DE 
#cluster1
data=cbind(count.data[,cluster5],normal.data)

normal.column=c(cluster[cluster5],normal)




colData=data.frame(normal.column)





dds=DESeqDataSetFromMatrix(countData=as.data.frame(data), colData=colData, 
design=~normal.column)
dds=DESeq(dds)
res=results(dds,contrast = c("normal.column","cluster5","normal"))
resOrdered=res[order(res$padj),]
#summary(resOrdered)
#print(sum(resOrdered$padj<0.05 & abs(resOrdered$log2FoldChange)>1.5, na.rm=TRUE))
#q=sum(resOrdered$padj<0.01 & abs(resOrdered$log2FoldChange)>1, na.rm=TRUE)





cluster5_deg=resOrdered
save(cluster5_deg,file='cluster5_deg.rds')
#length(which(cluster1_deg$padj<0.05 & abs(cluster3_deg$log2FoldChange)>1.5))

cluster5_l1000=as.matrix(cluster5_deg[which(cluster5_deg$padj<0.05 
& abs(cluster5_deg$log2FoldChange)>1.5), c(1,3,6)])
write.table(file="cluster5_l1000.txt",cluster5_l1000,sep="\t",quote=F)


cluster1=as.matrix(cluster1_deg[which(cluster1_deg$padj<0.05 
& abs(cluster1_deg$log2FoldChange)>1.5), c(2,5)])
write.table(file="l1000_c1.txt",cluster1,sep="\t",quote=F)



DEG1=as.data.frame(res)
DEG1=na.omit(DEG1)

logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   

DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))

this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))


g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +   
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   
ggsave(g,filename='gbm_cluster5_volcano.pdf')








'''
data=cbind(count.data[,cluster5],normal.data)

normal.column=c(cluster[cluster5],normal)

'''





#DE cluster3 vs rest and cluster5 vs rest


data=count.data


col_cluster3=c()
for (i in cluster1)
{
	col_cluster3[i]="rest"
}	

for (i in cluster2)
{
	col_cluster3[i]="rest"
}

for (i in cluster3)
{
	col_cluster3[i]="cluster3"
}

for (i in cluster4)
{
	col_cluster3[i]="rest"
}

for (i in cluster5)
{
	col_cluster3[i]="rest"
}


colData=data.frame(col_cluster3)


dds=DESeqDataSetFromMatrix(countData=as.data.frame(data), colData=colData, 
design=~col_cluster3)
dds=DESeq(dds)
res=results(dds,contrast = c("col_cluster3","cluster3","rest"))
resOrdered=res[order(res$padj),]




cluster3_rest_deg=resOrdered
save(cluster3_rest_deg,file='cluster3_rest_deg.rds')
#length(which(cluster1_deg$padj<0.05 & abs(cluster3_deg$log2FoldChange)>1.5))

cluster3_rest_deg_filter=as.matrix(cluster3_rest_deg[which(cluster3_rest_deg$padj<0.05 
& abs(cluster3_rest_deg$log2FoldChange)>1.5),])
write.table(file="cluster3_rest_deg_filter.txt",cluster3_rest_deg_filter,sep="\t",quote=F)



DEG1=as.data.frame(res)
DEG1=na.omit(DEG1)

logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   

DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))

this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))


g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +   
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   
ggsave(g,filename='cluster3_rest_volcano.pdf')











#cluster5
data=count.data


col_cluster3=c()
for (i in cluster1)
{
	col_cluster3[i]="rest"
}	

for (i in cluster2)
{
	col_cluster3[i]="rest"
}

for (i in cluster3)
{
	col_cluster3[i]="rest"
}

for (i in cluster4)
{
	col_cluster3[i]="rest"
}

for (i in cluster5)
{
	col_cluster3[i]="cluster5"
}


colData=data.frame(col_cluster3)


dds=DESeqDataSetFromMatrix(countData=as.data.frame(data), colData=colData, 
design=~col_cluster3)
dds=DESeq(dds)
res=results(dds,contrast = c("col_cluster3","cluster5","rest"))
resOrdered=res[order(res$padj),]




cluster5_rest_deg=resOrdered
save(cluster5_rest_deg,file='cluster5_rest_deg.rds')
#length(which(cluster1_deg$padj<0.05 & abs(cluster3_deg$log2FoldChange)>1.5))

cluster3_rest_deg_filter=as.matrix(cluster3_rest_deg[which(cluster3_rest_deg$padj<0.05 
& abs(cluster3_rest_deg$log2FoldChange)>1.5),])
write.table(file="cluster5_rest_deg_filter.txt",cluster5_rest_deg_filter,sep="\t",quote=F)



DEG1=as.data.frame(res)
DEG1=na.omit(DEG1)

logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   

DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))

this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))


g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +   
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   
ggsave(g,filename='cluster5_rest_volcano.pdf')











#DE function 
library(DESeq2)

DE=function(x){
data=cbind(count.data[,x],normal.data)

normal.column=c(cluster[x],normal)


colData=data.frame(normal.column)


dds=DESeqDataSetFromMatrix(countData=data, colData=colData, 
design=~normal.column)
dds=DESeq(dds)
res=results(dds)
resOrdered=res[order(res$padj),]
summary(resOrdered)
print(sum(resOrdered$padj<0.01 & abs(resOrdered$log2FoldChange)>1, na.rm=TRUE))

}














#differentional expression analysis
dim(rna.data)

library(DESeq2)


sample.de.01=c(grep("-01A",colnames(rna.data),value=T),grep("-01B",colnames(rna.data),value=T),grep("-01C",colnames(rna.data),value=T))
sample.de.02=c(grep("-02A",colnames(rna.data),value=T),grep("-02B",colnames(rna.data),value=T))
sample.de.11=grep("-11A",colnames(rna.data),value=T)


tumor.type=c()
for (i in 1:length(colnames(rna.data)))
{
	if (colnames(rna.data)[i] %in% sample.de.01)
	{tumor.type=append(tumor.type,"solid01")}
	if (colnames(rna.data)[i] %in% sample.de.02)
	{tumor.type=append(tumor.type,"recurrent02")}
	if (colnames(rna.data)[i] %in% sample.de.11)
	{tumor.type=append(tumor.type,"normal11")}

}

cancer=c()
for (i in 1:length(colnames(rna.data)))
{
	if (colnames(rna.data)[i] %in% sample.de.01)
	{cancer=append(cancer,"cancer")}
	if (colnames(rna.data)[i] %in% sample.de.02)
	{cancer=append(cancer,"cancer")}
	if (colnames(rna.data)[i] %in% sample.de.11)
	{cancer=append(cancer,"normal")}

}


colData=data.frame(cancer,tumor.type)
rownames(colData)=colnames(rna.data)


dds=DESeqDataSetFromMatrix(countData=rna.data, colData=colData, 
design=~cancer)
dds=DESeq(dds)
res=results(dds)
resOrdered=res[order(res$padj),]
summary(resOrdered)
print(sum(resOrdered$padj<0.01 & abs(resOrdered$log2FoldChange)>1, na.rm=TRUE))

pdf('gbm_ma.pdf')
plotMA(res,ylim=c(-2,2))
dev.off()

pdf('gbm_count.pdf')
plotCounts(dds,gene=which.min(res$padj),intgroup="cancer")
dev.off()




DEG1=as.data.frame(res)
DEG1=na.omit(DEG1)

logFC_cutoff=with(DEG1,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
logFC_cutoff
2^logFC_cutoff   

DEG1$change=as.factor(ifelse(DEG1$pvalue<0.05 &  abs(DEG1$log2FoldChange)>logFC_cutoff,
                             ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))

this_tile=paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                '\nThe number of up gene is ',nrow(DEG1[DEG1$change=='UP',]),
                '\nThe number of down gene is ',nrow(DEG1[DEG1$change=='DOWN',]))


g=ggplot(data=DEG1,
         aes(x=log2FoldChange,y=-log10(pvalue),   
             color=change)) +
  geom_point(alpha=0.4,size=1.75) +     
  theme_set(theme_set(theme_bw(base_size=20))) +
  xlab("log2 fold change")+ylab("-log10 pvalue") +   
  ggtitle(this_tile)+theme(plot.title=element_text(size=15,hjust=0.5)) +
  scale_color_manual(values=c('blue','black','red'))   
ggsave(g,filename='gbm_volcano.pdf')









dds.2=DESeqDataSetFromMatrix(countData=rna.data, colData=colData, 
design=~tumor.type)
dds.2=DESeq(dds.2)
res.2=results(dds.2)
resOrdered.2=res.2[order(res.2$padj),]
summary(resOrdered.2)
print(sum(resOrdered.2$padj<0.01 & abs(resOrdered.2$log2FoldChange)>1, na.rm=TRUE))
















#survival analysis
library("survival")
library("survminer")

clinical_gbm

clinical_gbm_matrix_filter

cluster1.name=colnames(count.data)[cluster1]
cluster2.name=colnames(count.data)[cluster2]
cluster3.name=colnames(count.data)[cluster3]
cluster4.name=colnames(count.data)[cluster4]
cluster5.name=colnames(count.data)[cluster5]


save(cluster1.name,file="cluster1_name.rds")
save(cluster2.name,file="cluster2_name.rds")
save(cluster3.name,file="cluster3_name.rds")
save(cluster4.name,file="cluster4_name.rds")
save(cluster5.name,file="cluster5_name.rds")

normal.sample=sample.de.11
clin=clinical_gbm_matrix_filter



'''
row.number=which(clin[,159] %in% c(cluster1.name,normal.sample))
cluster1.survival=clin[row.number,c(21,126)]
disease=c(rep("cluster1_gbm",length(cluster1)),rep("normal",5))


cluster.sur=cbind(cluster1.survival,disease)
cluster.sur=as.data.frame(cluster.sur)

cluster.sur$days_to_last_follow_up=as.numeric(cluster.sur$days_to_last_follow_up)
cluster.sur$vital_status=as.numeric(cluster.sur$vital_status)
cluster.sur$disease=as.numeric(cluster.sur$disease)
'''


row1.number=which(clin[,159] %in% cluster1.name)
row2.number=which(clin[,159] %in% cluster2.name)
row3.number=which(clin[,159] %in% cluster3.name)
row4.number=which(clin[,159] %in% cluster4.name)
row5.number=which(clin[,159] %in% cluster5.name)
row.number=c(row1.number,row2.number,row3.number,row4.number,row5.number)




gbm_cluster=c(rep("cluster1_gbm",length(cluster1)),
rep("cluster2_gbm",length(cluster2)),
rep("cluster3_gbm",length(cluster3)),
rep("cluster4_gbm",length(cluster4)),
rep("cluster5_gbm",length(cluster5)))




row.number=c(row2.number,row5.number)
gbm_cluster=c(rep("cluster2_gbm",length(cluster2)),
rep("cluster5_gbm",length(cluster5)))

cluster1.survival=clin[row.number,c(21,126)]


cluster.sur=cbind(cluster1.survival,gbm_cluster)

cluster.sur=as.data.frame(cluster.sur)
cluster.sur$days_to_last_follow_up=as.numeric(cluster.sur$days_to_last_follow_up)
cluster.sur$vital_status=as.numeric(cluster.sur$vital_status)
#cluster.gbm_cluster$gbm_cluster=as.numeric(cluster.sur$gbm_cluster)


fit=survfit(Surv(days_to_last_follow_up,vital_status) ~ gbm_cluster, data=cluster.sur)



pdf('c2_c5_survival.pdf')
ggsurvplot(fit, data = cluster.sur,
           pval = TRUE,
           pval.coord = c(0, 0.03), 
           surv.median.line = "hv", 
           legend.title="Cluster",
           ylab="Cumulative survival (percentage)",xlab = " Time (Days)", 
           censor.shape = 124,censor.size = 2,conf.int = FALSE, 
           break.x.by = 20,
           risk.table = TRUE,tables.height = 0.2,
           tables.theme = theme_cleantable(),
           ggtheme = theme_bw())

dev.off()



