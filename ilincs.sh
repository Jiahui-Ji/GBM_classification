cluster1=as.matrix(cluster1_deg[which(cluster1_deg$padj<0.05 
& abs(cluster1_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster1),invert=TRUE)
cluster1=cluster1[number,]

data.c1=as.data.frame(cbind(rownames(cluster1),cluster1[,2],cluster1[,5]))
write.table(file="l1000_c1.txt",data.c1,sep="\t",quote=F)





cluster2=as.matrix(cluster2_deg[which(cluster2_deg$padj<0.05 
& abs(cluster2_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster2),invert=TRUE)
cluster2=cluster2[number,]

data.c2=as.data.frame(cbind(rownames(cluster2),cluster2[,2],cluster2[,5]))
write.table(file="l1000_c2.txt",data.c2,sep="\t",quote=F)





cluster3=as.matrix(cluster3_deg[which(cluster3_deg$padj<0.05 
& abs(cluster3_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster3),invert=TRUE)
cluster3=cluster3[number,]

data.c3=as.data.frame(cbind(rownames(cluster3),cluster3[,2],cluster3[,5]))
write.table(file="l1000_c3.txt",data.c3,sep="\t",quote=F)





cluster4=as.matrix(cluster4_deg[which(cluster4_deg$padj<0.05 
& abs(cluster4_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster4),invert=TRUE)
cluster4=cluster4[number,]

data.c4=as.data.frame(cbind(rownames(cluster4),cluster4[,2],cluster4[,5]))
write.table(file="l1000_c4.txt",data.c4,sep="\t",quote=F)






cluster5=as.matrix(cluster5_deg[which(cluster5_deg$padj<0.05 
& abs(cluster5_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster5),invert=TRUE)
cluster5=cluster5[number,]

data.c5=as.data.frame(cbind(rownames(cluster5),cluster5[,2],cluster5[,5]))
write.table(file="l1000_c5.txt",data.c5,sep="\t",quote=F)






cluster3_rest=as.matrix(cluster3_rest_deg[which(cluster3_rest_deg$padj<0.05 
& abs(cluster3_rest_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster3_rest),invert=TRUE)
cluster3_rest=cluster3_rest[number,]


cluster3_up=cluster3_rest[which(cluster3_rest[,2]>0),]
cluster3_down=cluster3_rest[which(cluster3_rest[,2]<0),]

write.table(file='c3_up.csv',rownames(cluster3_up), sep=',',quote=F)
write.table(file='c3_down.csv',rownames(cluster3_down), sep=',',quote=F)

write.table(file="c3_rest.csv",cluster3_rest,sep=",",quote=F)





cluster5_rest=as.matrix(cluster5_rest_deg[which(cluster5_rest_deg$padj<0.05 
& abs(cluster5_rest_deg$log2FoldChange)>1.5), ])

number=grep('NA',rownames(cluster5_rest),invert=TRUE)
cluster5_rest=cluster5_rest[number,]

cluster5_up=cluster5_rest[which(cluster5_rest[,2]>0),]
cluster5_down=cluster5_rest[which(cluster5_rest[,2]<0),]

write.table(file='c5_up.csv',rownames(cluster5_up), sep=',',quote=F)
write.table(file='c5_down.csv',rownames(cluster5_down), sep=',',quote=F)

write.table(file="c5_rest.csv",cluster5_rest,sep=",",quote=F)








data=read.table('drug.txt',header=T,sep='\t')

venn.plot <- venn.diagram(
  list(Cluster1=data$Cluster1,Cluster1=data$Cluster2,Cluster3=data$Cluster3,Cluster4=data$Cluster4,Cluster5=data$Cluster5),
  filename = "out5venn.tiff",
  lty = "dotted",
  lwd = 2,
  col = "black",  #"transparent",
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  alpha = 0.60,
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)







#t-cell sceinc

ex0=read.table("exprMatD0.txt",header=T,sep='\t')
ex0=as.matrix(ex0)

ind = apply(ex0, 1, function(x) all(is.na(x)))

expr0=ex0[!ind,]









#output barchart
options(stringsAsFactors=F) 
library(ggplot2) 
library(gplots)

m2=read.table('c3_go.csv',header=T,sep=',')
head(m2)

inp=m2
inp$Stastistics=as.numeric(m2$Stastistics)
inp$logP=-log10(inp$Stastistics)

inp=inp[order(inp$logP,decreasing = F),]
head(inp)

pdf("c3_go.pdf", w=5, h=3.5)
par(las=1)

par(mar=c(5,20,4,0.5))

barplot2(inp$logP,horiz = T,names.arg=inp$PathwayName,xlim=c(0,70)
  , cex.names = 0.7,cex.axis=0.8, border = NA
, xlab = "")
abline(v=1.3013, col='blue', lty=2) # marks 0.05 FDR 
abline(v=2, col='red', lty=2) # marks 0.01 FDR 
dev.off()






rna.clusters

pr5=pam(rna.clusters)

corrected.t=t(corrected)
dis=dist(corrected.t)

si = silhouette(rcc3[[5]]$consensusClass,dis)
ssi = summary(si)

pdf('gbm_sihouette.pdf')
plot(si)
plot(si, col = c("red", "green","orange", "blue", "purple"))
dev.off()






sil=silhouette_SimilarityMatrix(rcc3[[5]]$consensusClass,rcc3[[5]]$consensusMatrix)


pdf('gbm_sihouette.pdf')
plot(sil)
plot(sil, col = c("CornflowerBlue","LightCoral","ForestGreen","DarkOrange","grey"))
dev.off()











#drug specificcity
clList = list()

files = list.files()[grep("csv",list.files())]


for (i in files){
    
    clList[[gsub("[.]csv","", i)]]<- read.csv(i, header = T)
    
}




# calculate relative ranks

for(i in 1:length(clList)){ 

    tmpMat = matrix(nrow = length(nms), ncol = 5); rownames(tmpMat) = nms
    
    nms =  clList[[i]]$ID
    select.k = c(1:5)
    
        for(j in 1:length(nms)){
            
            for (k in c(select.k[select.k != i])){
            
                tmpMat[j,k]<- clList[[i]]$Rank[clList[[i]]$ID == nms[j]]/ clList[[k]]$Rank[clList[[k]]$ID == nms[j]]
 
            }   
        
        }

    tmpMat = as.data.frame(tmpMat)
    tmpMat$mean = apply(tmpMat,1, mean, na.rm=TRUE)

    clList[[i]] = cbind(clList[[i]],tmpMat)

    }







tmp1 = clList[[1]][order(clList[[1]]$mean),]

pdf('cluster1-_v2.pdf')
options(repr.plot.width=4, repr.plot.height=4)

plot(tmp1$Score[tmp1$mean< 1], tmp1$mean[tmp1$mean< 1],
xlab="CMap score",ylab="Specificity score",pch=19,
ylim=c(0,1),xlim=c(-100,0))

abline(h=0.3)
abline(v=-90)
dev.off()



data=tmp1[tmp1$mean< 1,]
pdf('cluster1-_v2.pdf')


scale_y_continuous(expand=c(0,0))





tmp2 = clList[[2]][order(clList[[2]]$mean),]

pdf('cluster2-_v2.pdf')
options(repr.plot.width=4, repr.plot.height=4)

plot(tmp2$Score[tmp2$mean< 1], tmp2$mean[tmp2$mean< 1],
xlab="CMap score",ylab="Specificity score",pch=19,
ylim=c(0,1),xlim=c(-100,0))

abline(h=0.3)
abline(v=-90)
dev.off()


tmp3 = clList[[3]][order(clList[[3]]$mean),]

pdf('cluster3-_v2.pdf')
options(repr.plot.width=4, repr.plot.height=4)

plot(tmp3$Score[tmp3$mean< 1], tmp3$mean[tmp3$mean< 1],
xlab="CMap score",ylab="Specificity score",pch=19,
ylim=c(0,1),xlim=c(-100,0))

abline(h=0.3)
abline(v=-90)
dev.off()


tmp4 = clList[[4]][order(clList[[4]]$mean),]

pdf('cluster4-_v2.pdf')
options(repr.plot.width=4, repr.plot.height=4)

plot(tmp4$Score[tmp4$mean< 1], tmp4$mean[tmp4$mean< 1],
xlab="CMap score",ylab="Specificity score",pch=19,
ylim=c(0,1),xlim=c(-100,0))

abline(h=0.3)
abline(v=-90)
dev.off()




tmp5 = clList[[5]][order(clList[[5]]$mean),]

pdf('cluster5-_v2.pdf')
options(repr.plot.width=4, repr.plot.height=4)

plot(tmp5$Score[tmp5$mean< 1], tmp5$mean[tmp5$mean< 1],
xlab="CMap score",ylab="Specificity score",pch=19,
ylim=c(0,1),xlim=c(-100,0))

abline(h=0.3)
abline(v=-90)
dev.off()


tmp2[tmp2$mean<0.2 & tmp2$Score > 80,]





