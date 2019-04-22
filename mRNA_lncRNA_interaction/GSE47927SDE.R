library(limma)
setwd("F:\\yan\\GSE47927\\HSC_DEG") 
exp<-read.table("MEP_normal_bc_mRNA.txt",sep="\t",header=T,check.names=F)
rownames(exp)<-exp[,1]
exp<-exp[2:6]
exp<-as.matrix(exp)
samps<-factor(c(rep("case",3),rep("control",2)))
design <- model.matrix(~0+samps) ;
colnames(design) <- c("case","control")
fit <- lmFit(exp, design)
cont.matrix<-makeContrasts(case-control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
final<-topTable(fit2, coef=1, number=dim(exp)[1], adjust.method="BH", sort.by="B", resort.by="M")
write.table(final,paste("normal_bc_mRNA.txt"),quote=FALSE,sep="\t")

