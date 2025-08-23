library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)

setwd("C:\\Users\\ACER\\Desktop\\Autophagy\\13.risk")             
rt=read.table("riskTest.txt",sep="\t",header=T,row.names=1,check.names=F)       
rt=rt[order(rt$riskScore),]                                     

#Draw the risk curve
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore2.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Samples number",
     ylab="Risk score",
     col=c(rep("#00d293",lowLength),
           rep("#ebba37",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "Low risk"),bty="n",pch=19,col=c("#ebba37","#00d293"),cex=1.2)
dev.off()

#Draw a survival status map
color=as.vector(rt$fustat)
color[color==1]="#3e3a39"
color[color==0]="#ff91c2"
pdf(file="survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Samples number",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("#3e3a39","#ff91c2"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

####Pheatmap####
DEG <- read.table("exp.txt",header=T, sep="\t", check.names=F,row.names = 1)
annotation_col <- read.table("annotation_col.txt",header=T, sep="\t", check.names=F,row.names = 1)
DEG <- t(DEG)
DEG<- DEG[match(rownames(annotation_col), rownames(DEG)), ]#让DEG按annotation_col行名排列
DEG <- t(DEG)
#Convert to factor type
annotation_col[] <- lapply(annotation_col, as.factor)


# Create annotation colors
annotation_colors <- list(
  risk=c("high"="#ebba37","low"="#00d293"),
  Gender = c("Male"="#c9bc9c","Female"="#4665d9"),
  Fustat = c("Alive"="#ff91c2",
             "Dead"="#3e3a39"),
  Group = c("Metastases"="#911fb4","No_metastases"="#818001")
)

# Draw a heat map
pdf("model_geneheatmap3.pdf",width = 8.5,height =3.9)
p <- pheatmap(DEG,
              scale = "row",
              cluster_rows = T,
              cluster_cols = F,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(50),
              show_rownames = TRUE,
              show_colnames =FALSE,
              treeheight_row =20,  
              treeheight_col = 0 )

dev.off()

####survive plot####
inputFile="riskTest.txt"       
survFile="test-survival.pdf"        
rocFile="test-ROC.pdf"             

inputFile="riskTrain.txt"       
survFile="train-survival.pdf"        
rocFile="train-ROC.pdf"  


rt=read.table(inputFile,header=T,sep="\t")

diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}

#Sur Plot
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=T,
                   pval=pValue,
                   pval.size=5,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="Risk",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c('#EBBA37',"#8AD293"),
                   risk.table.height=.25,
                   title = "Merged Test")
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()

###Roc Plot
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=7,height=5)
plot(ROC_rt,time=1,col="#8AD293",title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='#EBBA37',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col="#4665D9",add=TRUE,title=FALSE,lwd=2)
title(main = "Merged Test", cex.main = 1.5)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
       col=c("#8AD293",'#EBBA37',"#4665D9"),lwd=2,bty = 'n')
dev.off()


