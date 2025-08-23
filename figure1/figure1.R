library(ggpubr)
library(survival)
library(forestplot)
library(rms)
library(Hmisc)
library(rmda)
library(regplot)
library(ggplot2)
library(survcomp)
library(sva)
library(limma)
library(FactoMineR)
library(factoextra)
library(VennDiagram)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
library(pheatmap)

####Comparison of clinical information####
inputFile="TARGET_clinical.txt"#input file
outFile="TARGET-boxplot.pdf"      #output file
inputFile="GEO_clinical.txt"
outFile="GEO-boxplot.pdf"

#reading a file
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
x=colnames(rt)[2]
y=colnames(rt)[3]
colnames(rt)=c("id","Group","Expression")

#Set up a comparison group
group=levels(factor(rt$Group))
rt$Group=factor(rt$Group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#Draw a boxplot
boxplot=ggboxplot(rt, x="Group", y="Expression", color="Group",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  palette = c("#F1515E","#1DBDE6"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "TARGET")
#stat_compare_means(comparisons = my_comparisons,method =  "wilcox.test")+labs(x='',y='Futime(years)',title = "GES21257")
#save
pdf(file=outFile,width=6,height=5)
print(boxplot)
dev.off()

####multicox analysis####
rm(list = ls())
clrs <- fpColors(box="black",line="black", summary="black")            
rt=read.table("clinical.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)


rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          
pdf(file="forest-multicox.pdf",onefile = FALSE,
    width = 6,             
    height = 4,           
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()

####Nomogram and calibration curve####
#Read the risk input file
rt=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
dd <- datadist(rt)
options(datadist = "dd")

colnames(rt)

#model building

fit <- coxph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, data = rt)              # BpM

pdf(file="BpM_Nomogram.pdf", width=6, height=6)
nom=regplot(fit,
            clickable=F,
            title="",
            points=TRUE,
            droplines=T,
            observation=NULL,
            rank=NULL,
            failtime = c(1,3,5),
            showP = F,
            prfail = F) 
dev.off()
#calibration curves#
#1 year
pdf(file="BpMcalibration.pdf", width=6, height=6)
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=1, m=28,B=1000)##  m = 28 indicates that 20% of the total sample size.
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=3, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5 year
fit1 <- cph(Surv(futime, fustat) ~ rcs(Age, 3) + Gender+Group, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(fit1, cmethod="KM", method="boot", u=5, m=28,B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()


##remove batch####
rm(list = ls())
rt=read.table("mergegroup.txt",sep="\t",header=T,check.names=F)
rt <- t(rt)
rt=as.matrix(rt)

exp <- rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

expr <- data
batch <- c(rep("GSE21257",53),rep("TARGET",86))
tissue <- c(rep("Non_metastasis",19),rep("Metastasis",34),rep("Non_metastasis",65),rep("Metastasis",21))
mode <- model.matrix(~as.factor(tissue))

####remove Batch Effect
limma_expr <- removeBatchEffect(expr,batch = batch,design = mode)

#PCA analysis without removing batch effects
pre.pca <- PCA(t(expr),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             title = "Before Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )

#Batch Effect-Corrected PCA Analysis
combat.pca <- PCA(t(limma_expr),graph = FALSE)
fviz_pca_ind(combat.pca,
             geom= "point",
             title = "After Batchremoving",
             col.ind = batch,
             addEllipses = TRUE,
             legend.title="Group"  )

####venn plot####
plot <- venn.diagram(
  x = list(set1,set2,set3),
  category.names = c("ARG", "TARGET","GSE21257"),
  filename =NULL,  
  output=FALSE,
  # Circle Properties:
  col = "black", 
  lty = 1, 
  lwd = 1, 
  fill = c("#C1C1C1","#7AAACB","#E89874"),
  alpha = 0.60, 
  label.col = "black",
  cex = .5, 
  fontfamily = "serif",
  fontface = "bold",
  
  # Collection Name Attribute:
  cat.col = c("#C1C1C1","#7AAACB","#E89874"),
  cat.cex = .6,
  cat.fontfamily = "serif"
)

####Volcano plot####
#input file
rm(list = ls())
df <- read.csv("limmaTab.csv",header = T,row.names = 1)
head(df)

df$Group <- factor(ifelse(df$P.Value < 0.05 & abs(df$logFC) >= 0.2,
                          ifelse(df$logFC >= 0.2, 'Up','Down'),'Stable'))
df[1:10,1:7]

table(df$Group)
df$gene <- row.names(df)

p <- ggplot(df, aes(x = logFC, y = -log10(P.Value),,colour = Group))+
  geom_point( shape = 19, size=2.5,stroke = 0.5)+
  
  scale_color_manual(values=c( "#1874CD",'gray',"#CD2626"))+
  ylab('-log10 (Pvalue)')+
  xlab('log2 (Fold Change)')+
  labs(title = "No_metastases vs Metastases")+
  #The gene name of the added focus point
  geom_text_repel(
    data = df[df$P.Value < 0.05 & abs(df$logFC) > 0.2,],
    aes(label = gene),
    size = 3.5,
    segment.color = NA )+ 
  geom_vline(xintercept = c(-0.2,0.2),lty = 2, col = "black", lwd = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5)+
  theme_bw(
    base_line_size = 1  
  )+
  guides(fill = guide_legend(override.aes = list(size =3)))+
  theme_bw()+
  theme(
    axis.title.x = element_text(hjust = 0.5),
    legend.position = c(0.08, 0.86)
  )


p + theme(  panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(
              hjust = 0.5, 
              size = 14,   
              face = "bold", 
              vjust = 1.5
            )
)

####Pheatmap####

# Load the expression data and annotation data
DEG <- read.table("diffExp.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_col <- read.table("annotation_col.txt",header=T, sep="\t", check.names=F,row.names = 1)

annotation_row <- read.table("annotation_row.txt",header=T, sep="\t", check.names=F,row.names = 1)
#Convert to factor type
annotation_col[] <- lapply(annotation_col, as.factor)

annotation_row[] <- lapply(annotation_row, as.factor)
#Set the name of the row for the comment
rownames(annotation_row) <- rownames(DEG)

rownames(annotation_col) <- colnames(DEG)

#Set annotation color
annotation_colors <- list(
  Regulation=c("Up"="#FC9F5B","Down"="#25998F"),
  Gender = c("Male"="#c9bc9c",
             "Female"="#4665d9"
  ),
  Fustat = c("Alive"="#ff91c2",
             "Dead"="#3e3a39"),
  Group = c("Metastases"="#911fb4","No_metastases"="#818001")
)

# Draw a heat map
p <- pheatmap(DEG,
              scale = "row",
              cluster_rows = F,
              cluster_cols = F,
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(50),
              show_rownames = TRUE,
              show_colnames =FALSE,
              treeheight_row = 0,
              treeheight_col = 0 )
p

dev.off()

