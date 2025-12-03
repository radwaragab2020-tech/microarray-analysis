###############################################
# - MODA Course, Microarrary Data Analysis    #
# - lec 4,5                                   #
# - Gene expression profilling of bladder     #
# - Benign & Malignant tissue (GSE7476)       #
# - 22/11/2024                                #
# - Copyright: Radwa Ragab                    # 
###############################################

#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22000)
#Matrix products: default
#locale:
#[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
################################################################################

if(!requireNamespace("BiocManager"))
   install.packages("BiocManager")
BiocManager::install(c('GenomicFeatures','AnnotationDbi', "multtest", "affy","affyPLM","hgu133a.db","hgu133plus2.db","hgu133acdf","hgu133plus2cdf"))
install.packages("org.Hs.eg.db")
if (!requireNamespace("affyPLM", quietly = TRUE)) {
  install.packages("affyPLM")
}
### Library LOading ####
library(readxl)
library(readr)
library(ggfortify)
library(matrixStats)
library(ComplexHeatmap)
library(genefilter)
library(multtest)
library(Biobase)
library(affy)
library(affyPLM) # for threestep function
#library(hgu133a.db)
library(hgu133plus2.db)
#ibrary(hgu133acdf)
#library(hgu133plus2cdf)
###############################################################################

tempdir()
# [1] "C:\Users\XYZ~1\AppData\Local\Temp\Rtmp86bEoJ\Rtxt32dcef24de2"
dir.create(tempdir())


#### Loading Cel files
celFilesDirectory= "D:/MODA/microarray/GSE7476"
cels= list.files ("D:/MODA/microarray/GSE7476", pattern = "CEL")
cels
affyData=ReadAffy(celfile.path = celFilesDirectory)
affyData  

##################Data exploration#########################################
class(affyData)
sampleNames(affyData)
featureNames(affyData)
head(featureNames(affyData))
tail(featureNames(affyData))
annotation(affyData)
dim(affyData)

#See how the RAW expression look like without processing : notice the large value
head(exprs(affyData))
head(Biobase::exprs(affyData))
View(Biobase::exprs(affyData))
Biobase::exprs(affyData)[1:3,1:5]
exp.raw=Biobase::exprs(affyData)

## DO pca on exp.raw ##

# 1.histogram
hist(affyData,main="Histogram of affy data")
cols=seq(1:length(sampleNames(affyData)))
legend("topright", sampleNames(affyData), col=cols,lty=1,lwd=3,cex = 0.3,ncol = 2)

boxplot(affyData,main="Boxplot",col=cols)
####################################################################################
#data pre-processing in one step by eset
#threestep(background correction, normalization, summarization)
eset = threestep(affyData,
                 background.method = "IdealMM",
                 normalize.method = "quantile",
                 summary.method = "average.log")

eset2 = expresso(affyData, bgcorrect.method = "rma",
                 normalize.method = "constant",
                 pmcorrect.method ="pmonly",
                 summary.method = "avgdiff")
# the expresso function doesn't do log transformation.
#don't forget to do it yourself... check the ranges
#View(exprs(eset2))
#range(exprs(eset2))

View(exprs(eset))
View(exprs(eset2))
range(exprs(eset))
range(exprs(eset2))
#eset2=log2(exprs(eset2))

hist(eset, main = "Histogram affyData")
boxplot(eset,main = "Box plot GSE7476",col=seq(1:23))
data=exprs(eset)
View(data)

data=cbind(probe_id=row.names(data),data)
colnames(data)

### mapping probe_id into gene symbol
ls("package:hgu133plus2.db")
hgu133plus2() #give the no. of known probe
mapper = hgu133plus2SYMBOL
mapper

map.df = as.data.frame(mapper)
head(map.df)

# merge the two data frames to have the symbole annotation in the data object
data2=merge(data,map.df,by="probe_id",all.x=T)
head(data2)

# drop the probe id column cause we don't need it again
data2=data2[ ,-1]

# remove nulls from data
data2=data2[! is.na(data2$symbol), ]

# check duplication of the gene symbols
x=duplicated(data2$symbol)
sum(x)

# if yes, you have to aggregate these duplicated genes
exp.data=data2[-dim(data2)[2]] #remove to the last col.(gene symbol)to do agg.
exp.data=apply(exp.data,2,as.numeric)

## remove duplication
exp.data.agg=aggregate(exp.data,by=list(data2$symbol),FUN=mean)
names(exp.data.agg)
row.names(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[-1]
names(exp.data.agg)=unlist(sapply(strsplit(names(exp.data.agg),"\\."),c)[1,]) #1 for col

# read the phenotable file 
phenotable = read.csv("D:/MODA/microarray/GSE7476/phenotable.csv")
phenotable$sample.type2=unlist(sapply(strsplit(phenotable$sample.type, "_"),c)[2,]) #2 for row

# save the final object in a RDATA object
exp=exp.data.agg
save(exp,phenotable, file= "GSE7476.RData")
##################################
###### Do differential expression analysis##########
groups=unique(phenotable$sample.type2)

group1=groups[1]
group2=groups[3]
group1.columns= phenotable[phenotable$sample.type2 ==group1, ]$sample.id
group2.columns= phenotable[phenotable$sample.type2 ==group2, ]$sample.id

exp=exp[ ,c(group1.columns,group2.columns)]

# apply non-specific filter to remove the 25% least variant genes before the differential expression analysis

#calculating the LFC ,1 is rows
lfc=apply(exp,1,function(x) mean(x[group2.columns])-mean(x[group1.columns]))
res=as.data.frame(lfc)

#calculating the pvalue
f=factor(c(rep(1, length(group1.columns)) , rep(2, length(group2.columns)) ))
f  # arrange sample in groups to do t test correctly & know which groups will be compared together(treatment,control)

pval=rowttests(as.matrix(exp), f)$p.value

# another way 
#adjust p value
adj.pval <- p.adjust(pval, method = "BH")

res=cbind(lfc,pval,adj.pval)

# Step 1: Calculate p-values using rowttests
pval <- rowttests(as.matrix(exp), f)$p.value

# Step 2: Define the function to adjust p-values
correctPvalueandReturnAll <- function(tt.pval, method) {
  adj.pval <- p.adjust(tt.pval, method = method)
  return(adj.pval)
}

# Step 3: Adjust the p-values using BH method
adj.pval <- correctPvalueandReturnAll(pval, "BH")

# Step 4: Combine log fold changes, p-values, and adjusted p-values
res <- cbind(lfc, pval, adj.pval)

# Step 5: Check the result
print(res)
##########################################################

### selection criteria for identifying DEGS
#degs.res=res[adj.pval<0.05,]#identify DEGs based on the significance level of variance
#degs.res=res[abs(lfc) > log(2) , ]#identify DEGs based on the lfc only for highexpression

res.degs=res[abs(lfc) > log(2) & adj.pval < 0.05 , ]

res.degs=res.degs[order(res.degs[,3]) , ]# according to col.3 (adj.pvl)

degs.genes= row.names(res.degs) #  names of high expressed genes
exp.degs=exp[degs.genes,] # most important genes

# export them for further analysis
write.table(degs.genes,file = "DEGs.txt",row.names = F,col=F,quote =F )

res.degs=as.data.frame(res.degs)
res.degs[ ,"regulation"]="Down"  # add col. called regulation = down for all rows

res.degs[res.degs$lfc > 0, ]$regulation="Up"

res.degs.up=res.degs[res.degs$regulation=="Up", ]
res.degs.down=res.degs[res.degs$regulation=="Down", ]

### creating volcano plot
plot(res[,1], -10*log10(res[,3]),pch=".",main="Volcano plot",
     xlim=c(-4,4),xlab="log fold-change",
     ylab= "-10*log10(FDR)",cex=1.5
)
abline(h= -10*log10(0.05),col="blue")
abline(v= -log2(2),col="blue")
abline(v= log2(2),col="blue")
text(-3,-10*log10(0.05) -0.3+0.2,"FDR<0.05",col="blue",cex=0.8)
grid()

#highlight only points
points(res.degs.up$lfc,-10*log10(res.degs.up$adj.pval),labels="",pch="o",cex=0.5,col="red")
points(res.degs.down$lfc,-10*log10(res.degs.down$adj.pval),labels="",pch="o",cex=0.5,col="green")

#highlight gene names
text(res.degs.up$lfc,-10*log10(res.degs.up$adj.pval),labels = rownames(res.degs.up))
text(res.degs.down$lfc,-10*log10(res.degs.down$adj.pval),labels = rownames(res.degs.down))
 
# creating a heatmap for the top 100 DEGs genes
top100 = head(res.degs, 100)       
top100.genes = rownames(top100)  
exp.top100 = exp[top100.genes, ]  

annotation_col <- data.frame(Group = factor(c(rep("Control", length(group1.columns)),
                                              rep("Treatment", length(group2.columns)))))
rownames(annotation_col) <- colnames(exp.top100)


if(!requireNamespace("pheatmap", quietly = TRUE)){
  install.packages("pheatmap")
}
library(pheatmap)

pheatmap(exp.top100,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue","white","red"))(100),
         main = "Top 100 DEGs Heatmap")


phenotable.sub=phenotable[phenotable$sample.type2 %in% c(group1,group2),]
column_ha= HeatmapAnnotation(sample.type= phenotable.sub$sample.type2)
Heatmap(exp.degs, row_names_gp = gpar(fontsize=3.5),name = "Exp" ,column_names_gp = gpar(fontsize = 6))
