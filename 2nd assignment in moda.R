###############################################
# - MODA Course, Microarrary Data Analysis    #
# - 2nd Assignment                            #
# - Gene expression profilling of bladder     #
# - Benign & Malignant tissue (GSE7476)       #
# - 25/11/2024                                #
# - Copyright: Radwa Ragab                    # 
###############################################

#R version 4.4.0 (2024-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 22000)
#Matrix products: default
#locale:
#[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
################################################################################

### Loading libraries  ####
library(readxl)
library(readr)
library(ggfortify)
library(matrixStats)
library("rgl")
library(ComplexHeatmap)
library(genefilter)
library(multtest)
library(Biobase)
library(affy)
library(affyPLM)
#library(hgu133a.db)
library(hgu133plus2.db)
#ibrary(hgu133acdf)
#library(hgu133plus2cdf)
BiocManager::install('EnhancedVolcano')
library('EnhancedVolcano')
############################################################################### 

#### Loading Cel files
celFilesDirectory= "C:/Users/Radwa.DESKTOP-77IVMB8/Documents/MODA/GSE7476"
cels= list.files ("C:/Users/Radwa.DESKTOP-77IVMB8/Documents/MODA/GSE7476", pattern = "CEL")
cels
affyData=ReadAffy(celfile.path = celFilesDirectory)
affyData 

# read the phenotable file 
phenotable = read.csv("C:/Users/Radwa.DESKTOP-77IVMB8/Documents/MODA/GSE7476/phenotable.csv")
phenotable$sample.type2=unlist(sapply(strsplit(phenotable$sample.type, "_"),c)[2,]) #2 for row

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

#Q1:
## pca on exp.raw ##
exp.pca=prcomp(t(exp.raw),scale=T)
autoplot(exp.pca, data = phenotable, colour = 'sample.type2',frame = F)
autoplot(exp.pca, data = phenotable, colour = 'sample.type2',frame = T, frame.type="norm")
autoplot(exp.pca, data = phenotable, colour = 'sample.type',frame = F, label = TRUE, label.size = 2)

# 1.histogram
hist(affyData,main="Histogram of affy data")
cols=seq(1:length(sampleNames(affyData)))
legend("topright", sampleNames(affyData), col=cols,lty=1,lwd=3,cex = 0.3,ncol = 2)

# 2.boxplot
boxplot(affyData,main="Boxplot",col=cols)
####################################################################################
#data pre-processing in one step by eset
#threestep(background correction, normalization, summarization)
eset = threestep(affyData,
                 background.method = "IdealMM",
                 normalize.method = "quantile",
                 summary.method = "average.log")

View(exprs(eset))
range(exprs(eset))

hist(eset, main = "Histogram after pre-processing for affyData")
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
exp.data=data2[-dim(data2)[2]]
exp.data=apply(exp.data,2,as.numeric)

## remove duplication
exp.data.agg=aggregate(exp.data,by=list(data2$symbol),FUN=mean)
names(exp.data.agg)
row.names(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[-1]
names(exp.data.agg)=unlist(sapply(strsplit(names(exp.data.agg),"\\."),c)[1,]) #1 for col

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
f

pval=rowttests(as.matrix(exp), f)$p.value

#adjust p value
correctPvalueandReturnAll <- function(tt.pval, method) {
adj.pval <- p.adjust(tt.pval, method = method)
return(adj.pval)
}
 adj.pval=correctPvalueandReturnAll(pval,"BH")
 
 # Assuming lfc and adj.pval are vectors with the same length
 res <- data.frame(lfc = lfc, adj.pval = adj.pval)
 
 # Check the structure of the resulting data frame
 str(res)
 
 print(res)
##########################################################

### selection criteria for identifying DEGS

res.degs=res[abs(lfc) > log(2) & adj.pval < 0.05 , ]

res.degs=res.degs[order(res.degs[,2]) , ]

degs.genes= row.names(res.degs)
exp.degs=exp[degs.genes,]

str(res.degs)
dim(res.degs)
colnames(res.degs)


# export them for further analysis
write.table(degs.genes,file = "DEGs.txt",row.names = F,col=F,quote =F )

# Ensure res.degs is a data frame
res.degs <- as.data.frame(res.degs)

# Initialize the 'regulation' column as "Down"
res.degs$regulation <- "Down"

# Update the 'regulation' column based on lfc values
res.degs$regulation[res.degs$lfc > 0] <- "Up"

# Check the result
table(res.degs$regulation)


res.degs=as.data.frame(res.degs)
res.degs[ ,"regulation"]="Down"
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
exp.degs=as.matrix(exp.degs)
phenotable.sub=phenotable[phenotable$sample.type2 %in% c(group1,group2),]
column_ha= HeatmapAnnotation(sample.type = phenotable$sample.type2)
Heatmap(exp.degs, row_names_gp = gpar(fontsize=3.5),name = "Exp" ,column_names_gp = gpar(fontsize = 6))
Heatmap(t(scale(t(exp.degs))),row_names_gp = gpar(fontsize = 2.5),name = "Z-score",column_names_gp = gpar(fontsize = 6))

save(exp,phenotable,res,res.degs,res.degs.up,res.degs.down,exp.degs ,file= "GSE7476.RData")
##############################################################################################
# 2D pca on degs:
exp2.pca=prcomp(t(exp.degs),scale=T)
autoplot(exp2.pca, data = phenotable, colour = 'sample.type2',frame = F)
autoplot(exp2.pca, data = phenotable, colour = 'sample.type2',frame = T, frame.type="norm")
autoplot(exp2.pca, data = phenotable, colour = 'sample.type2',frame = F, label = TRUE, label.size = 2)

# 3D pca on degs:
mycolors= rep("lightgreen",dim(phenotable)[1])
mycolors[which(phenotable$sample.type2=="PD")]="red"

plot3d(exp2.pca$x[,1:3], pch=20 )
plot3d(exp2.pca$x[,1:3], pch=30 ,col=mycolors , size = 12)

# volcano plot using the advanced method:
EnhancedVolcano(res.degs,
                lab = rownames(res.degs),
                x = 'lfc',
                y = 'pval',
                title = 'Volcano plot',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

EnhancedVolcano(res.degs,
                lab = rownames(res.degs),
                x = 'lfc',
                y = 'pval')
