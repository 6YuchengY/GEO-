library(GEOquery)
library(limma)
library(umap)
library(WebGestaltR)
library(annotate)
library(biomaRt)
library(openxlsx)
library(readxl)
library(writexl)



folder_dir <- "GEO/result/"
Gene <- "DIDO1"

# load series and platform data from GEO

gset <- getGEO("GSE156713", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15207", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

dat <- ex

#进行limma分析需要三个矩阵：
#1.1)表达矩阵(exp); 2)分组矩阵(design);3)差异比较矩阵(contrast.matrix)
#看下实验设计,决定谁是实验组，谁是对照组。
control <- c(0,0,0,1,1,1)
exp <- c(1,1,1,0,0,0)
# 将向量组合成矩阵
design <- cbind(control, exp)

# 设置列名
colnames(design) <- c("control", "exp")
rownames(design)=c(colnames(dat)[1:3],colnames(dat)[4:6])


# 差异表达矩阵
contrast.matrix<-makeContrasts(paste0(c("exp","control"),collapse = "-"),levels = design)
contrast.matrix

# limma包做差异分析
#1）step1:lmFit
fit <- lmFit(cbind(ex[,1:3],ex[,4:6]),design)
#2)step2:eBayes
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2) #default no trend 默认没有趋势
#3)step3:toTable
tempOutput = topTable(fit2, coef=1, n=Inf)#如果需要第一组A的差异基因，就用coef=1;第二组B的差异基因，coef=2
nrDEG = na.omit(tempOutput) 
head(nrDEG)

# 可以在调用 topTable 函数时使用 logbase 参数指定所需的对数底
# tempOutput = topTable(fit2, coef=1, n=Inf, logbase=10)

condition1 <- nrDEG$adj.P.Val < 0.05
# FC看情况调整，有些数据
condition2 <- nrDEG$logFC>log2(1.5)
condition <- condition1 & condition2
filtered_data <- nrDEG[condition,]
old_row_name <- rownames(filtered_data)

# 试试菜鸟团的方法来注释
## GPL15207 [PrimeView] Affymetrix Human Gene Expression Array
library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL15207', destdir=".")
colnames(Table(gpl)) ## [1] 49395 24
probe2symbol=Table(gpl)[,c(1,17)]



converted_data <- filtered_data
# 匹配并转换行名
matching_rows <- match(rownames(filtered_data), probe2symbol[, 1])
converted_rows <- probe2symbol[matching_rows, 2]

# 这个就是我们要输出的第一张表

S1 <- cbind(converted_rows,converted_data)
S1 <- na.omit(S1)
colnames(S1)[1] <- "Gene"
# 创建一个 Excel 工作簿
wb <- createWorkbook()
# 在工作簿中创建一个工作表
addWorksheet(wb, paste0(Gene,"_S1"))
# 将数据写入工作表
writeData(wb, paste0(Gene,"_S1"), S1)
# 保存工作簿为 xls 文件
saveWorkbook(wb, paste(folder_dir,paste0(Gene,".xlsx"),sep = ""), overwrite = TRUE)





# 准备进行通路富集分析
gene_list <- S1[,1]
enrichDatabase <- c("pathway_KEGG","geneontology_Biological_Process_noRedundant","pathway_Reactome","community-contributed_Hallmark50")[c(1,2,3,4)]
# refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
refFile <- "GEO/genome protein-coding.txt"

pathway_analysis <- WebGestaltR(
  enrichMethod = "ORA",  #目前只用这个
  organism = "hsapiens", #这个是人,mmusculus是老鼠,listOrganism()看所有
  enrichDatabase = enrichDatabase, #
  interestGene = gene_list,
  interestGeneType = "genesymbol",
  referenceGeneFile = refFile,
  referenceGeneType = "entrezgene_protein-coding",
  isOutput=TRUE,
  outputDirectory = folder_dir,
  projectName = Gene
)


addWorksheet(wb, paste0(Gene,"_S2"))
# 将数据写入工作表
writeData(wb, paste0(Gene,"_S2"), pathway_analysis)
# 保存工作簿为 xls 文件
saveWorkbook(wb, paste(folder_dir,paste0(Gene,".xlsx"),sep = ""), overwrite = TRUE)

