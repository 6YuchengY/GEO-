#PDLIM2
# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples
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
Gene <- "PDLIM2"

# load series and platform data from GEO

gset <- getGEO("GSE86845", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
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
control <- c(1, 0, 1, 0)
exp <- c(0, 1, 0, 1)
# 将向量组合成矩阵
design <- cbind(control, exp)

# 设置列名
colnames(design) <- c("control", "exp")

rownames(design)=colnames(dat)[1:4]


# 差异表达矩阵
contrast.matrix<-makeContrasts(paste0(c("exp","control"),collapse = "-"),levels = design)
contrast.matrix

# limma包做差异分析
#1）step1:lmFit
fit <- lmFit(ex[,1:4],design)
#2)step2:eBayes
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2) #default no trend 默认没有趋势
#3)step3:toTable
tempOutput = topTable(fit2, coef=1, n=Inf)#如果需要第一组A的差异基因，就用coef=1;第二组B的差异基因，coef=2
nrDEG = na.omit(tempOutput) 
head(nrDEG)

# 可以在调用 topTable 函数时使用 logbase 参数指定所需的对数底
# tempOutput = topTable(fit2, coef=1, n=Inf, logbase=10)

condition1 <- nrDEG$P.Value < 0.05
# FC看情况调整，有些数据
condition2 <- nrDEG$logFC>log2(1.5)
condition <- condition1 & condition2
filtered_data <- nrDEG[condition,]
old_row_name <- rownames(filtered_data)

# 将行名变回基因名
# 假设filtered_data的行名为原始编号，ids的第一列是原始编号，第二列是目的编号
# 创建一个新的数据框，用于存储转换后的数据
# ids是从对应的R包获取的

#芯片都这么干
library(biomaRt)
mart <- useMart("ensembl","hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
dataset = listDatasets(mart)
mydataset = useDataset("hsapiens_gene_ensembl",mart = mart)
hg_symbols <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
                    filters = 'ensembl_gene_id', values = gsub('_at','',rownames(filtered_data)), mart = mart)



# 空值放到后面处理
# 这里已经注释好了
# hg_symbols <- hg_symbols[which(hg_symbols$external_gene_name != ""),]
# 假设 filtered_data 的行名是原始名称
# hg_symbols 的两列分别是原始名称（entrezgene_id）和目标名称（external_gene_name）


# 这个是芯片的
# 找到匹配位置
#这一行是_at结尾才要
rownames(filtered_data) <- gsub('_at','',rownames(filtered_data))

matched_indices <- match(rownames(filtered_data), hg_symbols$ensembl_gene_id)
# 获取匹配到的目标名称，如果未找到匹配项，则为 NA
new_rownames <- hg_symbols$external_gene_name[matched_indices]
# 只保留找到匹配项的行
filtered_data_matched <- filtered_data[!is.na(matched_indices), ]
matched_new_rownames <- new_rownames[!is.na(matched_indices)]
# 找出那些目标名称不唯一的行
duplicate_targets <- duplicated(matched_new_rownames) | duplicated(matched_new_rownames, fromLast = TRUE)
# 移除目标名称不唯一的行
filtered_data_unique <- filtered_data_matched[!duplicate_targets, ]
# 现在可以安全地设置唯一的行名
rownames(filtered_data_unique) <- matched_new_rownames[!duplicate_targets]





# 这个是RNA-seq的
# 找到匹配位置
matched_indices <- match(rownames(filtered_data), hg_symbols$entrezgene_id)
# 获取匹配到的目标名称，如果未找到匹配项，则为 NA
new_rownames <- hg_symbols$external_gene_name[matched_indices]
# 只保留找到匹配项的行
filtered_data_matched <- filtered_data[!is.na(matched_indices), ]
matched_new_rownames <- new_rownames[!is.na(matched_indices)]
# 找出那些目标名称不唯一的行
duplicate_targets <- duplicated(matched_new_rownames) | duplicated(matched_new_rownames, fromLast = TRUE)
# 移除目标名称不唯一的行
filtered_data_unique <- filtered_data_matched[!duplicate_targets, ]
# 现在可以安全地设置唯一的行名
rownames(filtered_data_unique) <- matched_new_rownames[!duplicate_targets]







rownames(filtered_data) <- hg_symbols$external_gene_name








library(illuminaHumanv4.db)
ids <- toTable(illuminaHumanv4SYMBOL)


converted_data <- filtered_data_unique

converted_data <- filtered_data

# 匹配并转换行名
matching_rows <- match(rownames(filtered_data), ids[, 1])
converted_rows <- ids[matching_rows, 2]

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

