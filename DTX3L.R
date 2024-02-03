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
Gene <- "DTX3L"

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE133876", "file=GSE133876_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts >     
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)
dat <- ex

# #进行limma分析需要三个矩阵：
# #1.1)表达矩阵(exp); 2)分组矩阵(design);3)差异比较矩阵(contrast.matrix)
# #看下实验设计,决定谁是实验组，谁是对照组。
# control <- c(0,0,0,1,1,1)
# exp <- c(1,1,1,0,0,0)
# # 将向量组合成矩阵
# design <- cbind(control, exp)
# 
# # 设置列名
# colnames(design) <- c("control", "exp")
# rownames(design)=c(colnames(dat)[1:3],colnames(dat)[10:12])
# 
# 
# # 差异表达矩阵
# contrast.matrix<-makeContrasts(paste0(c("exp","control"),collapse = "-"),levels = design)
# contrast.matrix
# 
# # limma包做差异分析
# #1）step1:lmFit
# fit <- lmFit(cbind(dat[,1:3],dat[,10:12]),design)
# #2)step2:eBayes
# fit2<-contrasts.fit(fit,contrast.matrix)
# fit2<-eBayes(fit2) #default no trend 默认没有趋势
# #3)step3:toTable
# tempOutput = topTable(fit2, coef=1, n=Inf)#如果需要第一组A的差异基因，就用coef=1;第二组B的差异基因，coef=2
# nrDEG = na.omit(tempOutput) 
# head(nrDEG)
# 
# # 可以在调用 topTable 函数时使用 logbase 参数指定所需的对数底
# # tempOutput = topTable(fit2, coef=1, n=Inf, logbase=10)
# 
# condition1 <- nrDEG$adj.P.Val < 0.05
# # FC看情况调整，有些数据
# condition2 <- nrDEG$logFC>log10(1.5)
# condition <- condition1 & condition2
# filtered_data <- nrDEG[condition,]
# old_row_name <- rownames(filtered_data)

# 老方法
dat_processed <- dat[(rowSums(is.na(dat)) / (ncol(dat))) * 2 <= 0.24, ]
exp_group <- c(1,2,3)
con_group <- c(10,11,12)
log2_fc <- c()
p_value <-  c()


for (row in seq(1, nrow(dat_processed))) {
  tryCatch({
    # 样本量很小就用平均数替代中位数
    fc_temp <- mean(dat_processed[row, exp_group] - dat_processed[row, con_group])
    log2_fc <- c(log2_fc, fc_temp)
    
    # 算完FC算p值，强行做t检验？
    p_value_temp <- t.test(dat_processed[row, exp_group], dat_processed[row, con_group])
    p_value <- c(p_value, p_value_temp$p.value)
    
    print(row)
  }, error = function(e) {
    # 捕获到错误时，给 fc_temp 和 p_value_temp 空值
    fc_temp <- NA
    p_value_temp <- NA
    log2_fc <- c(log2_fc, fc_temp)
    p_value <- c(p_value, p_value_temp)
  })
}
# 调整p值
adj_p <- p.adjust(p_value, method = "BH", n = length(p_value))
total_data <- cbind(dat_processed,p_value,adj_p,log2_fc)

df_total_data <- data.frame(total_data)



condition1 <- df_total_data$p_value < 0.05
# FC看情况调整，有些数据
condition2 <- df_total_data$log2_fc>log10(1.5)
condition <- condition1 & condition2
filtered_data <- df_total_data[condition,]
filtered_data <- na.omit(filtered_data)
old_row_name <- rownames(filtered_data)



#RNA-seq
mart <- useMart("ensembl","hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
mart <- useMart("ensembl","hsapiens_gene_ensembl") #镜像出问题了就用这个

dataset = listDatasets(mart)
mydataset = useDataset("hsapiens_gene_ensembl",mart = mart)
hg_symbols <- getBM(attributes=c('entrezgene_id','external_gene_name'), 
                    filters = 'entrezgene_id', values = old_row_name, mart = mart)

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


converted_data <- filtered_data_unique



# 这个就是我们要输出的第一张表

converted_data <- converted_data[, (ncol(converted_data)-2):ncol(converted_data)]

S1 <- cbind(rownames(converted_data),converted_data)
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
gene_list <- rownames(converted_data)
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

