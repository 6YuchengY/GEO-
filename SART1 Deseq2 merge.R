#SART1
library(GEOquery)
library(limma)
library(DESeq2)
library(dplyr)
library(umap)
library(WebGestaltR)
library(annotate)
library(biomaRt)
library(openxlsx)
library(readxl)
library(writexl)


folder_dir <- "GEO/result/"
Gene <- "SART1"

# Differential expression analysis with DESeq2
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE59376", "file=GSE59376_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)


#下载不了，暂时用这个读取一下
destfile <- paste0(folder_dir,"","GSE59376_raw_counts_GRCh38.p13_NCBI.tsv.gz")
tbl <- as.matrix(data.table::fread(destfile, header=T, colClasses="integer"), rownames=1)

# 过滤低表达基因
# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]


# 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
exp_group <- c(9,10)
con_group <- c(1,2)
#要拿出来的数据的索引
total_group <- c(exp_group,con_group)
colnames(tbl)[total_group]#看一下拿出来的组对不对
#分组信息，固定实验组在前，对照组在后
group <- c(rep("exp",2),rep("con",2))
colData <- data.frame(row.names = colnames(tbl[,total_group]),
                      group = group)
colData$group <- factor(colData$group, levels = c("exp", "con"))


# 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
dds <- DESeqDataSetFromMatrix(countData = tbl[,total_group], # 需要分析的表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息

# 构建dds矩阵需要：
# 
# 表达矩阵，即上述代码中的tbl[,total_group]，就是我们前面构建的表达矩阵，
# 行为基因，列为样本，中间为计算reads或者fragment得到的整数。
# 
# 样品信息矩阵，即上述代码中的colData，它的类型是一个dataframe（数据框），
# 行名为样本名，第一列是样品的处理情况（对照还是处理、肿瘤还是正常等），即group，
# condition的类型是一个factor。
# 
# 差异比较矩阵，即上述代码中的design。
# 差异比较矩阵就是告诉差异分析函数是要从要分析哪些变量间的差异，简单说就是说明哪些是对照哪些是处理。

# 进行差异表达分析
dds <- DESeq(dds)
# 查看结果的名称
resultsNames(dds)

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", "exp", "con"))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框,这个就是我们想要的，但要注意fc是log2的
df_total_data <- as.data.frame(resOrdered)



condition1 <- df_total_data$padj < 0.05
# FC看情况调整，有些数据
condition2 <- df_total_data$log2FoldChange>log2(1.5)
condition <- condition1 & condition2
filtered_data <- df_total_data[condition,]
filtered_data <- na.omit(filtered_data)
old_row_name <- rownames(filtered_data)

# 将行名变回基因名
# 假设filtered_data的行名为原始编号，ids的第一列是原始编号，第二列是目的编号
# 创建一个新的数据框，用于存储转换后的数据
# ids是从对应的R包获取的

#RNA-seq
mart <- useMart("ensembl","hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
# mart <- useMart("ensembl","hsapiens_gene_ensembl")

dataset = listDatasets(mart)
mydataset = useDataset("hsapiens_gene_ensembl",mart = mart)
hg_symbols <- getBM(attributes=c('entrezgene_id','external_gene_name'), 
                    filters = 'entrezgene_id', values = old_row_name, mart = mart)


#注释
data_to_be_merge <- tibble::rownames_to_column(df_total_data,"ID")
converted_data <- merge(hg_symbols,data_to_be_merge,by.x = "entrezgene_id", by.y = "ID",all=F)
converted_data <- converted_data[,2:ncol(converted_data)]
converted_data_ordered <- converted_data[order(converted_data$pvalue),]
#去重
converted_data_ordered_uni <-converted_data_ordered[!duplicated(converted_data_ordered$external_gene_name),]


S1 <- cbind(rownames(converted_data_ordered_uni),converted_data_ordered_uni)
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
gene_list <- converted_data_ordered_uni$external_gene_name
enrichDatabase <- c("pathway_KEGG","geneontology_Biological_Process_noRedundant","pathway_Reactome","community-contributed_Hallmark50")[c(1,2,3,4)]
# refFile2 <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
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
  projectName = Gene,
  fdrThr = 1
)
pathway_analysis <- pathway_analysis[pathway_analysis$pValue<0.05,]

addWorksheet(wb, paste0(Gene,"_S2"))
# 将数据写入工作表
writeData(wb, paste0(Gene,"_S2"), pathway_analysis)
# 保存工作簿为 xls 文件
saveWorkbook(wb, paste(folder_dir,paste0(Gene,".xlsx"),sep = ""), overwrite = TRUE)

