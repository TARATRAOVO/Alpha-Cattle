# 安装biomaRt
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

# 使用biomaRt查询SNP位置
library(biomaRt)
snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
snp_info <- getBM(
  attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'allele', 'consequence_type'),
  filters = 'snp_filter',
  values = 'SS105143724',  # 输入你的SNP编号
  mart = snp_mart
)
print(snp_info)