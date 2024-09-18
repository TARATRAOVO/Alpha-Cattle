# 安装biomaRt
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# 使用biomaRt查询SNP位置
library(biomaRt)

# 选择SNP数据库
snp_mart <- useEnsembl(biomart = "snp", dataset = "btaurus_snp")

# 使用ss编号进行查询
snp_info <- getBM(
  attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'allele'),
  filters = 'snp_filter',
  values = 'ss105143725',  # 替换为你的ss编号
  mart = snp_mart
)

# 打印结果
print(snp_info)
# 移除重复的打印
# print(snp_info)  # 此行已被移除