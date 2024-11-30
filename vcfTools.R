# 在脚本的开头加载库
library(vcfR)
library(stringr)

# 定义一个函数来读取VCF文件并提取基因型矩阵
read_vcf_to_matrix <- function(vcf_file_path) {
  # 读取VCF文件
  vcf <- read.vcfR(vcf_file_path)
  
  # 提取基因型数据（GT字段）
  genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  
  # 将基因型数据转换为矩阵
  genotype_matrix <- as.matrix(genotype_matrix)
  
  # 返回矩阵
  return(genotype_matrix)
}

# 

# 定义一个函数来计算每个位点的REF百分比
calculate_ref_percent <- function(genotype_matrix) {
  # 获取矩阵的行数（即位点数）
  num_loci <- nrow(genotype_matrix)
  
  # 初始化向量存储每个位点的REF百分比
  ref_percent <- numeric(num_loci)
  
  # 遍历每个位点，计算REF等位基因的百分比
  for (i in 1:num_loci) {
    # 获取当前位点的所有样本基因型
    genotypes <- genotype_matrix[i, ]
    
    # 统计REF等位基因（0）的出现次数，参考等位基因可以是"0|0", "0/0", "0|1", "0/1"
    ref_count <- sum(genotypes == "0|0" | genotypes == "0/0" | genotypes == "0|1" | genotypes == "0/1")
    
    # 计算REF百分比 (参考等位基因出现的次数 / 总样本数 * 100)
    ref_percent[i] <- (ref_count / length(genotypes)) * 100
  }
  
  # 返回REF百分比
  return(ref_percent)
}

# 定义一个函数来查询特定位点的REF百分比
query_ref_percent <- function(ref_percent, chromosome, position) {
  # 构造位点名称
  locus_name <- paste(chromosome, position, sep = "_")
  
  # 检查位点名称是否存在于ref_percent向量的名称中
  if (locus_name %in% names(ref_percent)) {
    return(ref_percent[locus_name])
  } else {
    return(NA)  # 如果位点不在ref_percent中，返回NA
  }
}

# # 示例调用函数
# vcf_file <- "path_to_your_vcf_file.vcf"  # 替换为实际VCF文件路径
# genotype_matrix <- read_vcf_to_matrix(vcf_file)

# # 计算每个位点的REF百分比
# ref_percent <- calculate_ref_percent(genotype_matrix)

# # 查询特定位点的REF百分比
# specific_ref_percent <- query_ref_percent(ref_percent, "chr1", 12345)
# print(specific_ref_percent)