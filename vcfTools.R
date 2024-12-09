# 加载必要的库
library(vcfR)
library(stringr)
library(data.table) # 用于更快的数据处理

# 定义一个函数来读取VCF文件并提取基因型矩阵
read_vcf_to_matrix <- function(vcf_file_path) {
  # 验证文件是否存在
  if (!file.exists(vcf_file_path)) {
    stop("VCF文件不存在：", vcf_file_path)
  }

  # 使用tryCatch处理可能的错误
  tryCatch(
    {
      # 读取VCF文件
      vcf <- read.vcfR(vcf_file_path, verbose = FALSE)

      # 添加基本信息验证
      message(sprintf("VCF文件包含 %d 个变异位点", nrow(vcf@fix)))
      message(sprintf("VCF文件包含 %d 个样本", ncol(vcf@gt) - 1))

      # 提取基因型数据（GT字段）
      genotype_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

      # 验证基因型矩阵
      message("基因型矩阵维度:")
      message(paste(dim(genotype_matrix), collapse = " x "))

      # 将基因型数据转换为矩阵
      genotype_matrix <- as.matrix(genotype_matrix)

      # 构造行名并验证
      chr <- getCHROM(vcf)
      pos <- getPOS(vcf)
      rownames(genotype_matrix) <- paste(chr, pos, sep = "_")

      message("前几个位点的染色体号:")
      message(paste(head(chr), collapse = ", "))
      message("前几个位点的位置:")
      message(paste(head(pos), collapse = ", "))

      return(genotype_matrix)
    },
    error = function(e) {
      stop("读取VCF文件时发生错误：", e$message)
    }
  )
}

# 修改计算单个位点REF比例的函数
calculate_single_position_ref_percent <- function(vcf_file_path, chromosome, position) {
  # 构造临时输出文件路径
  output_file <- tempfile(pattern = "position_", fileext = ".vcf")
  
  # 构造tabix命令
  tabix_command <- sprintf("tabix %s %s:%s-%s > %s",
                          vcf_file_path,
                          chromosome,
                          position,
                          position,
                          output_file)
  
  # 执行tabix命令
  system_result <- system(tabix_command)
  
  if (system_result != 0) {
    warning(sprintf("tabix提取位点失败：%s:%s", chromosome, position))
    return(NA)
  }
  
  # 读取提取的VCF文件
  tryCatch({
    vcf_content <- readLines(output_file)
    if (length(vcf_content) == 0) {
      message(sprintf("位点 %s:%s 未找到数据", chromosome, position))
      return(NA)
    }
    
    # 解析VCF行
    fields <- strsplit(vcf_content[length(vcf_content)], "\t")[[1]]
    genotypes <- fields[10:length(fields)]  # 从第10列开始是基因型数据
    
    # 计算REF比例
    ref_count <- sum(genotypes %in% c("0|0", "0/0", "0|1", "1|0", "0/1", "1/0"))
    total_count <- length(genotypes)
    ref_percent <- (ref_count / total_count) * 100
    
    message(sprintf("位点 %s:%s 的REF比例: %.2f%%", chromosome, position, ref_percent))
    return(ref_percent)
    
  }, error = function(e) {
    message(sprintf("处理位点 %s:%s 时发生错误: %s", chromosome, position, e$message))
    return(NA)
  }, finally = {
    # 清理临时文件
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
  })
}

# 修改查询函数
query_ref_percent <- function(vcf_file_path, chromosome, position) {
  message(sprintf("查询位点: %s:%s", chromosome, position))
  result <- calculate_single_position_ref_percent(vcf_file_path, chromosome, position)
  
  if (!is.na(result)) {
    message(sprintf("找到位点，REF比例为: %.2f%%", result))
  } else {
    message(sprintf("位点 %s:%s 未找到或处理失败", chromosome, position))
  }
  
  return(result)
}

# 修改批量查询函数
batch_query_ref_percent <- function(vcf_file_path, query_df) {
  required_cols <- c("Chr.", "start")
  if (!all(required_cols %in% colnames(query_df))) {
    stop("查询数据框必须包含 'Chr.' 和 'start' 列")
  }
  
  results <- sapply(1:nrow(query_df), function(i) {
    query_ref_percent(vcf_file_path, query_df$Chr.[i], query_df$start[i])
  })
  
  return(results)
}

# 使用示例：
# vcf_file <- "path_to_your_vcf_file.vcf"
# genotype_matrix <- read_vcf_to_matrix(vcf_file)
# ref_percent <- calculate_ref_percent(genotype_matrix)
#
# # 单个查询
# specific_ref_percent <- query_ref_percent(ref_percent, "chr1", 12345)
#
# # 批量查询
# snp_info <- read.csv("data_base/meg_snp_info/snp_info.csv")
# batch_results <- batch_query_ref_percent(ref_percent, snp_info)
