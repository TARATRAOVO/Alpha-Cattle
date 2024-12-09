# 安装并加载所需的包
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}
library(rentrez)

# 通过 accession number 查询序列信息并提取染色体编号
query_chr <- function(accession) {
  # 使用 entrez_fetch 获取序列信息
  seq_info <- entrez_fetch(
    db = "nucleotide",
    id = accession,
    rettype = "gb",
    retmode = "text"
  )

  # 使用正则表达式提取染色体信息
  chr_match <- regexpr("/chromosome=\"[^\"]+\"", seq_info)
  if (chr_match > 0) {
    chr_info <- regmatches(seq_info, chr_match)
    chr_num <- gsub("/chromosome=\"([^\"]+)\"", "\\1", chr_info)
    return(chr_num)
  } else {
    return(NA)
  }
}

# 读取输入文件
input_data <- read.table("/srv/shiny-server/liftover/output.bed", sep = "\t", header = FALSE)
accessions <- input_data$V1

# 查询每个序列的染色体信息并直接写入文件
output_file <- "chromosome_results.txt"
file.create(output_file) # 创建新文件

for (acc in accessions) {
  chr <- query_chr(acc)
  if (!is.na(chr)) {
    cat(chr, "\n", file = output_file, append = TRUE)
  } else {
    cat("NA\n", file = output_file, append = TRUE)
  }
}
