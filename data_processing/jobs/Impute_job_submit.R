#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied (path to parameter file)", call. = FALSE)
}

# 加载参数
par_file <- file.path(args[1], "par.RData")
if (!file.exists(par_file)) {
  stop("Parameter file not found: ", par_file)
}
load(par_file)

# 执行填充脚本
script_path <- file.path(res$temp_dir, "run_imputation.sh")
if (!file.exists(script_path)) {
  stop("Imputation script not found: ", script_path)
}
system(script_path)

# 检查结果
imputed_file <- file.path(res$dir, "imputed.vcf.gz")
if (!file.exists(imputed_file)) {
  stop("Imputation failed - output file not found: ", imputed_file)
}

# 创建结果下载链接
current_month <- format(Sys.Date(), "%Y%m")
download_url <- paste0(
  "https://xxx.com/data/results/",
  current_month, "/",
  basename(res$dir), "/",
  "imputed.vcf.gz"
)

# 发送邮件
library(mailR)
send.mail(
  from = "noreply@xxx.com",
  to = res$email,
  subject = "Imputation Job Completed",
  body = paste0(
    "Dear user,\n\n",
    "Your imputation job has been completed successfully.\n\n",
    "Results are available at:\n",
    download_url, "\n\n",
    "Note: The results will be available for 30 days.\n\n",
    "Best regards,\n",
    "HEGS Team"
  ),
  smtp = list(host.name = "smtp.xxx.com", port = 25),
  authenticate = FALSE,
  send = TRUE
)

# 清理临时文件
unlink(res$temp_dir, recursive = TRUE)
