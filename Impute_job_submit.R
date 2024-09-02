library(data.table)
library(pedigree)
library(tidyverse)
library(emayili)
library(stringr)

options <- commandArgs(trailingOnly = TRUE)
path = options[1]
breed<-"Holstein"

#定义变量
#path<-"/disk195/zz/shinyApp/HEGS/demo/tb_demo"
#path<-"/home/zhaow/HEGS/joint_breed_demo"
setwd(path)

report_path<-"/srv/shiny-server"
load("par.RData")

dir<-res$dir
res$gen_file<-gsub(".chip","",res$gen_file) ###先去掉demo里geno文件名里的chip
gen_file<-res$gen_file
density<-res$density

plink<-"/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink"
from_email<-"zzsjtu1988@163.com"
to_email<-res$email
task<-res$task

fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
    ### 过滤SNP MAF:0.01 HWE:1e-10
    num_of_variant<-rep(0,length(geno_names))
    for (i in 1:length(geno_names)){
      system(paste0(plink," --bfile ",
                    geno_names[i],
                    " --maf 0.01 --hwe 1e-10 --out ",geno_names[i],".cln --make-bed"))
      
      con  <- file(paste0(geno_names[i],".cln.log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant.+people pass)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant.+people pass)"))
        }
        
      }
      
      close(con)
    }
    
    num_of_variant_final<-sum(num_of_variant)
    ###提取出panel里的个体和位点
    ###分染色体的时候应该命名为chr1 chr2
    ###准备一个id和品种对应的文件  也可以直接每个品种准备一套数据包括表型基因型系谱
    ref_names=vector()
    if(density=="150K_SNP"){
    for(i in 1:length(geno_names)){
      Success=0
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    geno_names[i],
                    " --export vcf  --out ",geno_names[i]))
      system(paste0("bgzip ",geno_names[i],".vcf"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
      ref_names[i]=str_extract(geno_names[i], "chr\\d")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.vcf.gz"))
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".conform.vcf.gz ne=100 chrom=",i," out=",geno_names[i],".beagle.150K"))
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --vcf ",geno_names[i],".beagle.150K.vcf.gz --make-bed --out ",geno_names[i],".beagle.150K"))
      setwd(path)
    }
   }else if(density=="Whole_Genome"){
    for(i in 1:length(geno_names)){
      Success=1
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    geno_names[i],
                    " --export vcf  --out ",geno_names[i]))
      system(paste0("bgzip ",geno_names[i],".vcf"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
      ref_names[i]=str_extract(geno_names[i], "chr\\d")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform.150K"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.150K.vcf.gz"))
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".conform.150K.vcf.gz ne=100 chrom=",i," out=",geno_names[i],".beagle.150K"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".beagle.150K.vcf.gz"))
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=/public/home/chenjc/ltw/nyjh/fill/imputed1000/genotype_beagle_",ref_names[i],".vcf gt=",geno_names[i],".beagle.150K.vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform.WGS"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.WGS.vcf.gz"))
      system(paste0("java -Xmx60g -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=/public/home/chenjc/ltw/nyjh/fill/imputed1000/genotype_beagle_",ref_names[i],".vcf gt=",geno_names[i],".conform.WGS.vcf.gz ne=100 chrom=",i," out=",geno_names[i],".beagle.WGS"))
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --vcf ",geno_names[i],".beagle.WGS.vcf.gz --make-bed --out ",geno_names[i],".beagle.WGS"))
      setwd(path)
    }
  }else{Success=2}


# 打包附件

if(Success==0){
    res_file1=vector()
    res_file2=vector()
    res_file3=vector()
    res_file4=vector()
    for(i in 1:length(geno_names)){
    res_file1[i]=paste0(geno_names[i],".beagle.150K.vcf.gz")
    res_file2[i]=paste0(geno_names[i],".beagle.150K.bed")
    res_file3[i]=paste0(geno_names[i],".beagle.150K.bim")
    res_file4[i]=paste0(geno_names[i],".beagle.150K.fam")
}
    res_file=c(res_file1,res_file2,res_file3,res_file4)
}

if(Success==1){
    res_file1=vector()
    res_file2=vector()
    res_file3=vector()
    res_file4=vector()
    for(i in 1:length(geno_names)){
    res_file1[i]=paste0(geno_names[i],".beagle.WGS.vcf.gz")
    res_file2[i]=paste0(geno_names[i],".beagle.WGS.bed")
    res_file3[i]=paste0(geno_names[i],".beagle.WGS.bim")
    res_file4[i]=paste0(geno_names[i],".beagle.WGS.fam")
}
    res_file=c(res_file1,res_file2,res_file3,res_file4)
}

system(paste0("zip hegs_res.zip ", paste(res_file,collapse = " ")))

# 发送邮件
Success<-as.character(Success)
mail_text<-switch(Success,
                  "0"="Dear user,\n\nYour run of 150K_SNP Imputation has sucessfully finished. Please find results in the attachment. If you have any question, please do not hesitate to contact 1015887056@qq.com.\n\nBest\n\nHEGS Team",
                  "1"="Dear user,\n\nYour run of Whole Genome Imputation has sucessfully finished. Please find results in the attachment. If you have any question, please do not hesitate to contact 1015887056@qq.com.\n\nBest\n\nHEGS Team",
                  "2"="Dear user,\n\nYour run of Imputation has run into some errors. Please check your data format once more. If you have any question, please do not hesitate to contact 1015887056@qq.com.\n\nBest\n\nHEGS Team")

email <- envelope(
  to = to_email,
  from = from_email,
  subject = paste0("HEGS report for job No.", task),
  text = mail_text
) %>% 
  attachment(path="hegs_res.zip")


smtp <- server(host = "smtp.163.com",
               port = 25,
               username = from_email,
               password = "zhangzhe3056436")

tryCatch(
  {
    smtp(email, verbose = F)
    setwd("..")
    system(paste0("rm -rf ", task))
  },
  
  error=function(cond){
    message("[ERROR:] Emailing the results failed!")
    message("Here's the original error message:")
    message(cond)
    return(NA)
  },
  warning=function(cond){
    message("[WARNING:] 	ing the results caused a wanning!\n")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  }
)












