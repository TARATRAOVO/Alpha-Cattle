library(data.table)
library(pedigree)
library(tidyverse)
library(emayili)
library(stringr)
library(plyr)
library(dplyr)

options <- commandArgs(trailingOnly = TRUE)
path = options[1]

#定义变量
#path<-"/disk195/zz/shinyApp/HEGS/demo/tb_demo"
#path<-"/home/zhaow/HEGS/joint_breed_demo"
setwd(path)

report_path<-"/srv/shiny-server"
load("par.RData")
#res$method<-"ssGBLUP"

method<-res$method
model<-res$model
dir<-res$dir
phe_file<-res$phe_file
ped_file<-res$ped_file
res$gen_file<-gsub(".chip","",res$gen_file) ###先去掉demo里geno文件名里的chip
gen_file<-res$gen_file
breed<-"Holstein"
type<-res$type

DMU1<-"/srv/shiny-server/dmu/dmu1"
DMUAI<-"/srv/shiny-server/dmu/dmuai"
DMU4<-"/srv/shiny-server/dmu/dmu4"
plink<-"/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink"
gcta<-"/public/home/chenjc/miniconda3/pkgs/gcta-1.93.2beta-h9ee0642_1/bin/gcta64"
from_email<-"zzsjtu1988@163.com"
to_email<-res$email
task<-res$task

if(type=="Body_Type_Traits"){

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
}


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)

if(method == "ssGBLUP"){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped<-ped[,-4] ###暂时去掉order一列
}

fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)

phe[[1]]=as.character(phe[[1]])
ped[[1]]=as.character(ped[[1]])
ped[[2]]=as.character(ped[[2]])
ped[[3]]=as.character(ped[[3]])
fam[[1]]=as.character(fam[[1]])
fam[[2]]=as.character(fam[[2]])

data_summary_check<-TRUE
if(model_check_pass & method == "ssGBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass &method == "GBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           res$num_of_snps,
           nrow(fam))
  )
}else{
  data_summary_check<-FALSE
}

## 检查表型、系谱、基因型文件ID的一致性
check_id<-function(method){
  res<-list()
  msg<-""
  res$phe_not_in_ped<-character(0)
  res$fam_not_in_ped<-character(0)
  res$phe_not_in_fam<-character(0)
  if(method=="GBLUP" ){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" ){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    fam_not_in_ped<-fam[[1]][!fam[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    res$fam_not_in_ped<-fam_not_in_ped
    if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)==0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. All genotypes are in pedigree.\n", length(phe_not_in_ped))
    }else if(length(phe_not_in_ped)==0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file. All individuals with phenotype are in pedigree.\n", length(fam_not_in_ped))
    }else if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. %d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file.\n", 
                   length(phe_not_in_ped), length(fam_not_in_ped))
      
    }else{
      msg<-"all individuals with phenotype and all genotyped individuals are in pedigree.\n"
    }
  }
  
  res$msg<-msg
  
  if(length(res$phe_not_in_ped)>0){
    write.table(res$phe_not_in_ped,file="phe_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$fam_not_in_ped)>0){
    write.table(res$fam_not_in_ped,file="fam_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$phe_not_in_fam)>0){
    write.table(res$phe_not_in_fam,file="phe_not_in_fam.txt",col.names = F,row.names = F,quote = F)
  }
  
  return(res)
}

if(model_check_pass & data_summary_check){
  check_id_res<-check_id(method)
}

## 处理掉不用的个体
### 不用的表型
phenotype_check_pass<-TRUE
phenotype_summary_msg<-"The summary of each trait in the filtered phenotype data is as follows:\n#+ echo=FALSE\nknitr::kable(phe_summary)"

if(model_check_pass & data_summary_check){
  phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_ped,]
  if(method=="GBLUP" ){
    phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_fam,]
  }
  
  
  phe_summary_num<-map_dbl(model$traits,function(x)sum(!is.na(phe[[x]])))
  phe_summary_mean<-map_dbl(model$traits,function(x)mean(phe[[x]],na.rm=T))
  phe_summary_sd<-map_dbl(model$traits,function(x)sd(phe[[x]],na.rm=T))
  
  phe_summary<-data.frame(
    Trait=model$traits,
    Size=phe_summary_num,
    Mean=phe_summary_mean,
    SD=phe_summary_sd
  )
  
  if(nrow(phe)<=0){    ##不用参考群体时要加 | sum(phe_summary_num)==0
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  fam_ref<-read.table(paste0("../",breed,"_ref/pos_merge.fam"),h=F)
  fam_new=setDF(fam_new)
  fam_ref[,1]=as.character(fam_ref[,1])
  fam_ref[,2]=as.character(fam_ref[,2])
  fam_new[,1]=as.character(fam_new[,1])
  fam_new[,2]=as.character(fam_new[,2])
  fam_all<-rbind(fam_new,fam_ref)##后续计算grm时用到
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_all) ##这个ord是gblup直接用
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i], " --make-bed"))
      }
    }
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
    for(i in 1:length(geno_names)){
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    geno_names[i],
                    " --export vcf  --out ",geno_names[i]))
      system(paste0("bgzip ",geno_names[i],".vcf"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
      ref_names[i]=str_extract(geno_names[i], "chr\\d")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.vcf.gz"))
      system(paste0("/srv/bcftools-1.14/bcftools isec -p ",geno_names[i],"_dir ../",breed,"_ref/",ref_names[i],".vcf.gz ",geno_names[i],".conform.vcf.gz"))
      setwd(paste0(geno_names[i],"_dir"))
      system("bgzip 0001.vcf")
      system("bgzip 0002.vcf")
      system("bgzip 0003.vcf")
      system("/srv/bcftools-1.14/bcftools index 0001.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0002.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0003.vcf.gz")
      system("/srv/bcftools-1.14/bcftools concat 0001.vcf.gz 0003.vcf.gz -o 13_concat.vcf")
      system("/srv/bcftools-1.14/bcftools sort 13_concat.vcf -o 13_sort.vcf")
      system("bgzip 13_sort.vcf")
      system("/srv/bcftools-1.14/bcftools index 13_sort.vcf.gz")
      system("/srv/bcftools-1.14/bcftools merge 13_sort.vcf.gz 0002.vcf.gz -o final.vcf")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../../",breed,"_ref/",ref_names[i],".vcf.gz gt=final.vcf ne=100 chrom=",i," out=merge.beagle"))
      system("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --vcf merge.beagle.vcf.gz --make-bed --out merge.beagle")
      setwd(path)
    }
    
    
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ./", geno_names[i],"_dir/merge.beagle --maf 0 --make-grm --autosome-num 29 --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg <- "[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg <- "[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}


## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & method == "ssGBLUP" & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ped_ref<-read.table(paste0("../",breed,"_ref/ped_ref.txt"),h=T)
    ped=setDF(ped)
    ped_ref[,1]=as.character(ped_ref[,1])
    ped_ref[,2]=as.character(ped_ref[,2])
    ped_ref[,3]=as.character(ped_ref[,3])
    ped[,1]=as.character(ped[,1])
    ped[,2]=as.character(ped[,2])
    ped[,3]=as.character(ped[,3])
    ped<-rbind(ped,ped_ref)
    ord<-1:nrow(ped)  ##这个order是ssblup调用的
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]####这样写要求所有系谱中的个体都在第一列中，跟nadiv包里要求的一样
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999

    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F,na="-9999")
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}



# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-c(model$fix_fct,model$rnd_non)  ##先合并再转换成数字
  phe_fct<-map_dfc(fct_vec,function(x)(phe[[x]]))
  
  phe_ref<-read.table(paste0("../",breed,"_ref/growth.txt"),h=T)
  phe_fct_ref<-map_dfc(fct_vec,function(x)(phe_ref[[x]]))
  phe_fct_all<-rbind(phe_fct,phe_fct_ref)
  names(phe_fct_all)<-fct_vec
  
  fct_all<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe_fct_all[[x]])))  ##把HYS效应变为从1开始编码这样
  names(fct_all)<-fct_vec

  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe_fct_all[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })    ##上一步的对应关系
  names(fct_name)<-fct_vec
  
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  rl_ref<-select(phe_ref,c(model$fix_cov,model$traits))
  rl_all<-rbind(rl,rl_ref) ##和本地群体合并
  
  ## 编码遗传随机效应
  phe_id<-select(phe,model$rnd_gen)
  phe_id_ref<-select(phe_ref,model$rnd_gen)
  phe_id=setDF(phe_id)
  phe_id[,1]=as.character(phe_id[,1])
  phe_id_ref[,1]=as.character(phe_id_ref[,1])
  phe_id_all<-rbind(phe_id,phe_id_ref)
  if(method == "ssGBLUP"){

    rnd_gen<-map_dfc(model$rnd_gen,function(x){##把phe的id换成数字
      ord[match(phe_id_all[[x]],ped[[1]])] ##这个ped已经是两个群体的
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{

    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe_id_all[[x]],fam_all[[1]])]##这里的fam要两个群体的geno id
    })
    rnd_gen_name<-data.frame(Name=fam_all[[1]],Code=ord) 
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<- -9999
  fct_all[is.na(fct_all)]<- -9999
  rl_all[is.na(rl_all)]<- -9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_all){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_all){
    size=nrow(fam)
    grm<-fread(infile,h=F)
    grm_mat<-matrix(0,size,size)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    grm_mat<-t(grm_mat)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    gmat<-grm_mat+0.01*diag(nrow(grm_mat))
    ginv<-solve(gmat)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    #fam_id<-1:nrow(fam)
    ginv_dat<-data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
    return(ginv_dat)
  }
      ### ssGBLUP
  if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl_all,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct_all,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 1 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 PGMIX 2 ASCII ../DMU_PED ../GENO_ID ../GMAT 0.2 G-ADJUST

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name),length(rl_tmp),
                  paste(int_name,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov), paste(1:length(model$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### GBLUP
  }else if(method=="GBLUP"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl_all,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct_all,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 31 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 GREL ASCII ../GINV

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name),length(rl_tmp),
                  paste(int_name,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov), paste(1:length(model$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
  }
  
  dmu_log<-map(model$traits,function(x){
    outdir<-toupper(x)
    setwd(outdir)
    dmu1_log<-system(paste0(DMU1," <DMU.DIR>DMU.lst"), intern=T)
    if(!str_detect(method,"VAR$")){
      dmuai_log<-system(paste0(DMUAI," >>DMU.lst"), intern=T)
    }else{
      dmuai_log<-system(paste0(DMU4," >>DMU.lst"), intern=T)
    }
    
    setwd("../")
    return(sum(attr(dmu1_log,'status'),attr(dmuai_log,'status')))
  })
  names(dmu_log)<-model$traits
  
  if(any(dmu_log>0)){
    Success=1
    error_trait<-names(dmu_log)[dmu_log>0]
  }
  
  if(all(dmu_log>0)){
    Success=2
  }
  
  
  ## 先检验SOL文件是否生成，再检验模型是否收敛
  
  if(!str_detect(method,"VAR$")){
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      lst<-paste0(outdir,"/DMU.lst")
      lst_con<-file(lst, open = "r")
      converged<-FALSE
      while (length(oneLine <- readLines(lst_con, n = 1, warn = FALSE)) > 0) {
        if(str_detect(oneLine, "Converged on")){
          converged<-TRUE
        }
      }
      close(lst_con)
      converge_msg<-""
      if(sol_exists & converged){
        converge_msg<-paste0("Variance components estimation for ",x," has successfully converged.\n")
      }
      if(sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & converged){
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg,SOL=sol_exists, CON=converged)
      return(res)
    })
  }
  
  if(Success!=2){
    ## 生成几个结果的表格
    ### Variance component
    
    if(!str_detect(method,"VAR$")){
      var_comp<-as.data.frame(matrix(NA,length(model$traits),length(rnd_code)+2))
      colnames(var_comp)<-c("Trait",int_name[rnd_code],"Residual")
      var_comp$Trait<-model$traits
      
      walk(model$traits,function(x){
        outdir<-toupper(x)
        if(converge_msg$SOL[converge_msg$TRT==x] & converge_msg$CON[converge_msg$TRT==x]){
          var_com<-read.table(paste0(outdir,"/PAROUT"),h=F)
          var_comp[var_comp$Trait==x,2:ncol(var_comp)]<<-var_com$V4
        }
      })
    }else{
      var_comp<-cbind(model$traits,var_comp_input)
      colnames(var_comp)<-c("Trait",int_name[rnd_code],"Residual")
    }
    write.table(var_comp,file="Variance_component.txt",row.names = F,quote=F,sep="\t")
    
    get_estimate<-function(trait){
      res<-list()
      res$fix_eff<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      res$genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      if(length(model$rnd_non)>0){res$non_genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)}
      outdir<-toupper(trait)
      if(converge_msg$SOL[converge_msg$TRT==trait] & converge_msg$CON[converge_msg$TRT==trait]){
        ### 读取solution
        sol<-fread(paste0(outdir,"/SOL"))
        fix_sol<-sol[sol$V1==2,]
        gen_sol<-sol[sol$V1==4,]
        nongen_sol<-sol[sol$V1==3,]
        if(method=="GBLUP"){
        gen_sol<-sol[sol$V1==3,]
        nongen_sol<-sol[sol$V1==4,]}
        if(length(model$fix_cov)>0){cov_sol<-sol[sol$V1==1,]}
        
        ### fixed factor effect
        fixed_fct_name<-int_name[-rnd_code]
        
        fix_fct_sol<-map_dfr(fixed_fct_name,function(x){
          code<-fct_name[[x]]
          sol_dat<-fix_sol[fix_sol$V4==which(fixed_fct_name==x),]
          sol_dat2<-sol_dat[match(code$Code,sol_dat$V5),]
          sol_res<-data.frame(Name=x,Value=code$Name,Code=code$Code,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
          return(sol_res)
        })
        
        ### Covariate effect
        if(length(model$fix_cov)>0){
          cov_name<-model$fix_cov
          fix_cov_sol<-data.frame(Name=cov_name,Value=1,Code=1:length(cov_name),Estimate=cov_sol$V8,SE=cov_sol$V9)
        }else{
          fix_cov_sol<-NULL
        }
        
        ### Genetic effect
        if(method =="ssGBLUP"){
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fam_all[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name[[x]]$Name,Code=fct_name[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          non_genetic_sol<-non_genetic_sol[complete.cases(non_genetic_sol),]
          
        }
        
        
        
        res$fix_eff<-rbind(fix_fct_sol,fix_cov_sol)
        res$genetic_sol<-genetic_sol
        if(length(model$rnd_non)>0){res$non_genetic_sol<-non_genetic_sol}
      }
      return(res)
    }
    
    all_estimate<-map(model$traits,get_estimate)
    names(all_estimate)<-model$traits
    
    fix_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$fix_eff),idcol="Trait")
    fix_eff$Trait<-model$traits[fix_eff$Trait]
    
    gen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$genetic_sol),idcol="Trait")
    gen_eff$Trait<-model$traits[gen_eff$Trait]
    
    if(length(model$rnd_non)>0){
      nongen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$non_genetic_sol),idcol="Trait")
      nongen_eff$Trait<-model$traits[nongen_eff$Trait]
    }
    
    fwrite(fix_eff, file="fixed_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    fwrite(gen_eff, file="genetic_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    if(length(model$rnd_non)>0){fwrite(nongen_eff, file="nongenetic_random_effect_estimate.txt",row.names=F,quote=F,sep="\t")}
  }
}else{
  Success<-2
}


# 输出html格式报告

con<-file(paste0(report_path,"/report_template.R"), open = "r")
con_out<-file("report.R", open = "w")
model_formula<-paste(unlist(lapply(model$traits,function(x)paste0("\\\\text{",x,"}"))),
                     paste0("=",paste0("\\\\text{",paste(c(model$fix_fct,model$xp,model$rnd_gen,model$rnd_non),
                                                         collapse = "} + \\\\text{"),"}")),collapse = "\\\\\\\\")
model_formula<-paste0("The model is as follows:\n#' $$\n#' ",model_formula,"\n#' $$")

data_summary_entry<-"The basic information of your input data is as follows:\n#+ echo=FALSE\nknitr::kable(data_summary)"

if(length(model$rnd_non)>0){
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
#' 
#' 
#' 
#' ### Non-genetic random effects of first 10 records
#' 
#+ echo=FALSE
nongen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))

"
  
}else{
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
"
}

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
  if(model_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] The model must contain at least one fixed and one random effects. Please resubmit your job!\\*\\*", 
                         model_formula)
  }
  if(data_summary_check){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Please ensure that you have properly uploaded pedigree or genotype files!\\*\\*", 
                         data_summary_entry)
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no data check message is shown. Please check your data first!\\*\\*", 
                         paste0("After initial comparison, ",check_id_res$msg))
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no phenotype data summary message is shown. Please check your data first!\\*\\*", 
                         phenotype_summary_msg)
  }
  
  if(data_summary_check & phenotype_check_pass){
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         genotype_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no genotype data summary message is shown. Please check your data first!**")
  }
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass){
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         pedigree_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no pedigree data summary message is shown. Please check your data first!**")
  }
  
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/all.css ."))
knitr::spin("report.R",knit=F)
rmarkdown::render("report.Rmd")
system("rm report.Rmd")

# 打包附件
res_files<-map_chr(model$traits,function(x){
  system(paste0("cp ",toupper(x),"/DMU.lst", " ./",x,".dmu_run.lst"))
  return(paste0(x,".dmu_run.lst"))
})

if(Success!=2){
  if(length(model$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt","nongenetic_random_effect_estimate.txt")
  }else{
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt")
  }
}

system(paste0("zip hegs_res.zip ", paste(res_files,collapse = " ")))

# 发送邮件
Success<-as.character(Success)
mail_text<-switch(Success,
                  "0"="Dear user,\n\nYour run of HEGS has sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "1"=paste0("Dear user,\n\nSome traits (",paste(error_trait,collapse = ", "),") of your run had not sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team"),
                  "2"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team")

email <- envelope(
  to = to_email,
  from = from_email,
  subject = paste0("HEGS report for job No.", task),
  text = mail_text
) %>% 
  attachment(path="report.html") %>%
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
}else if(type=="Milk_Production_Traits"){

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
}


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)
phe=phe[phe[["DIM"]]>4 & phe[["DIM"]]<336]

if(method == "ssGBLUP"){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped<-ped[,-4] ###暂时去掉order一列
}

fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)

phe[[1]]=as.character(phe[[1]])
ped[[1]]=as.character(ped[[1]])
ped[[2]]=as.character(ped[[2]])
ped[[3]]=as.character(ped[[3]])
fam[[1]]=as.character(fam[[1]])
fam[[2]]=as.character(fam[[2]])

data_summary_check<-TRUE
if(model_check_pass & method == "ssGBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass &method == "GBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           res$num_of_snps,
           nrow(fam))
  )
}else{
  data_summary_check<-FALSE
}

## 检查表型、系谱、基因型文件ID的一致性
check_id<-function(method){
  res<-list()
  msg<-""
  res$phe_not_in_ped<-character(0)
  res$fam_not_in_ped<-character(0)
  res$phe_not_in_fam<-character(0)
  if(method=="GBLUP" ){
    phe_not_in_fam<-unique(phe[[1]][!phe[[1]]%in%fam$V1])
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }

  if(method=="ssGBLUP" ){
    phe_not_in_ped<-unique(phe[[1]][!phe[[1]]%in%ped[[1]]])
    fam_not_in_ped<-fam[[1]][!fam[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    res$fam_not_in_ped<-fam_not_in_ped
    if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)==0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. All genotypes are in pedigree.\n", length(phe_not_in_ped))
    }else if(length(phe_not_in_ped)==0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file. All individuals with phenotype are in pedigree.\n", length(fam_not_in_ped))
    }else if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. %d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file.\n", 
                   length(phe_not_in_ped), length(fam_not_in_ped))
      
    }else{
      msg<-"all individuals with phenotype and all genotyped individuals are in pedigree.\n"
    }
  }

  res$msg<-msg
  
  if(length(res$phe_not_in_ped)>0){
    write.table(res$phe_not_in_ped,file="phe_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$fam_not_in_ped)>0){
    write.table(res$fam_not_in_ped,file="fam_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$phe_not_in_fam)>0){
    write.table(res$phe_not_in_fam,file="phe_not_in_fam.txt",col.names = F,row.names = F,quote = F)
  }
  
  return(res)
}

if(model_check_pass & data_summary_check){
  check_id_res<-check_id(method)
}


## 处理掉不用的个体
### 不用的表型
phenotype_check_pass<-TRUE
phenotype_summary_msg<-"The summary of each trait in the filtered phenotype data is as follows:\n#+ echo=FALSE\nknitr::kable(phe_summary)"

if(model_check_pass & data_summary_check){
  phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_ped,]
  if(method=="GBLUP" ){
    phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_fam,]
  }
  
  
  phe_summary_num<-map_dbl(model$traits,function(x)sum(!is.na(phe[[x]])))
  phe_summary_mean<-map_dbl(model$traits,function(x)mean(phe[[x]],na.rm=T))
  phe_summary_sd<-map_dbl(model$traits,function(x)sd(phe[[x]],na.rm=T))
  
  phe_summary<-data.frame(
    Trait=model$traits,
    Size=phe_summary_num,
    Mean=phe_summary_mean,
    SD=phe_summary_sd
  )
  
  if(nrow(phe)<=0){    ##不用参考群体时要加 | sum(phe_summary_num)==0
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  fam_ref<-read.table(paste0("../",breed,"_ref/pos_merge.fam"),h=F)
  fam_new=setDF(fam_new)
  fam_ref[,1]=as.character(fam_ref[,1])
  fam_ref[,2]=as.character(fam_ref[,2])
  fam_new[,1]=as.character(fam_new[,1])
  fam_new[,2]=as.character(fam_new[,2])
  fam_all<-rbind(fam_new,fam_ref)##后续计算grm时用到
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_all) ##这个ord是gblup直接用
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i], " --make-bed"))
      }
    }
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
    for(i in 1:length(geno_names)){
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    geno_names[i],
                    " --export vcf  --out ",geno_names[i]))
      system(paste0("bgzip ",geno_names[i],".vcf"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
      ref_names[i]=str_extract(geno_names[i], "chr\\d")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.vcf.gz"))
      system(paste0("/srv/bcftools-1.14/bcftools isec -p ",geno_names[i],"_dir ../",breed,"_ref/",ref_names[i],".vcf.gz ",geno_names[i],".conform.vcf.gz"))
      setwd(paste0(geno_names[i],"_dir"))
      system("bgzip 0001.vcf")
      system("bgzip 0002.vcf")
      system("bgzip 0003.vcf")
      system("/srv/bcftools-1.14/bcftools index 0001.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0002.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0003.vcf.gz")
      system("/srv/bcftools-1.14/bcftools concat 0001.vcf.gz 0003.vcf.gz -o 13_concat.vcf")
      system("/srv/bcftools-1.14/bcftools sort 13_concat.vcf -o 13_sort.vcf")
      system("bgzip 13_sort.vcf")
      system("/srv/bcftools-1.14/bcftools index 13_sort.vcf.gz")
      system("/srv/bcftools-1.14/bcftools merge 13_sort.vcf.gz 0002.vcf.gz -o final.vcf")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../../",breed,"_ref/",ref_names[i],".vcf.gz gt=final.vcf ne=100 chrom=",i," out=merge.beagle"))
      system("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --vcf merge.beagle.vcf.gz --make-bed --out merge.beagle")
      setwd(path)
    }
    
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ./", geno_names[i],"_dir/merge.beagle --maf 0 --make-grm --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg <- "[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg <- "[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}


## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & method == "ssGBLUP" & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ped_ref<-read.table(paste0("../",breed,"_ref/ped_ref.txt"),h=T)
    ped=setDF(ped)
    ped_ref[,1]=as.character(ped_ref[,1])
    ped_ref[,2]=as.character(ped_ref[,2])
    ped_ref[,3]=as.character(ped_ref[,3])
    ped[,1]=as.character(ped[,1])
    ped[,2]=as.character(ped[,2])
    ped[,3]=as.character(ped[,3])
    ped<-rbind(ped,ped_ref)
    ord<-1:nrow(ped)  ##这个order是ssblup调用的
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]####这样写要求所有系谱中的个体都在第一列中，跟nadiv包里要求的一样
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999

    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F,na="-9999")
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}


# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
  
  ## 编码固定效应中的因子与非遗传随机效应，这一部分提取固定效应，要求growth.txt中的固定效应包含用户上传的，所以要确定名字一致
  fct_vec<-model$fix_fct ##先合并再转换成数字
  phe_fct<-map_dfc(fct_vec,function(x)(phe[[x]]))
  
  phe_ref<-read.table(paste0("../",breed,"_ref/tdm_growth.txt"),h=T)   ##规定growth.txt里有哪些固定效应
  phe_fct_ref<-map_dfc(fct_vec,function(x)(phe_ref[[x]]))
  phe_fct_all<-rbind(phe_fct,phe_fct_ref)  ##如果用户选择了growth.txt中没有的名字，那么这里就会报错
  names(phe_fct_all)<-fct_vec

  fct_all<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe_fct_all[[x]])))   ##转换为数字
  names(fct_all)<-fct_vec

  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe_fct_all[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })
  names(fct_name)<-fct_vec

  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  rl_ref<-select(phe_ref,c(model$fix_cov,model$traits))
  rl_all<-rbind(rl,rl_ref) ##和本地群体合并
  
  ## 编码遗传随机效应
  phe_id<-select(phe,model$rnd_gen)
  phe_id_ref<-select(phe_ref,model$rnd_gen)
  phe_id=setDF(phe_id)
  phe_id[,1]=as.character(phe_id[,1])
  phe_id_ref[,1]=as.character(phe_id_ref[,1])
  phe_id_all<-rbind(phe_id,phe_id_ref)

  ##编码lg01234
  lg<-select(phe,"DIM")   ##DIM不能为空
  lg<-setDF(lg)
  lg[,1]=as.numeric(lg[,1])    ##DIM
  lg[,1][is.na(lg[,1])]=0
  lg[,2]=-1+2*(lg[,1]-5)/(335-5)  ##ts
  lg[,3]=ceiling((lg[,1]-4.99)/15)   ##DIMph
  lg[,4]=0.7071   ##lg0
  lg[,5]=1.2247*lg[,2]
  lg[,6]=(-0.7906+2.3717*lg[,2]^2)
  lg[,7]=3.5^0.5*(0.5*(5*lg[,2]^3-3*lg[,2]))
  lg[,8]=4.5^0.5*(0.125*(35*lg[,2]^4-30*lg[,2]^2+3))
  names(lg)=c("DIM","ts","DIMph","lg0","lg1","lg2","lg3","lg4")
  lg_ref<-select(phe_ref,c("lg0","lg1","lg2","lg3","lg4"))
  lg_ref<-setDF(lg_ref)
  lg_all<-rbind(lg[,4:8],lg_ref)
  
  lg_DIMph=data.frame(lg[,3])
  names(lg_DIMph)="DIMph"
  lg_all_DIMph=rbind(lg_DIMph,setDF(select(phe_ref,"DIMph")))
  fct_all_prenames=names(fct_all)
  fct_all<-cbind(fct_all,lg_all_DIMph)
  names(fct_all)=c(fct_all_prenames,"DIMph")

  rl_lg_all_p1=select(rl_all,model$fix_cov)
  rl_lg_all_p2=lg_all
  rl_lg_all_p3=select(rl_all,model$traits)
  rl_lg_all=cbind(rl_lg_all_p1,rl_lg_all_p2,rl_lg_all_p3)
  
  if(method == "ssGBLUP"){

    rnd_gen<-map_dfc(model$rnd_gen,function(x){##把phe的id换成数字
      ord[match(phe_id_all[[x]],ped[[1]])] ##这个ped已经是两个群体的
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{

    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe_id_all[[x]],fam_all[[1]])]##这里的fam要两个群体的geno id
    })
    rnd_gen_name<-data.frame(Name=fam_all[[1]],Code=ord) 
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<- -9999
  fct_all[is.na(fct_all)]<- -9999
  rl_all[is.na(rl_all)]<- -9999

  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_name1<-c(model$rnd_gen,model$fix_fct,"DIMph")
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))

  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_all){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_all){
    size=nrow(fam)
    grm<-fread(infile,h=F)
    grm_mat<-matrix(0,size,size)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    grm_mat<-t(grm_mat)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    gmat<-grm_mat+0.01*diag(nrow(grm_mat))
    ginv<-solve(gmat)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    #fam_id<-1:nrow(fam)
    ginv_dat<-data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
    return(ginv_dat)
  }
      ### ssGBLUP
  if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct_all,rl_lg_all)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 1 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 PGMIX 2 ASCII ../DMU_PED ../GENO_ID ../GMAT 0.2 G-ADJUST

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name1),length(rl_tmp),
                  paste(int_name1,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    ##paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " ")
    ### GBLUP
  }else if(method=="GBLUP"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct_all,rl_lg_all)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 31 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 GREL ASCII ../GINV

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name1),length(rl_tmp),
                  paste(int_name1,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
  }
  
  dmu_log<-map(model$traits,function(x){
    outdir<-toupper(x)
    setwd(outdir)
    dmu1_log<-system(paste0(DMU1," <DMU.DIR>DMU.lst"), intern=T)
    if(!str_detect(method,"VAR$")){
      dmuai_log<-system(paste0(DMUAI," >>DMU.lst"), intern=T)
    }else{
      dmuai_log<-system(paste0(DMU4," >>DMU.lst"), intern=T)
    }
    
    setwd("../")
    return(sum(attr(dmu1_log,'status'),attr(dmuai_log,'status')))
  })
  names(dmu_log)<-model$traits
  
  if(any(dmu_log>0)){
    Success=1
    error_trait<-names(dmu_log)[dmu_log>0]
  }
  
  if(all(dmu_log>0)){
    Success=2
  }  

  ## 先检验SOL文件是否生成，再检验模型是否收敛

  if(!str_detect(method,"VAR$")){
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      lst<-paste0(outdir,"/DMU.lst")
      lst_con<-file(lst, open = "r")
      converged<-FALSE
      while (length(oneLine <- readLines(lst_con, n = 1, warn = FALSE)) > 0) {
        if(str_detect(oneLine, "Converged on")){
          converged<-TRUE
        }
      }
      close(lst_con)
      converge_msg<-""
      if(sol_exists & converged){
        converge_msg<-paste0("Variance components estimation for ",x," has successfully converged.\n")
      }
      if(sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & converged){
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg,SOL=sol_exists, CON=converged)
      return(res)
    })
  }

  if(Success!=2){
    ## 生成几个结果的表格
    ### Variance component
     if(!str_detect(method,"VAR$")){
      var_comp<-as.data.frame(matrix(NA,length(model$traits),length(rnd_code)*5+2))
      colnames(var_comp)<-c("Trait",paste(int_name[rnd_code[1]],"*","lg0",sep=""),paste(int_name[rnd_code[1]],"*","lg1",sep=""),paste(int_name[rnd_code[1]],"*","lg2",sep=""),
      paste(int_name[rnd_code[1]],"*","lg3",sep=""),paste(int_name[rnd_code[1]],"*","lg4",sep=""),paste(int_name[rnd_code[2]],"*","lg0",sep=""),paste(int_name[rnd_code[2]],"*","lg1",sep=""),
      paste(int_name[rnd_code[2]],"*","lg2",sep=""),paste(int_name[rnd_code[2]],"*","lg3",sep=""),paste(int_name[rnd_code[2]],"*","lg4",sep=""),"Residual")
      var_comp$Trait<-model$traits

       walk(model$traits,function(x){
        outdir<-toupper(x)
        if(converge_msg$SOL[converge_msg$TRT==x] & converge_msg$CON[converge_msg$TRT==x]){
          var_com<-read.table(paste0(outdir,"/PAROUT"),h=F)
          var_com=var_com[var_com[,2]==var_com[,3],]
          var_comp[var_comp$Trait==x,2:ncol(var_comp)]<<-var_com$V4
        }
      })
    }else{
      var_comp<-cbind(model$traits,var_comp_input)
      colnames(var_comp)<-c("Trait",int_name[rnd_code],"Residual")
    }
    write.table(var_comp,file="Variance_component.txt",row.names = F,quote=F,sep="\t")

   get_estimate<-function(trait){
      res<-list()
      res$fix_eff<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      res$genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      if(length(model$rnd_non)>0){res$non_genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)}
      outdir<-toupper(trait)
      if(converge_msg$SOL[converge_msg$TRT==trait] & converge_msg$CON[converge_msg$TRT==trait]){
        ### 读取solution
        sol<-fread(paste0(outdir,"/SOL"))
        fix_sol<-sol[sol$V1==2,]
        gen_sol<-sol[sol$V1==4,]
        nongen_sol<-sol[sol$V1==3,]
        if(method=="GBLUP"){
        gen_sol<-sol[sol$V1==3,]
        nongen_sol<-sol[sol$V1==4,]}
        if(length(model$fix_cov)>0){cov_sol<-sol[sol$V1==1,]}
        
        ### fixed factor effect
        fixed_fct_name<-fct_vec
        
        fix_fct_sol<-map_dfr(fixed_fct_name,function(x){
          code<-fct_name[[x]]
          sol_dat<-fix_sol[fix_sol$V4==which(fixed_fct_name==x),]
          sol_dat2<-sol_dat[match(code$Code,sol_dat$V5),]
          sol_res<-data.frame(Name=x,Value=code$Name,Code=code$Code,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
          return(sol_res)
        })

        ### Covariate effect
        if(length(model$fix_cov)>0){
          cov_name<-model$fix_cov
          fix_cov_sol<-data.frame(Name=cov_name,Value=1,Code=1:length(cov_name),Estimate=cov_sol$V8,SE=cov_sol$V9)
        }else{
          fix_cov_sol<-NULL
        }
        
        ### Genetic effect
        meantdm=function(x){
            mean=0.707*x[5]-0.1261812*x[4]-0.09011002*x[3]-0.1276164*x[2]-0.03028584*x[1]
            return(mean)}
        if(method =="ssGBLUP"){
          genetic_sol<-map_dfr(model$rnd_gen,function(x){

           sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
           sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
           sol_dat3<-aggregate(V8 ~ V5, data = sol_dat2, meantdm)
           sol_res<-data.frame(Name=x,Value=ped[,1],Code=ord,Estimate=sol_dat3$V2)
           return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
            sol_dat3<-aggregate(V8 ~ V5, data = sol_dat2, meantdm)
            sol_res<-data.frame(Name=x,Value=fam_all[[1]],Code=ord,Estimate=sol_dat3$V2)
            return(sol_res)
          })
        }

        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          fct_vec1<-c(model$fix_fct,model$rnd_non)
          phe_fct1<-map_dfc(fct_vec1,function(x)(phe[[x]]))
          phe_fct_ref1<-map_dfc(fct_vec1,function(x)(phe_ref[[x]]))
          phe_fct_all1<-rbind(phe_fct1,phe_fct_ref1)
          names(phe_fct_all1)<-fct_vec1
          fct_name1<-map(fct_vec1,function(x){
            name<-levels(as.factor(phe_fct_all1[[x]]))
            data.frame(Name=name,Code=1:length(name))
            })
          names(fct_name1)<-fct_vec1
          non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
            sol_dat3<-aggregate(V8 ~ V5, data = sol_dat2, meantdm)
            sol_res<-data.frame(Name=x,Value=fct_name1[[x]]$Name,Code=fct_name1[[x]]$Code,
                                Estimate=sol_dat3$V2)
            return(sol_res)
          })
          non_genetic_sol<-non_genetic_sol[complete.cases(non_genetic_sol),]
          
        }

        res$fix_eff<-rbind(fix_fct_sol,fix_cov_sol)  ##fix_cov_sol有warning，第二列会变成NA
        res$genetic_sol<-genetic_sol
        if(length(model$rnd_non)>0){res$non_genetic_sol<-non_genetic_sol}
      }
      return(res)
    }

    all_estimate<-map(model$traits,get_estimate)
    names(all_estimate)<-model$traits
    
    fix_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$fix_eff),idcol="Trait")
    fix_eff$Trait<-model$traits[fix_eff$Trait]
    
    gen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$genetic_sol),idcol="Trait")
    gen_eff$Trait<-model$traits[gen_eff$Trait]
    
    if(length(model$rnd_non)>0){
      nongen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$non_genetic_sol),idcol="Trait")
      nongen_eff$Trait<-model$traits[nongen_eff$Trait]
    }
    
    fwrite(fix_eff, file="fixed_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    fwrite(gen_eff, file="genetic_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    if(length(model$rnd_non)>0){fwrite(nongen_eff, file="nongenetic_random_effect_estimate.txt",row.names=F,quote=F,sep="\t")}
  }
}else{
  Success<-2
}

options (warn = -1)
con<-file(paste0(report_path,"/report_template.R"), open = "r")
con_out<-file("report.R", open = "w")
model_formula<-paste(unlist(lapply(model$traits,function(x)paste0("\\\\text{",x,"}"))),
                     paste0("=",paste0("\\\\text{",paste(c(model$fix_fct,model$xp,model$rnd_gen,"lg0","lg1","lg2","lg3","lg4",model$rnd_non),
                                                         collapse = "} + \\\\text{"),"}")),collapse = "\\\\\\\\")
model_formula<-paste0("The model is as follows:\n#' $$\n#' ",model_formula,"\n#' $$")

data_summary_entry<-"The basic information of your input data is as follows:\n#+ echo=FALSE\nknitr::kable(data_summary)"

if(length(model$rnd_non)>0){
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
#' 
#' 
#' 
#' ### Non-genetic random effects of first 10 records
#' 
#+ echo=FALSE
nongen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))

"
  
}else{
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
"
}

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
  if(model_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] The model must contain at least one fixed and one random effects. Please resubmit your job!\\*\\*", 
                         model_formula)
  }
  if(data_summary_check){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Please ensure that you have properly uploaded pedigree or genotype files!\\*\\*", 
                         data_summary_entry)
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no data check message is shown. Please check your data first!\\*\\*", 
                         paste0("After initial comparison, ",check_id_res$msg))
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no phenotype data summary message is shown. Please check your data first!\\*\\*", 
                         phenotype_summary_msg)
  }
  
  if(data_summary_check & phenotype_check_pass){
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         genotype_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no genotype data summary message is shown. Please check your data first!**")
  }
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass){
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         pedigree_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no pedigree data summary message is shown. Please check your data first!**")
  }
  
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/all.css ."))
knitr::spin("report.R",knit=F)
rmarkdown::render("report.Rmd")
system("rm report.Rmd")

# 打包附件
res_files<-map_chr(model$traits,function(x){
  system(paste0("cp ",toupper(x),"/DMU.lst", " ./",x,".dmu_run.lst"))
  return(paste0(x,".dmu_run.lst"))
})

if(Success!=2){
  if(length(model$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt","nongenetic_random_effect_estimate.txt")
  }else{
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt")
  }
}

system(paste0("zip hegs_res.zip ", paste(res_files,collapse = " ")))

# 发送邮件
Success<-as.character(Success)
mail_text<-switch(Success,
                  "0"="Dear user,\n\nYour run of HEGS has sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "1"=paste0("Dear user,\n\nSome traits (",paste(error_trait,collapse = ", "),") of your run had not sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team"),
                  "2"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team")

email <- envelope(
  to = to_email,
  from = from_email,
  subject = paste0("HEGS report for job No.", task),
  text = mail_text
) %>% 
  attachment(path="report.html") %>%
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
}else if(type=="Reproductive_Traits"){

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
}


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)

if(method == "ssGBLUP"){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped<-ped[,-4] ###暂时去掉order一列
}

fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)

phe[[1]]=as.character(phe[[1]])
ped[[1]]=as.character(ped[[1]])
ped[[2]]=as.character(ped[[2]])
ped[[3]]=as.character(ped[[3]])
fam[[1]]=as.character(fam[[1]])
fam[[2]]=as.character(fam[[2]])

data_summary_check<-TRUE
if(model_check_pass & method == "ssGBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass &method == "GBLUP"){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           res$num_of_snps,
           nrow(fam))
  )
}else{
  data_summary_check<-FALSE
}

## 检查表型、系谱、基因型文件ID的一致性
check_id<-function(method){
  res<-list()
  msg<-""
  res$phe_not_in_ped<-character(0)
  res$fam_not_in_ped<-character(0)
  res$phe_not_in_fam<-character(0)
  if(method=="GBLUP" ){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" ){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    fam_not_in_ped<-fam[[1]][!fam[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    res$fam_not_in_ped<-fam_not_in_ped
    if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)==0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. All genotypes are in pedigree.\n", length(phe_not_in_ped))
    }else if(length(phe_not_in_ped)==0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file. All individuals with phenotype are in pedigree.\n", length(fam_not_in_ped))
    }else if(length(phe_not_in_ped)>0 & length(fam_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file. %d genotyped individual(s) not in pedigree were removed. You can find removed ID in `fam_id_not_in_ped.txt` file.\n", 
                   length(phe_not_in_ped), length(fam_not_in_ped))
      
    }else{
      msg<-"all individuals with phenotype and all genotyped individuals are in pedigree.\n"
    }
  }
  
  res$msg<-msg
  
  if(length(res$phe_not_in_ped)>0){
    write.table(res$phe_not_in_ped,file="phe_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$fam_not_in_ped)>0){
    write.table(res$fam_not_in_ped,file="fam_id_not_in_ped.txt",col.names = F,row.names = F,quote = F)
  }
  if(length(res$phe_not_in_fam)>0){
    write.table(res$phe_not_in_fam,file="phe_not_in_fam.txt",col.names = F,row.names = F,quote = F)
  }
  
  return(res)
}

if(model_check_pass & data_summary_check){
  check_id_res<-check_id(method)
}

## 处理掉不用的个体
### 不用的表型
phenotype_check_pass<-TRUE
phenotype_summary_msg<-"The summary of each trait in the filtered phenotype data is as follows:\n#+ echo=FALSE\nknitr::kable(phe_summary)"

if(model_check_pass & data_summary_check){
  phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_ped,]
  if(method=="GBLUP" ){
    phe<-phe[!phe[[1]]%in%check_id_res$phe_not_in_fam,]
  }
  
  
  phe_summary_num<-map_dbl(model$traits,function(x)sum(!is.na(phe[[x]])))
  phe_summary_mean<-map_dbl(model$traits,function(x)mean(phe[[x]],na.rm=T))
  phe_summary_sd<-map_dbl(model$traits,function(x)sd(phe[[x]],na.rm=T))
  
  phe_summary<-data.frame(
    Trait=model$traits,
    Size=phe_summary_num,
    Mean=phe_summary_mean,
    SD=phe_summary_sd
  )
  
  if(nrow(phe)<=0){    ##不用参考群体时要加 | sum(phe_summary_num)==0
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  fam_ref<-read.table(paste0("../",breed,"_ref/pos_merge.fam"),h=F)
  fam_new=setDF(fam_new)
  fam_ref[,1]=as.character(fam_ref[,1])
  fam_ref[,2]=as.character(fam_ref[,2])
  fam_new[,1]=as.character(fam_new[,1])
  fam_new[,2]=as.character(fam_new[,2])
  fam_all<-rbind(fam_new,fam_ref)##后续计算grm时用到
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_all) ##这个ord是gblup直接用
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i], " --make-bed"))
      }
    }
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
    for(i in 1:length(geno_names)){
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    geno_names[i],
                    " --export vcf  --out ",geno_names[i]))
      system(paste0("bgzip ",geno_names[i],".vcf"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
      ref_names[i]=str_extract(geno_names[i], "chr\\d")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform"))
      system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.vcf.gz"))
      system(paste0("/srv/bcftools-1.14/bcftools isec -p ",geno_names[i],"_dir ../",breed,"_ref/",ref_names[i],".vcf.gz ",geno_names[i],".conform.vcf.gz"))
      setwd(paste0(geno_names[i],"_dir"))
      system("bgzip 0001.vcf")
      system("bgzip 0002.vcf")
      system("bgzip 0003.vcf")
      system("/srv/bcftools-1.14/bcftools index 0001.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0002.vcf.gz")
      system("/srv/bcftools-1.14/bcftools index 0003.vcf.gz")
      system("/srv/bcftools-1.14/bcftools concat 0001.vcf.gz 0003.vcf.gz -o 13_concat.vcf")
      system("/srv/bcftools-1.14/bcftools sort 13_concat.vcf -o 13_sort.vcf")
      system("bgzip 13_sort.vcf")
      system("/srv/bcftools-1.14/bcftools index 13_sort.vcf.gz")
      system("/srv/bcftools-1.14/bcftools merge 13_sort.vcf.gz 0002.vcf.gz -o final.vcf")
      system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../../",breed,"_ref/",ref_names[i],".vcf.gz gt=final.vcf ne=100 chrom=",i," out=merge.beagle"))
      system("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --vcf merge.beagle.vcf.gz --make-bed --out merge.beagle")
      setwd(path)
    }
    
    
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ./", geno_names[i],"_dir/merge.beagle --maf 0 --make-grm --autosome-num 29 --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg <- "[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg <- "[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}


## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & method == "ssGBLUP" & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ped_ref<-read.table(paste0("../",breed,"_ref/ped_ref.txt"),h=T)
    ped=setDF(ped)
    ped_ref[,1]=as.character(ped_ref[,1])
    ped_ref[,2]=as.character(ped_ref[,2])
    ped_ref[,3]=as.character(ped_ref[,3])
    ped[,1]=as.character(ped[,1])
    ped[,2]=as.character(ped[,2])
    ped[,3]=as.character(ped[,3])
    ped<-rbind(ped,ped_ref)
    ord<-1:nrow(ped)  ##这个order是ssblup调用的
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]####这样写要求所有系谱中的个体都在第一列中，跟nadiv包里要求的一样
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999

    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F,na="-9999")
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}



# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-c(model$fix_fct,model$rnd_non)  ##先合并再转换成数字
  phe_fct<-map_dfc(fct_vec,function(x)(phe[[x]]))
  
  phe_ref<-read.table(paste0("../",breed,"_ref/growth_production.txt"),h=T)
  phe_fct_ref<-map_dfc(fct_vec,function(x)(phe_ref[[x]]))
  phe_fct_all<-rbind(phe_fct,phe_fct_ref)
  names(phe_fct_all)<-fct_vec
  
  fct_all<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe_fct_all[[x]])))  ##把HYS效应变为从1开始编码这样
  names(fct_all)<-fct_vec

  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe_fct_all[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })    ##上一步的对应关系
  names(fct_name)<-fct_vec
  
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  rl_ref<-select(phe_ref,c(model$fix_cov,model$traits))
  rl_all<-rbind(rl,rl_ref) ##和本地群体合并
  
  ## 编码遗传随机效应
  phe_id<-select(phe,model$rnd_gen)
  phe_id_ref<-select(phe_ref,model$rnd_gen)
  phe_id=setDF(phe_id)
  phe_id[,1]=as.character(phe_id[,1])
  phe_id_ref[,1]=as.character(phe_id_ref[,1])
  phe_id_all<-rbind(phe_id,phe_id_ref)
  if(method == "ssGBLUP"){

    rnd_gen<-map_dfc(model$rnd_gen,function(x){##把phe的id换成数字
      ord[match(phe_id_all[[x]],ped[[1]])] ##这个ped已经是两个群体的
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{

    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe_id_all[[x]],fam_all[[1]])]##这里的fam要两个群体的geno id
    })
    rnd_gen_name<-data.frame(Name=fam_all[[1]],Code=ord) 
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<- -9999
  fct_all[is.na(fct_all)]<- -9999
  rl_all[is.na(rl_all)]<- -9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_all){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_all){
    size=nrow(fam)
    grm<-fread(infile,h=F)
    grm_mat<-matrix(0,size,size)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    grm_mat<-t(grm_mat)
    grm_mat[upper.tri(grm_mat,diag = T)]<-grm$V4
    gmat<-grm_mat+0.01*diag(nrow(grm_mat))
    ginv<-solve(gmat)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    #fam_id<-1:nrow(fam)
    ginv_dat<-data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
    return(ginv_dat)
  }
      ### ssGBLUP
  if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl_all,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct_all,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 1 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 PGMIX 2 ASCII ../DMU_PED ../GENO_ID ../GMAT 0.2 G-ADJUST

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name),length(rl_tmp),
                  paste(int_name,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov), paste(1:length(model$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### GBLUP
  }else if(method=="GBLUP"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl_all,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct_all,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F,na="-9999")
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 31 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1 1 0 0 0
0
%d 0 %d %s
%d %s
%d %s
0

$VAR_STR 1 GREL ASCII ../GINV

$RESIDUALS ASCII

$DMUAI
10  Emstep       Number of steps before full weight on EM in imet = 2
1.0d-7  Conv_ndelta  Convergence criteria for norm of the update vector 
1.0d-6  Conv_gnorm   Convergence criteria for norm of the gradient vector (AI)
1 Printout     1 -> Solution vector is printed/written to file SOL
0 Fspopt       Time (0) or memory (1) optimised FSPAK
0 P_neval      Restart an analysis from evaluation p_neval.
",method,x,length(int_name),length(rl_tmp),
                  paste(int_name,collapse = " "),
                  paste(rl_tmp,collapse = " "),
                  length(rl_tmp),length(int_code),
                  paste(int_code,collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov), paste(1:length(model$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
  }
  
  dmu_log<-map(model$traits,function(x){
    outdir<-toupper(x)
    setwd(outdir)
    dmu1_log<-system(paste0(DMU1," <DMU.DIR>DMU.lst"), intern=T)
    if(!str_detect(method,"VAR$")){
      dmuai_log<-system(paste0(DMUAI," >>DMU.lst"), intern=T)
    }else{
      dmuai_log<-system(paste0(DMU4," >>DMU.lst"), intern=T)
    }
    
    setwd("../")
    return(sum(attr(dmu1_log,'status'),attr(dmuai_log,'status')))
  })
  names(dmu_log)<-model$traits
  
  if(any(dmu_log>0)){
    Success=1
    error_trait<-names(dmu_log)[dmu_log>0]
  }
  
  if(all(dmu_log>0)){
    Success=2
  }
  
  
  ## 先检验SOL文件是否生成，再检验模型是否收敛
  
  if(!str_detect(method,"VAR$")){
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      lst<-paste0(outdir,"/DMU.lst")
      lst_con<-file(lst, open = "r")
      converged<-FALSE
      while (length(oneLine <- readLines(lst_con, n = 1, warn = FALSE)) > 0) {
        if(str_detect(oneLine, "Converged on")){
          converged<-TRUE
        }
      }
      close(lst_con)
      converge_msg<-""
      if(sol_exists & converged){
        converge_msg<-paste0("Variance components estimation for ",x," has successfully converged.\n")
      }
      if(sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & converged){
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & !converged){
        converge_msg<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg,SOL=sol_exists, CON=converged)
      return(res)
    })
  }
  
  if(Success!=2){
    ## 生成几个结果的表格
    ### Variance component
    
    if(!str_detect(method,"VAR$")){
      var_comp<-as.data.frame(matrix(NA,length(model$traits),length(rnd_code)+2))
      colnames(var_comp)<-c("Trait",int_name[rnd_code],"Residual")
      var_comp$Trait<-model$traits
      
      walk(model$traits,function(x){
        outdir<-toupper(x)
        if(converge_msg$SOL[converge_msg$TRT==x] & converge_msg$CON[converge_msg$TRT==x]){
          var_com<-read.table(paste0(outdir,"/PAROUT"),h=F)
          var_comp[var_comp$Trait==x,2:ncol(var_comp)]<<-var_com$V4
        }
      })
    }else{
      var_comp<-cbind(model$traits,var_comp_input)
      colnames(var_comp)<-c("Trait",int_name[rnd_code],"Residual")
    }
    write.table(var_comp,file="Variance_component.txt",row.names = F,quote=F,sep="\t")
    
    get_estimate<-function(trait){
      res<-list()
      res$fix_eff<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      res$genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      if(length(model$rnd_non)>0){res$non_genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)}
      outdir<-toupper(trait)
      if(converge_msg$SOL[converge_msg$TRT==trait] & converge_msg$CON[converge_msg$TRT==trait]){
        ### 读取solution
        sol<-fread(paste0(outdir,"/SOL"))
        fix_sol<-sol[sol$V1==2,]
        gen_sol<-sol[sol$V1==4,]
        nongen_sol<-sol[sol$V1==3,]
        if(method=="GBLUP"){
        gen_sol<-sol[sol$V1==3,]
        nongen_sol<-sol[sol$V1==4,]}
        if(length(model$fix_cov)>0){cov_sol<-sol[sol$V1==1,]}
        
        ### fixed factor effect
        fixed_fct_name<-c(model$fix_fct)
        
        fix_fct_sol<-map_dfr(fixed_fct_name,function(x){
          code<-fct_name[[x]]
          sol_dat<-fix_sol[fix_sol$V4==which(fixed_fct_name==x),]
          sol_dat2<-sol_dat[match(code$Code,sol_dat$V5),]
          sol_res<-data.frame(Name=x,Value=code$Name,Code=code$Code,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
          return(sol_res)
        })
        
        ### Covariate effect
        if(length(model$fix_cov)>0){
          cov_name<-model$fix_cov
          fix_cov_sol<-data.frame(Name=cov_name,Value=1,Code=1:length(cov_name),Estimate=cov_sol$V8,SE=cov_sol$V9)
        }else{
          fix_cov_sol<-NULL
        }
        
        ### Genetic effect
        if(method =="ssGBLUP"){
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fam_all[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name[[x]]$Name,Code=fct_name[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          non_genetic_sol<-non_genetic_sol[complete.cases(non_genetic_sol),]
          
        }
        
        
        
        res$fix_eff<-rbind(fix_fct_sol,fix_cov_sol)
        res$genetic_sol<-genetic_sol
        if(length(model$rnd_non)>0){res$non_genetic_sol<-non_genetic_sol}
      }
      return(res)
    }
    
    all_estimate<-map(model$traits,get_estimate)
    names(all_estimate)<-model$traits
    
    fix_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$fix_eff),idcol="Trait")
    fix_eff$Trait<-model$traits[fix_eff$Trait]
    
    gen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$genetic_sol),idcol="Trait")
    gen_eff$Trait<-model$traits[gen_eff$Trait]
    
    if(length(model$rnd_non)>0){
      nongen_eff<-rbindlist(map(model$traits,function(x)all_estimate[[x]]$non_genetic_sol),idcol="Trait")
      nongen_eff$Trait<-model$traits[nongen_eff$Trait]
    }
    
    fwrite(fix_eff, file="fixed_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    fwrite(gen_eff, file="genetic_effect_estimate.txt",row.names=F,quote=F,sep="\t")
    if(length(model$rnd_non)>0){fwrite(nongen_eff, file="nongenetic_random_effect_estimate.txt",row.names=F,quote=F,sep="\t")}
  }
}else{
  Success<-2
}


# 输出html格式报告

con<-file(paste0(report_path,"/report_template.R"), open = "r")
con_out<-file("report.R", open = "w")
model_formula<-paste(unlist(lapply(model$traits,function(x)paste0("\\\\text{",x,"}"))),
                     paste0("=",paste0("\\\\text{",paste(c(model$fix_fct,model$xp,model$rnd_gen,model$rnd_non),
                                                         collapse = "} + \\\\text{"),"}")),collapse = "\\\\\\\\")
model_formula<-paste0("The model is as follows:\n#' $$\n#' ",model_formula,"\n#' $$")

data_summary_entry<-"The basic information of your input data is as follows:\n#+ echo=FALSE\nknitr::kable(data_summary)"

if(length(model$rnd_non)>0){
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
#' 
#' 
#' 
#' ### Non-genetic random effects of first 10 records
#' 
#+ echo=FALSE
nongen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))

"
  
}else{
  model_res_entry<-"**Note that you can find all estimates and predictions in supplementary files.**
#'
#' ### Variance components
#+ echo=FALSE
knitr::kable(var_comp)
#' 
#' 
#' ### Fixed effect estimates
#' 
#+ echo=FALSE
DT::datatable(fix_eff,rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
              options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                                                             lengthChange = FALSE))
#'
#'
#'
#' ### Estimated breeding values of first 10 records
#' 
#+ echo=FALSE
gen_eff %>% 
  group_by(Trait) %>%
  do(head(.,10)) %>%
  DT::datatable(rownames=F,filter='top',selection=\"multiple\", escape=FALSE, 
                options = list(pageLength = 10,sDom  = '<\"top\">flrt<\"bottom\">ip',
                               lengthChange = FALSE))
"
}

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  
  if(model_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] The model must contain at least one fixed and one random effects. Please resubmit your job!\\*\\*", 
                         model_formula)
  }
  if(data_summary_check){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Please ensure that you have properly uploaded pedigree or genotype files!\\*\\*", 
                         data_summary_entry)
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no data check message is shown. Please check your data first!\\*\\*", 
                         paste0("After initial comparison, ",check_id_res$msg))
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, no phenotype data summary message is shown. Please check your data first!\\*\\*", 
                         phenotype_summary_msg)
  }
  
  if(data_summary_check & phenotype_check_pass){
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         genotype_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "genotype_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no genotype data summary message is shown. Please check your data first!**")
  }
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass){
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         pedigree_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "pedigree_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no pedigree data summary message is shown. Please check your data first!**")
  }
  
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass ){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/all.css ."))
knitr::spin("report.R",knit=F)
rmarkdown::render("report.Rmd")
system("rm report.Rmd")

# 打包附件
res_files<-map_chr(model$traits,function(x){
  system(paste0("cp ",toupper(x),"/DMU.lst", " ./",x,".dmu_run.lst"))
  return(paste0(x,".dmu_run.lst"))
})

if(Success!=2){
  if(length(model$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt","nongenetic_random_effect_estimate.txt")
  }else{
    res_files<-c(res_files,"genetic_effect_estimate.txt","fixed_effect_estimate.txt",
                 "Variance_component.txt")
  }
}

system(paste0("zip hegs_res.zip ", paste(res_files,collapse = " ")))

# 发送邮件
Success<-as.character(Success)
mail_text<-switch(Success,
                  "0"="Dear user,\n\nYour run of HEGS has sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "1"=paste0("Dear user,\n\nSome traits (",paste(error_trait,collapse = ", "),") of your run had not sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team"),
                  "2"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team")

email <- envelope(
  to = to_email,
  from = from_email,
  subject = paste0("HEGS report for job No.", task),
  text = mail_text
) %>% 
  attachment(path="report.html") %>%
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
}








