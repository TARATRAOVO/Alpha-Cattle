library(data.table)
library(pedigree)
library(tidyverse)
library(emayili)
library(stringr)
library(plyr)
library(dplyr)

options(warn=-1)
options <- commandArgs(trailingOnly = TRUE)
path = options[1]

#定义变量
#path<-"/disk195/zz/shinyApp/HEGS/demo/tb_demo"
setwd(path)

report_path<-"/srv/shiny-server"
load("par.RData")

#method<-"BLUP_VAR"
method<-res$method
model<-res$model
dir<-res$dir
phe_file<-res$phe_file
ped_file<-res$ped_file
gen_file<-res$gen_file
DMU1<-"/srv/shiny-server/dmu/dmu1"
DMUAI<-"/srv/shiny-server/dmu/dmuai"
DMU4<-"/srv/shiny-server/dmu/dmu4"
plink<-"/home/nyjh/miniconda2/pkgs/plink-1.90b6.21-hec16e2b_2/bin"
gcta<-"/home/nyjh/miniconda2/pkgs/gcta-1.93.2beta-h9ee0642_1/bin/gcta64"
from_email<-"852322127@qq.com"
to_email<-res$email
task<-res$task
type<-res$type

if(type=="Body_Type_Traits"){

# 定义一个确定是否需要基因型或者系谱数据的变量
need_genotype<-TRUE
if(method=="BLUP" | method=="BLUP_VAR"){
  need_genotype<-FALSE
}

need_pedigree<-TRUE
if(method=="GBLUP" | method=="GBLUP_VAR"){
  need_pedigree<-FALSE
}

# 定义一个是否需要方差组分文件的变量
need_var_comp<-TRUE
if(!str_detect(method,"VAR$")){
  need_var_comp<-FALSE
}

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
} ###这里如果model check不通过，是不是会没有报错信息。


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)
phe[[1]]=as.character(phe[[1]])

if(need_pedigree){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped[[1]]=as.character(ped[[1]])
  ped[[2]]=as.character(ped[[2]])
  ped[[3]]=as.character(ped[[3]])
}

if(need_genotype){
  fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)
  fam[[1]]=as.character(fam[[1]])
  fam[[2]]=as.character(fam[[2]])
}




if(need_var_comp){
  var_comp_input<-fread(res$var_file,h=F)
}

data_summary_check<-TRUE
if(model_check_pass & need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass & !need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped))
  )
}else if(model_check_pass & need_genotype & !need_pedigree){
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
  if(method=="BLUP" | method=="BLUP_VAR"){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    if(length(phe_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file.\n", length(phe_not_in_ped))
    }else{
      msg<-"all individuals with phenotype are in pedigree.\n"
    }
  }
  
  if(method=="GBLUP" | method=="GBLUP_VAR"){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" | method=="ssGBLUP_VAR"){
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
  if(method=="GBLUP" | method=="GBLUP_VAR"){
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
  
  if(nrow(phe)<=0 | sum(phe_summary_num)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass & need_genotype & phenotype_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_new)
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i]," --make-bed"))
      }
    }
    
    ### 过滤SNP MAF:0.01 HWE:1e-10
    num_of_variant<-rep(0,length(geno_names))
    for (i in 1:length(geno_names)){
      system(paste0(plink," --bfile ",
                    geno_names[i],
                    " --maf 0.01 --hwe 1e-10 --out ", geno_names[i], ".cln --make-bed"))
      
      con  <- file(paste0(geno_names[i],".cln.log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant.+people pass)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant.+people pass)"))
        }
    
      }
      
      close(con)
    }
    
    num_of_variant_final<-sum(num_of_variant)
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ", geno_names[i], ".cln --maf 0 --make-grm --autosome-num 29  --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg<-"[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg<-"[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}

## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & need_pedigree & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ord<-1:nrow(ped)
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999
    
    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F)
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}

## 检查方差组分的维度是否与性状需要的数量相对应
var_check_pass<-TRUE
var_summary_msg<-""

if(model_check_pass & need_var_comp & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
  if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)==(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of rows in variance component file is not equal to the number of traits, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)==length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of columns in variance component file is not equal to the number of random effects in the model, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] The dimension of variance component file is wrong, so the whole analysis is stopped."
  }
}

# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-c(model$fix_fct,model$rnd_non)
  fct<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe[[x]])))
  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })
  names(fct_name)<-fct_vec
  names(fct)<-fct_vec
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],ped[[1]])]
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],fam_new[[1]])]
    })
    rnd_gen_name<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<--9999
  fct[is.na(fct)]<--9999
  rl[is.na(rl)]<--9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_new){
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
  
  ### BLUP
  if(method=="BLUP"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

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
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
  # BLUP_VAR
  }else if(method=="BLUP_VAR"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

$PRIOR VAR_INPUT

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
    
    ### ssGBLUP_VAR
  }else if(method=="ssGBLUP_VAR"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
    
    ### GBLUP_VAR
  }else if(method=="GBLUP_VAR"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
  }else{
    
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
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
        if(need_pedigree){
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
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
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
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         var_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no variance component data summary message is shown. Please check your data first!**")
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/Rmd.css ."))
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


smtp <- server(host = "smtp.qq.com",
               port = 465,
               username = from_email,
               password = "cwdqmmxjoimhbbac")

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
    message("[WARNING:] Emailing the results caused a wanning!\n")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  }
)
}else if(type=="Milk_Production_Traits_inTDM"){

# 定义一个确定是否需要基因型或者系谱数据的变量
need_genotype<-TRUE
if(method=="BLUP" | method=="BLUP_VAR"){
  need_genotype<-FALSE
}

need_pedigree<-TRUE
if(method=="GBLUP" | method=="GBLUP_VAR"){
  need_pedigree<-FALSE
}

# 定义一个是否需要方差组分文件的变量
need_var_comp<-TRUE
if(!str_detect(method,"VAR$")){
  need_var_comp<-FALSE
}

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
} ###这里如果model check不通过，是不是会没有报错信息。


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)
phe=phe[phe[["DIM"]]>4 & phe[["DIM"]]<336]
phe[[1]]=as.character(phe[[1]])

if(need_pedigree){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped[[1]]=as.character(ped[[1]])
  ped[[2]]=as.character(ped[[2]])
  ped[[3]]=as.character(ped[[3]])
}

if(need_genotype){
  fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)
  fam[[1]]=as.character(fam[[1]])
  fam[[2]]=as.character(fam[[2]])
}




if(need_var_comp){
  var_comp_input<-fread(res$var_file,h=F)
}

data_summary_check<-TRUE
if(model_check_pass & need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass & !need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped))
  )
}else if(model_check_pass & need_genotype & !need_pedigree){
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
  if(method=="BLUP" | method=="BLUP_VAR"){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    if(length(phe_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file.\n", length(phe_not_in_ped))
    }else{
      msg<-"all individuals with phenotype are in pedigree.\n"
    }
  }
  
  if(method=="GBLUP" | method=="GBLUP_VAR"){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" | method=="ssGBLUP_VAR"){
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
  if(method=="GBLUP" | method=="GBLUP_VAR"){
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
  
  if(nrow(phe)<=0 | sum(phe_summary_num)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass & need_genotype & phenotype_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_new)
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i]," --make-bed"))
      }
    }
    
    ### 过滤SNP MAF:0.01 HWE:1e-10
    num_of_variant<-rep(0,length(geno_names))
    for (i in 1:length(geno_names)){
      system(paste0(plink," --bfile ",
                    geno_names[i],
                    " --maf 0.01 --hwe 1e-10 --out ", geno_names[i], ".cln --make-bed"))
      
      con  <- file(paste0(geno_names[i],".cln.log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant.+people pass)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant.+people pass)"))
        }
    
      }
      
      close(con)
    }
    
    num_of_variant_final<-sum(num_of_variant)
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ", geno_names[i], ".cln --maf 0 --make-grm --autosome-num 29  --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg<-"[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg<-"[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}

## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & need_pedigree & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ord<-1:nrow(ped)
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999
    
    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F)
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}

## 检查方差组分的维度是否与性状需要的数量相对应
var_check_pass<-TRUE
var_summary_msg<-""

if(model_check_pass & need_var_comp & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
  if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)==(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of rows in variance component file is not equal to the number of traits, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)==length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of columns in variance component file is not equal to the number of random effects in the model, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] The dimension of variance component file is wrong, so the whole analysis is stopped."
  }
}

# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-model$fix_fct
  fct<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe[[x]])))
  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })
  names(fct_name)<-fct_vec
  names(fct)<-fct_vec
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],ped[[1]])]
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],fam_new[[1]])]
    })
    rnd_gen_name<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
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
  lg_seg=lg[,4:8]
  
  lg_DIMph=data.frame(lg[,3])
  names(lg_DIMph)="DIMph"
  fct_prenames=names(fct)
  fct<-cbind(fct,lg_DIMph)
  names(fct)=c(fct_prenames,"DIMph")

  rl_lg_p1=select(rl,model$fix_cov)
  rl_lg_p2=lg_seg
  rl_lg_p3=select(rl,model$traits)
  rl_lg=cbind(rl_lg_p1,rl_lg_p2,rl_lg_p3)

  rnd_gen[is.na(rnd_gen)]<--9999
  fct[is.na(fct)]<--9999
  rl[is.na(rl)]<--9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_name1<-c(model$rnd_gen,model$fix_fct,"DIMph")
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  int_code_p1<-c(match(model$fix_fct,int_name))
  int_code_p2<-paste(paste(match(model$rnd_gen,int_name),"(","0",")",sep=""),paste(match(model$rnd_gen,int_name),"(","0",")",sep=""),collapse=" ")
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_new){
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
  
  ### BLUP
  if(method=="BLUP"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### GBLUP
  }else if(method=="GBLUP"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
  # BLUP_VAR
  }else if(method=="BLUP_VAR"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

$PRIOR VAR_INPUT

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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### ssGBLUP_VAR
  }else if(method=="ssGBLUP_VAR"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
                  length(rnd_code),paste(1:length(rnd_code),collapse = " "),
                  length(model$fix_cov)+15, ifelse(length(model$fix_cov)==0,paste(c(paste(c(1,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "),
    paste(c(paste(1:length(model$fix_cov),collapse=" "),paste(c(2,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(3,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(4,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(5,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = ""),paste(c(6,"(",paste(c("0",which(int_code==1)),collapse =" "),")"),collapse = "")),collapse = " "))),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### GBLUP_VAR
  }else if(method=="GBLUP_VAR"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,"lg0","lg1","lg2","lg3","lg4",x)
      dmu_dat<-cbind(rnd_gen,fct,rl_lg)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
                  paste(c(int_code_p1,int_code_p2),collapse = " "),
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
  }else{
    
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
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
        if(need_pedigree){
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
            sol_dat3<-aggregate(V8 ~ V5, data = sol_dat2, meantdm)
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,Estimate=sol_dat3$V8)
            return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model$rnd_gen==x),]
            sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
            sol_dat3<-aggregate(V8 ~ V5, data = sol_dat2, meantdm)
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat3$V8)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          fct_vec1<-c(model$fix_fct,model$rnd_non)
          fct1<-map_dfc(fct_vec1,function(x)as.numeric(as.factor(phe[[x]])))
          fct_name1<-map(fct_vec1,function(x){
            name<-levels(as.factor(phe[[x]]))
            data.frame(Name=name,Code=1:length(name))
            })
          names(fct_name1)<-fct_vec1
          non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat%>%arrange(match(V5,ord), desc(V1),desc(V2),desc(V3),desc(V4),desc(V6),desc(V7),desc(V8),desc(V9))
            sol_dat3<-aggregate(V8 ~ V7, data = sol_dat2, meantdm)
            if(length(fct_name1[[x]]$Code)==length(sol_dat3$V8)){
            sol_res<-data.frame(Name=x,Value=fct_name1[[x]]$Name,Code=fct_name1[[x]]$Code,
                                Estimate=sol_dat3$V8)
            }else if(length(fct_name1[[x]]$Code)>length(sol_dat3$V8)){
            non_gen_estimate=sol_dat3$V8
            non_gen_estimate[setdiff(fct_name1[[x]]$Code,sol_dat3$V7)]=NA
            sol_res<-data.frame(Name=x,Value=fct_name1[[x]]$Name,Code=fct_name1[[x]]$Code,
                                Estimate=non_gen_estimate)
            }else if(length(fct_name1[[x]]$Code)<length(sol_dat3$V8)){
            non_gen_estimate=sol_dat3$V8
            non_gen_estimate=non_gen_estimate[intersect(sol_dat3$V7,fct_name1[[x]]$Code)]
            sol_res<-data.frame(Name=x,Value=fct_name1[[x]]$Name,Code=fct_name1[[x]]$Code,
                                Estimate=non_gen_estimate)
            }
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
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         var_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no variance component data summary message is shown. Please check your data first!**")
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/Rmd.css ."))
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


smtp <- server(host = "smtp.qq.com",
               port = 465,
               username = from_email,
               password = "cwdqmmxjoimhbbac")

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
    message("[WARNING:] Emailing the results caused a wanning!\n")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  }
)
}else if(type=="Reproductive_Traits"){
# 定义一个确定是否需要基因型或者系谱数据的变量
need_genotype<-TRUE
if(method=="BLUP" | method=="BLUP_VAR"){
  need_genotype<-FALSE
}

need_pedigree<-TRUE
if(method=="GBLUP" | method=="GBLUP_VAR"){
  need_pedigree<-FALSE
}

# 定义一个是否需要方差组分文件的变量
need_var_comp<-TRUE
if(!str_detect(method,"VAR$")){
  need_var_comp<-FALSE
}

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
} ###这里如果model check不通过，是不是会没有报错信息。


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)
phe[[1]]=as.character(phe[[1]])

if(need_pedigree){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped[[1]]=as.character(ped[[1]])
  ped[[2]]=as.character(ped[[2]])
  ped[[3]]=as.character(ped[[3]])
}

if(need_genotype){
  fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)
  fam[[1]]=as.character(fam[[1]])
  fam[[2]]=as.character(fam[[2]])
}




if(need_var_comp){
  var_comp_input<-fread(res$var_file,h=F)
}

data_summary_check<-TRUE
if(model_check_pass & need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass & !need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped))
  )
}else if(model_check_pass & need_genotype & !need_pedigree){
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
  if(method=="BLUP" | method=="BLUP_VAR"){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    if(length(phe_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file.\n", length(phe_not_in_ped))
    }else{
      msg<-"all individuals with phenotype are in pedigree.\n"
    }
  }
  
  if(method=="GBLUP" | method=="GBLUP_VAR"){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" | method=="ssGBLUP_VAR"){
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
  if(method=="GBLUP" | method=="GBLUP_VAR"){
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
  
  if(nrow(phe)<=0 | sum(phe_summary_num)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass & need_genotype & phenotype_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_new)
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i]," --make-bed"))
      }
    }
    
    ### 过滤SNP MAF:0.01 HWE:1e-10
    num_of_variant<-rep(0,length(geno_names))
    for (i in 1:length(geno_names)){
      system(paste0(plink," --bfile ",
                    geno_names[i],
                    " --maf 0.01 --hwe 1e-10 --out ", geno_names[i], ".cln --make-bed"))
      
      con  <- file(paste0(geno_names[i],".cln.log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant.+people pass)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant.+people pass)"))
        }
    
      }
      
      close(con)
    }
    
    num_of_variant_final<-sum(num_of_variant)
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ", geno_names[i], ".cln --maf 0 --make-grm --autosome-num 29  --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg<-"[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg<-"[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}

## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & need_pedigree & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ord<-1:nrow(ped)
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999
    
    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F)
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}

## 检查方差组分的维度是否与性状需要的数量相对应
var_check_pass<-TRUE
var_summary_msg<-""

if(model_check_pass & need_var_comp & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
  if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)==(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of rows in variance component file is not equal to the number of traits, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)==length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of columns in variance component file is not equal to the number of random effects in the model, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] The dimension of variance component file is wrong, so the whole analysis is stopped."
  }
}

# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-c(model$fix_fct,model$rnd_non)
  fct<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe[[x]])))
  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })
  names(fct_name)<-fct_vec
  names(fct)<-fct_vec
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe,c(model$fix_cov,model$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],ped[[1]])]
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe[[x]],fam_new[[1]])]
    })
    rnd_gen_name<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<--9999
  fct[is.na(fct)]<--9999
  rl[is.na(rl)]<--9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_new){
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
  
  ### BLUP
  if(method=="BLUP"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

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
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
  # BLUP_VAR
  }else if(method=="BLUP_VAR"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

$PRIOR VAR_INPUT

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
    
    ### ssGBLUP_VAR
  }else if(method=="ssGBLUP_VAR"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
    
    ### GBLUP_VAR
  }else if(method=="GBLUP_VAR"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
  }else{
    
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
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
        if(need_pedigree){
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
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          if(model$rnd_non=="ID"){
           non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          }else{
            non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name[[x]]$Name,Code=fct_name[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          }
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
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         var_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no variance component data summary message is shown. Please check your data first!**")
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/Rmd.css ."))
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


smtp <- server(host = "smtp.qq.com",
               port = 465,
               username = from_email,
               password = "cwdqmmxjoimhbbac")

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
    message("[WARNING:] Emailing the results caused a wanning!\n")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  }
)
}else if(type=="Milk_Production_Traits_inAVE"){
need_genotype<-TRUE
if(method=="BLUP" | method=="BLUP_VAR"){
  need_genotype<-FALSE
}

need_pedigree<-TRUE
if(method=="GBLUP" | method=="GBLUP_VAR"){
  need_pedigree<-FALSE
}

# 定义一个是否需要方差组分文件的变量
need_var_comp<-TRUE
if(!str_detect(method,"VAR$")){
  need_var_comp<-FALSE
}

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model$fix_fct)==0 | length(model$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model$rnd_gen<-model$rnd_gen[1]
} ###这里如果model check不通过，是不是会没有报错信息。


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)
phe[[1]]=as.character(phe[[1]])

if(need_pedigree){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
  ped[[1]]=as.character(ped[[1]])
  ped[[2]]=as.character(ped[[2]])
  ped[[3]]=as.character(ped[[3]])
}

if(need_genotype){
  fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)
  fam[[1]]=as.character(fam[[1]])
  fam[[2]]=as.character(fam[[2]])
}




if(need_var_comp){
  var_comp_input<-fread(res$var_file,h=F)
}

data_summary_check<-TRUE
if(model_check_pass & need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree", "Number of SNPs", "Number of genotyped individuals"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped), res$num_of_snps,
           nrow(fam))
  )
}else if(model_check_pass & !need_genotype & need_pedigree){
  data_summary<-data.frame(
    Item=c("Number of phenotype records"," number of phenotyped individuals",
           "Number of animals in pedigree"),
    Size=c(nrow(phe), length(unique(phe[[1]])),
           nrow(ped))
  )
}else if(model_check_pass & need_genotype & !need_pedigree){
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
  if(method=="BLUP" | method=="BLUP_VAR"){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    if(length(phe_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file.\n", length(phe_not_in_ped))
    }else{
      msg<-"all individuals with phenotype are in pedigree.\n"
    }
  }
  
  if(method=="GBLUP" | method=="GBLUP_VAR"){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
    }
  }
  
  if(method=="ssGBLUP" | method=="ssGBLUP_VAR"){
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
  if(method=="GBLUP" | method=="GBLUP_VAR"){
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
  
  if(nrow(phe)<=0 | sum(phe_summary_num)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
}

## 对基因型的处理
genotype_check_pass<-TRUE
genotype_summary_msg<-""
if(model_check_pass & need_genotype & phenotype_check_pass){
  ### 去除不用的基因型个体
  geno_names<-str_extract(res$gen_file,".+(?=\\.bed$)")
  geno_names<-geno_names[!is.na(geno_names)]
  fam_new<-fam[!fam$V1%in%check_id_res$fam_not_in_ped,]
  
  if(nrow(fam_new)>0){
    ord<-1:nrow(fam_new)
    if(length(check_id_res$fam_not_in_ped)>0){
      for (i in 1:length(geno_names)){
        system(paste0(plink," --bfile ",
                      geno_names[i]," --remove-fam fam_id_not_in_ped.txt --out ", 
                      geno_names[i]," --make-bed"))
      }
    }
    
    ### 过滤SNP MAF:0.01 HWE:1e-10
    num_of_variant<-rep(0,length(geno_names))
    for (i in 1:length(geno_names)){
      system(paste0(plink," --bfile ",
                    geno_names[i],
                    " --maf 0.01 --hwe 1e-10 --out ", geno_names[i], ".cln --make-bed"))
      
      con  <- file(paste0(geno_names[i],".cln.log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant.+people pass)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant.+people pass)"))
        }
    
      }
      
      close(con)
    }
    
    num_of_variant_final<-sum(num_of_variant)
    
    if(num_of_variant_final>0){
      genotype_summary_msg<-paste0(num_of_variant_final," SNPS in ", nrow(fam_new),
                                   " individuals were used in further analysis")
      
      ## 计算GRM
      if(length(geno_names)>=1){
        for (i in 1:length(geno_names)){
          system(paste0(gcta, " --bfile ", geno_names[i], ".cln --maf 0 --make-grm --autosome-num 29  --out ", geno_names[i], " --make-grm-alg 1"))
        }
      }
      write.table(geno_names,file="grm.lst",row.names = F,col.names = F,quote = F)
      system(paste0(gcta, " --mgrm grm.lst --out final --make-grm-gz"))
    }else{
      genotype_check_pass<-FALSE
      genotype_summary_msg<-"[ERROR:] no SNP is left in the data, so the whole analysis is stopped."
    }
    
  }else{
    genotype_summary_msg<-"[ERROR:] no individual is left in the data, so the whole analysis is stopped."
    genotype_check_pass<-FALSE
  }
}

## 将系谱按照从小到大编码并产生DMU格式的系谱
pedigree_check_pass<-TRUE
pedigree_summary_msg<-""
if(model_check_pass & need_pedigree & phenotype_check_pass & genotype_check_pass){
  if(nrow(ped)>0){
    ord<-1:nrow(ped)
    ped_ind<-ord
    sire_ind<-ord[match(ped[[2]],ped[[1]])]
    dam_ind<-ord[match(ped[[3]],ped[[1]])]
    dmu_ped<-data.frame(ped_ind,sire_ind,dam_ind,ord)
    dmu_ped[is.na(dmu_ped)]<--9999
    
    write.table(dmu_ped,file="DMU_PED",row.names=F,col.names=F,quote=F)
  }else{
    pedigree_check_pass<-FALSE
    pedigree_summary_msg<-"[ERROR:] No individual is left in the pedigree, so the whole analysis is stopped."
  }
  
}

## 检查方差组分的维度是否与性状需要的数量相对应
var_check_pass<-TRUE
var_summary_msg<-""

if(model_check_pass & need_var_comp & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
  if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)==(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of rows in variance component file is not equal to the number of traits, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)==length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] Number of columns in variance component file is not equal to the number of random effects in the model, so the whole analysis is stopped."
  }else if(nrow(var_comp_input)!=length(model$traits) & ncol(var_comp_input)!=(length(model$rnd_gen)+length(model$rnd_non)+1)){
    var_check_pass<-FALSE
    var_summary_msg<-"[ERROR:] The dimension of variance component file is wrong, so the whole analysis is stopped."
  }
}

# 上述所有数据检验通过之后，开始进入数据分析流程
Success<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
  ##对表型文件的处理
  names(phe)[str_detect(names(phe),fixed("parity",ignore_case = T))]="Parity"
  names(phe)[str_detect(names(phe),fixed("dim",ignore_case = T))]="DIM"
  phe=phe[phe[["DIM"]]>4 & phe[["DIM"]]<336]
  if(!is.null(model$rnd_non)){
     walk(model$traits,function(x){
        if(str_detect(x,fixed("milkyield",ignore_case = T))){
           phe1=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],x))}
           phe1=phe1[order(phe1[[1]],phe1$Parity,phe1$DIM),]
           phe1$ID_Parity=paste(setDF(phe1)[[1]],setDF(phe1)$Parity,sep="")
           ID_Parity_all=unique(phe1$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("mydf",i),phe1[phe1$ID_Parity==ID_Parity_all[i],])
               a2=sum(get(paste0("mydf",i))$DIM<305)
               if(get(paste0("mydf",i))$Parity[1]==1){
                   milk1=get(paste0("mydf",i))$DIM[1]*get(paste0("mydf",i))[[x]][1]*(0.605+0.0435*sqrt(get(paste0("mydf",i))$DIM[1]))         ##第一个矫正时间段
                 if(get(paste0("mydf",i))$DIM[1]>40){
                             for(j in 2:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                 }else if(get(paste0("mydf",i))$DIM[1]<40){
                       if(get(paste0("mydf",i))$DIM[2]>40){
                             milk2=0.5*(0.998+1.23/((get(paste0("mydf",i))$DIM[1])^2)+0.0113*((get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])+1)/(get(paste0("mydf",i))$DIM[1]))*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                       }else if(get(paste0("mydf",i))$DIM[2]<40){
                             milk2=0.5*(0.998+1.23/((get(paste0("mydf",i))$DIM[2])^2)+0.0113*((get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])+1)/(get(paste0("mydf",i))$DIM[2]))*(get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])*(get(paste0("mydf",i))[[x]][3]+get(paste0("mydf",i))[[x]][2])+0.5*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:a2-1){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j+1]-get(paste0("mydf",i))$DIM[j])*(get(paste0("mydf",i))[[x]][j+1]+get(paste0("mydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2),get))))
                       }
                   }
                 }
               }else if(get(paste0("mydf",i))$Parity[1]>1){
                   milk1=get(paste0("mydf",i))$DIM[1]*get(paste0("mydf",i))[[x]][1]*(0.635+0.0435*sqrt(get(paste0("mydf",i))$DIM[1]))         ##第一个矫正时间段
                 if(get(paste0("mydf",i))$DIM[1]>40){
                             for(j in 2:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                 }else if(get(paste0("mydf",i))$DIM[1]<40){
                       if(get(paste0("mydf",i))$DIM[2]>40){
                             milk2=0.5*(1.001-0.00042*(get(paste0("mydf",i))$DIM[1])+0.0109*((get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])+1)/(get(paste0("mydf",i))$DIM[1]))*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                 assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                       }else if(get(paste0("mydf",i))$DIM[2]<40){
                             milk2=0.5*(1.001-0.00042*(get(paste0("mydf",i))$DIM[2])+0.0109*((get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])+1)/(get(paste0("mydf",i))$DIM[2]))*(get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])*(get(paste0("mydf",i))[[x]][3]+get(paste0("mydf",i))[[x]][2])+0.5*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:a2-1){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j+1]-get(paste0("mydf",i))$DIM[j])*(get(paste0("mydf",i))[[x]][j+1]+get(paste0("mydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                 assign(paste0("milk",a2),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("mydf",i))[[model$rnd_gen]][a2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("mydf",i))[[model$fix_fct[k]]][a2])}
               assign(paste0("Parity",i),get(paste0("mydf",i))[["Parity"]][a2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("mydf",i))[[model$fix_cov[k]]][a2])}
               if(!(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non))){
                     a3=model$rnd_non[!(model$rnd_non=="ID")]
                     a4=a3[!(a3=="id")]
                     for(k in 1:(length(a4))){
                          assign(paste0(a4[k],i),get(paste0("mydf",i))[[a4[k]]][a2])}
               }
           }
           phe1_my305=data.table(unlist(lapply(paste0("my_305all",1:length(ID_Parity_all)),get)))
           phe1_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe1_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe1_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe1_fixfct=cbind(phe1_fixfct,phe1_fixfct1)}
           phe1_fixfct=phe1_fixfct[,-1]
           phe1_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe1_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe1_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe1_fixcov=cbind(phe1_fixcov,phe1_fixcov1)}
           phe1_fixcov=phe1_fixcov[,-1]
           if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
               phe11=cbind(phe1_rndgen,phe1_fixfct,phe1_Parity,phe1_fixcov,phe1_my305)
           }else{
               a3=model$rnd_non[!(model$rnd_non=="ID")]
               a4=a3[!(a3=="id")]
               phe1_rndnon=data.table(rep(NA,length(ID_Parity_all)))
               for(k in 1:(length(a4))){
                    phe1_rndnon1=data.table(unlist(lapply(paste0(a4[i],1:length(ID_Parity_all)),get)))
                    phe1_rndnon=cbind(phe1_rndnon,phe1_rndnon1)}
               phe1_rndnon=phe1_rndnon[,-1]
               phe11=cbind(phe1_rndgen,phe1_fixfct,phe1_Parity,phe1_fixcov,phe1_rndnon,phe1_my305)
           }
        }
        if(str_detect(x,fixed("fatyield",ignore_case = T))){
           phe2=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],x))}
           phe2=phe2[order(phe2[[1]],phe2$Parity,phe2$DIM),]
           phe2$ID_Parity=paste(setDF(phe2)[[1]],setDF(phe2)$Parity,sep="")
           ID_Parity_all=unique(phe2$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("fydf",i),phe2[phe2$ID_Parity==ID_Parity_all[i],])
               b2=sum(get(paste0("fydf",i))$DIM<305)
               if(get(paste0("fydf",i))$Parity[1]==1){
                   fat1=get(paste0("fydf",i))$DIM[1]*get(paste0("fydf",i))[[x]][1]*(0.235+0.239*sqrt(get(paste0("fydf",i))$DIM[1])-0.0225*get(paste0("fydf",i))$DIM[1]+0.000069*(get(paste0("fydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("fydf",i))$DIM[1]>40){
                             for(j in 2:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                 }else if(get(paste0("fydf",i))$DIM[1]<40){
                       if(get(paste0("fydf",i))$DIM[2]>40){
                             fat2=0.5*(0.998+1.23/((get(paste0("fydf",i))$DIM[1])^2)+0.0113*((get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])+1)/(get(paste0("fydf",i))$DIM[1]))*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                       }else if(get(paste0("fydf",i))$DIM[2]<40){
                             fat2=0.5*(0.998+1.23/((get(paste0("fydf",i))$DIM[2])^2)+0.0113*((get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])+1)/(get(paste0("fydf",i))$DIM[2]))*(get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])*(get(paste0("fydf",i))[[x]][3]+get(paste0("fydf",i))[[x]][2])+0.5*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:b2-1){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j+1]-get(paste0("fydf",i))$DIM[j])*(get(paste0("fydf",i))[[x]][j+1]+get(paste0("fydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2),get))))
                       }
                   }
                 }
               }else if(get(paste0("fydf",i))$Parity[1]>1){
                   fat1=get(paste0("fydf",i))$DIM[1]*get(paste0("fydf",i))[[x]][1]*(0.476+0.146*sqrt(get(paste0("fydf",i))$DIM[1])-0.0115*get(paste0("fydf",i))$DIM[1]+0.000038*(get(paste0("fydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("fydf",i))$DIM[1]>40){
                             for(j in 2:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                 }else if(get(paste0("fydf",i))$DIM[1]<40){
                       if(get(paste0("fydf",i))$DIM[2]>40){
                             fat2=0.5*(1.001-0.00042*(get(paste0("fydf",i))$DIM[1])+0.0109*((get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])+1)/(get(paste0("fydf",i))$DIM[1]))*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                 assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                       }else if(get(paste0("fydf",i))$DIM[2]<40){
                             fat2=0.5*(1.001-0.00042*(get(paste0("fydf",i))$DIM[2])+0.0109*((get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])+1)/(get(paste0("fydf",i))$DIM[2]))*(get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])*(get(paste0("fydf",i))[[x]][3]+get(paste0("fydf",i))[[x]][2])+0.5*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:b2-1){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j+1]-get(paste0("fydf",i))$DIM[j])*(get(paste0("fydf",i))[[x]][j+1]+get(paste0("fydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                 assign(paste0("fat",b2),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("fydf",i))[[model$rnd_gen]][b2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("fydf",i))[[model$fix_fct[k]]][b2])}
               assign(paste0("Parity",i),get(paste0("fydf",i))[["Parity"]][b2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("fydf",i))[[model$fix_cov[k]]][b2])}
               if(!(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non))){
                     b3=model$rnd_non[!(model$rnd_non=="ID")]
                     b4=b3[!(b3=="id")]
                     for(k in 1:(length(b4))){
                          assign(paste0(b4[k],i),get(paste0("fydf",i))[[b4[k]]][b2])}
               }
           }
           phe2_fy305=data.table(unlist(lapply(paste0("fy_305all",1:length(ID_Parity_all)),get)))
           phe2_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe2_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe2_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe2_fixfct=cbind(phe2_fixfct,phe2_fixfct1)}
           phe2_fixfct=phe2_fixfct[,-1]
           phe2_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe2_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe2_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe2_fixcov=cbind(phe2_fixcov,phe2_fixcov1)}
           phe2_fixcov=phe2_fixcov[,-1]
           if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
               phe22=cbind(phe2_rndgen,phe2_fixfct,phe2_Parity,phe2_fixcov,phe2_fy305)
           }else{
               b3=model$rnd_non[!(model$rnd_non=="ID")]
               b4=b3[!(b3=="id")]
               phe2_rndnon=data.table(rep(NA,length(ID_Parity_all)))
               for(k in 1:(length(b4))){
                    phe2_rndnon1=data.table(unlist(lapply(paste0(b4[i],1:length(ID_Parity_all)),get)))
                    phe2_rndnon=cbind(phe2_rndnon,phe2_rndnon1)}
               phe2_rndnon=phe2_rndnon[,-1]
               phe22=cbind(phe2_rndgen,phe2_fixfct,phe2_Parity,phe2_fixcov,phe2_rndnon,phe2_fy305)
           }
        }
        if(str_detect(x,fixed("proyield",ignore_case = T))){
           phe3=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],x))}
           phe3=phe3[order(phe3[[1]],phe3$Parity,phe3$DIM),]
           phe3$ID_Parity=paste(setDF(phe3)[[1]],setDF(phe3)$Parity,sep="")
           ID_Parity_all=unique(phe3$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("pydf",i),phe3[phe3$ID_Parity==ID_Parity_all[i],])
               c2=sum(get(paste0("pydf",i))$DIM<305)
               if(get(paste0("pydf",i))$Parity[1]==1){
                   pro1=get(paste0("pydf",i))$DIM[1]*get(paste0("pydf",i))[[x]][1]*(0.136+0.316*sqrt(get(paste0("pydf",i))$DIM[1])-0.0351*get(paste0("pydf",i))$DIM[1]+0.00013*(get(paste0("pydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("pydf",i))$DIM[1]>40){
                             for(j in 2:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                 }else if(get(paste0("pydf",i))$DIM[1]<40){
                       if(get(paste0("pydf",i))$DIM[2]>40){
                             pro2=0.5*(0.998+1.23/((get(paste0("pydf",i))$DIM[1])^2)+0.0113*((get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])+1)/(get(paste0("pydf",i))$DIM[1]))*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                       }else if(get(paste0("pydf",i))$DIM[2]<40){
                             pro2=0.5*(0.998+1.23/((get(paste0("pydf",i))$DIM[2])^2)+0.0113*((get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])+1)/(get(paste0("pydf",i))$DIM[2]))*(get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])*(get(paste0("pydf",i))[[x]][3]+get(paste0("pydf",i))[[x]][2])+0.5*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:c2-1){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j+1]-get(paste0("pydf",i))$DIM[j])*(get(paste0("pydf",i))[[x]][j+1]+get(paste0("pydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2),get))))
                       }
                   }
                 }
               }else if(get(paste0("pydf",i))$Parity[1]>1){
                   pro1=get(paste0("pydf",i))$DIM[1]*get(paste0("pydf",i))[[x]][1]*(0.177+0.324*sqrt(get(paste0("pydf",i))$DIM[1])-0.0366*get(paste0("pydf",i))$DIM[1]+0.000141*(get(paste0("pydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("pydf",i))$DIM[1]>40){
                             for(j in 2:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                 }else if(get(paste0("pydf",i))$DIM[1]<40){
                       if(get(paste0("pydf",i))$DIM[2]>40){
                             pro2=0.5*(1.001-0.00042*(get(paste0("pydf",i))$DIM[1])+0.0109*((get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])+1)/(get(paste0("pydf",i))$DIM[1]))*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                 assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                       }else if(get(paste0("pydf",i))$DIM[2]<40){
                             pro2=0.5*(1.001-0.00042*(get(paste0("pydf",i))$DIM[2])+0.0109*((get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])+1)/(get(paste0("pydf",i))$DIM[2]))*(get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])*(get(paste0("pydf",i))[[x]][3]+get(paste0("pydf",i))[[x]][2])+0.5*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:c2-1){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j+1]-get(paste0("pydf",i))$DIM[j])*(get(paste0("pydf",i))[[x]][j+1]+get(paste0("pydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                 assign(paste0("pro",c2),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("pydf",i))[[model$rnd_gen]][c2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("pydf",i))[[model$fix_fct[k]]][c2])}
               assign(paste0("Parity",i),get(paste0("pydf",i))[["Parity"]][c2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("pydf",i))[[model$fix_cov[k]]][c2])}
               if(!(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non))){
                     c3=model$rnd_non[!(model$rnd_non=="ID")]
                     c4=c3[!(c3=="id")]
                     for(k in 1:(length(c4))){
                          assign(paste0(c4[k],i),get(paste0("pydf",i))[[c4[k]]][c2])}
               }
           }
           phe3_py305=data.table(unlist(lapply(paste0("py_305all",1:length(ID_Parity_all)),get)))
           phe3_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe3_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe3_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe3_fixfct=cbind(phe3_fixfct,phe3_fixfct1)}
           phe3_fixfct=phe3_fixfct[,-1]
           phe3_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe3_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe3_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe3_fixcov=cbind(phe3_fixcov,phe3_fixcov1)}
           phe3_fixcov=phe3_fixcov[,-1]
           if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
               phe33=cbind(phe3_rndgen,phe3_fixfct,phe3_Parity,phe3_fixcov,phe3_py305)
           }else{
               c3=model$rnd_non[!(model$rnd_non=="ID")]
               c4=c3[!(c3=="id")]
               phe3_rndnon=data.table(rep(NA,length(ID_Parity_all)))
               for(k in 1:(length(c4))){
                    phe3_rndnon1=data.table(unlist(lapply(paste0(c4[i],1:length(ID_Parity_all)),get)))
                    phe3_rndnon=cbind(phe3_rndnon,phe3_rndnon1)}
               phe3_rndnon=phe3_rndnon[,-1]
               phe33=cbind(phe3_rndgen,phe3_fixfct,phe3_Parity,phe3_fixcov,phe3_rndnon,phe3_py305)
           }
        }
        if(str_detect(x,fixed("SCS",ignore_case = T))){
           phe4=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],x))}
           phe4=phe4[order(phe4[[1]],phe4$Parity,phe4$DIM),]
           phe4$ID_Parity=paste(setDF(phe4)[[1]],setDF(phe4)$Parity,sep="")
           ID_Parity_all=unique(phe4$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("scdf",i),phe4[phe4$ID_Parity==ID_Parity_all[i],])
               d2=sum(get(paste0("scdf",i))$DIM<305)
               assign(paste0("scc",i),2^(get(paste0("scdf",i))[[x]]-3)*100)
               assign(paste0("scc_305",i),mean(get(paste0("scc",i))))              ##SCS与SCC直接取平均值
               assign(paste0("scs_305",i),log2(get(paste0("scc_305",i))/100)+3)
               assign(paste0(model$rnd_gen,i),get(paste0("scdf",i))[[model$rnd_gen]][d2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("scdf",i))[[model$fix_fct[k]]][d2])}
               assign(paste0("Parity",i),get(paste0("scdf",i))[["Parity"]][d2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("scdf",i))[[model$fix_cov[k]]][d2])}
               if(!(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non))){
                     d3=model$rnd_non[!(model$rnd_non=="ID")]
                     d4=d3[!(d3=="id")]
                     for(k in 1:(length(d4))){
                          assign(paste0(d4[k],i),get(paste0("scdf",i))[[d4[k]]][d2])}
               }
           }
           phe4_scc305=data.table(unlist(lapply(paste0("scc_305",1:length(ID_Parity_all)),get)))
           phe4_scs305=data.table(unlist(lapply(paste0("scs_305",1:length(ID_Parity_all)),get)))
           phe4_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe4_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe4_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe4_fixfct=cbind(phe4_fixfct,phe4_fixfct1)}
           phe4_fixfct=phe4_fixfct[,-1]
           phe4_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe4_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe4_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe4_fixcov=cbind(phe4_fixcov,phe4_fixcov1)}
           phe4_fixcov=phe4_fixcov[,-1]
           if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
               phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_scc305,phe4_scs305)
           }else{
               d3=model$rnd_non[!(model$rnd_non=="ID")]
               d4=d3[!(d3=="id")]
               phe4_rndnon=data.table(rep(NA,length(ID_Parity_all)))
               for(k in 1:(length(d4))){
                    phe4_rndnon1=data.table(unlist(lapply(paste0(d4[i],1:length(ID_Parity_all)),get)))
                    phe4_rndnon=cbind(phe4_rndnon,phe4_rndnon1)}
               phe4_rndnon=phe4_rndnon[,-1]
               phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_rndnon,phe4_scc305,phe4_scs305)
           }
        }
        if(str_detect(x,fixed("SCC",ignore_case = T))){
           phe4=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],x))}
           phe4=phe4[order(phe4[[1]],phe4$Parity,phe4$DIM),]
           phe4$ID_Parity=paste(setDF(phe4)[[1]],setDF(phe4)$Parity,sep="")
           ID_Parity_all=unique(phe4$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("scdf",i),phe4[phe4$ID_Parity==ID_Parity_all[i],])
               d2=sum(get(paste0("scdf",i))$DIM<305)
               assign(paste0("scc_305",i),mean(get(paste0("scdf",i))[[x]]))              ##SCS与SCC直接取平均值
               assign(paste0("scs_305",i),log2(get(paste0("scc_305",i))/100)+3)
               assign(paste0(model$rnd_gen,i),get(paste0("scdf",i))[[model$rnd_gen]][d2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("scdf",i))[[model$fix_fct[k]]][d2])}
               assign(paste0("Parity",i),get(paste0("scdf",i))[["Parity"]][d2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("scdf",i))[[model$fix_cov[k]]][d2])}
               if(!(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non))){
                     d3=model$rnd_non[!(model$rnd_non=="ID")]
                     d4=d3[!(d3=="id")]
                     for(k in 1:(length(d4))){
                          assign(paste0(d4[k],i),get(paste0("scdf",i))[[d4[k]]][d2])}
               }
           }
           phe4_scc305=data.table(unlist(lapply(paste0("scc_305",1:length(ID_Parity_all)),get)))
           phe4_scs305=data.table(unlist(lapply(paste0("scs_305",1:length(ID_Parity_all)),get)))
           phe4_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe4_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe4_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe4_fixfct=cbind(phe4_fixfct,phe4_fixfct1)}
           phe4_fixfct=phe4_fixfct[,-1]
           phe4_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe4_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe4_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe4_fixcov=cbind(phe4_fixcov,phe4_fixcov1)}
           phe4_fixcov=phe4_fixcov[,-1]
           if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
               phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_scc305,phe4_scs305)
           }else{
               d3=model$rnd_non[!(model$rnd_non=="ID")]
               d4=d3[!(d3=="id")]
               phe4_rndnon=data.table(rep(NA,length(ID_Parity_all)))
               for(k in 1:(length(d4))){
                    phe4_rndnon1=data.table(unlist(lapply(paste0(d4[i],1:length(ID_Parity_all)),get)))
                    phe4_rndnon=cbind(phe4_rndnon,phe4_rndnon1)}
               phe4_rndnon=phe4_rndnon[,-1]
               phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_rndnon,phe4_scc305,phe4_scs305)
           }
        }
   )}
   if(exists("phe11") & exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   phe_final$propercent=phe_final[[ncol(phe11)+2]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","proyield","SCC","SCS","fatpercent","propercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","fatyiled","proyield","SCC","SCS","fatpercent","propercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","proyield","SCC","SCS","fatpercent","propercent")
   }else if(exists("phe11") & exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe33[[ncol(phe33)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   phe_final$propercent=phe_final[[ncol(phe11)+2]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","proyield","fatpercent","propercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","fatyiled","proyield","fatpercent","propercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","proyield","fatpercent","propercent")
   }else if(exists("phe11") & exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","fatpercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","fatyiled","fatpercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","fatpercent")
   }else if(exists("phe11") & !exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe33[[ncol(phe33)]])
   phe_final$propercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","proyield","propercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","proyield","propercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","proyiled","propercent")
   }else if(exists("phe11") & exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","SCC","SCS","fatpercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","fatyield","SCC","SCS","fatpercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","SCC","SCS","fatpercent")
   }else if(exists("phe11") & !exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$propercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","proyiled","SCC","SCS","propercent")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","proyield","SCC","SCS","propercent")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","proyiled","SCC","SCS","propercent")
   }else if(exists("phe11") & !exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","SCC","SCS")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled","SCC","SCS")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","SCC","SCS")
   }else if(exists("phe11") & !exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=phe11
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"milkyiled")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits="milkyiled"
   }else if(!exists("phe11") & exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe22,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","SCC","SCS")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"fatyiled","SCC","SCS")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","SCC","SCS")
   }else if(!exists("phe11") & exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=phe22
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"fatyiled")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled")
   }else if(!exists("phe11") & !exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe33,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"proyiled","SCC","SCS")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"proyiled","SCC","SCS")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("proyiled","SCC","SCS")
   }else if(!exists("phe11") & !exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=phe33
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"proyiled")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"proyiled")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("proyiled")
   }else if(!exists("phe11") & exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe22,phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","proyield","SCC","SCS")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"fatyiled","proyield","SCC","SCS")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","proyield","SCC","SCS")
   }else if(!exists("phe11") & exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe22,phe33[[ncol(phe33)]])
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","proyield")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"fatyiled","proyield")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","proyield")
   }else if(!exists("phe11") & !exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=phe44
   if(length(model$rnd_non)==1 & ("ID"%in%model$rnd_non | "id"%in%model$rnd_non)){
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"SCC","SCS")
   }else{
       names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,model$rnd_non[-which(model$rnd_non=="ID" | model$rnd_non=="id")],"SCC","SCS")
   }
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("SCC","SCS")
   }
}else if(is.null(model$rnd_non)){
   walk(model$traits,function(x){
        if(str_detect(x,fixed("milkyield",ignore_case = T))){
           phe1=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,x))}
           phe1=phe1[order(phe1[[1]],phe1$Parity,phe1$DIM),]
           phe1$ID_Parity=paste(setDF(phe1)[[1]],setDF(phe1)$Parity,sep="")
           ID_Parity_all=unique(phe1$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("mydf",i),phe1[phe1$ID_Parity==ID_Parity_all[i],])
               a2=sum(get(paste0("mydf",i))$DIM<305)
               if(get(paste0("mydf",i))$Parity[1]==1){
                   milk1=get(paste0("mydf",i))$DIM[1]*get(paste0("mydf",i))[[x]][1]*(0.605+0.0435*sqrt(get(paste0("mydf",i))$DIM[1]))         ##第一个矫正时间段
                 if(get(paste0("mydf",i))$DIM[1]>40){
                             for(j in 2:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                 }else if(get(paste0("mydf",i))$DIM[1]<40){
                       if(get(paste0("mydf",i))$DIM[2]>40){
                             milk2=0.5*(0.998+1.23/((get(paste0("mydf",i))$DIM[1])^2)+0.0113*((get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])+1)/(get(paste0("mydf",i))$DIM[1]))*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                       }else if(get(paste0("mydf",i))$DIM[2]<40){
                             milk2=0.5*(0.998+1.23/((get(paste0("mydf",i))$DIM[2])^2)+0.0113*((get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])+1)/(get(paste0("mydf",i))$DIM[2]))*(get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])*(get(paste0("mydf",i))[[x]][3]+get(paste0("mydf",i))[[x]][2])+0.5*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:a2-1){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j+1]-get(paste0("mydf",i))$DIM[j])*(get(paste0("mydf",i))[[x]][j+1]+get(paste0("mydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.071*(305-get(paste0("mydf",i))$DIM[a2])/(48.3-0.071*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2),get))))
                       }
                   }
                 }
               }else if(get(paste0("mydf",i))$Parity[1]>1){
                   milk1=get(paste0("mydf",i))$DIM[1]*get(paste0("mydf",i))[[x]][1]*(0.635+0.0435*sqrt(get(paste0("mydf",i))$DIM[1]))         ##第一个矫正时间段
                 if(get(paste0("mydf",i))$DIM[1]>40){
                             for(j in 2:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                  assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                 }else if(get(paste0("mydf",i))$DIM[1]<40){
                       if(get(paste0("mydf",i))$DIM[2]>40){
                             milk2=0.5*(1.001-0.00042*(get(paste0("mydf",i))$DIM[1])+0.0109*((get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])+1)/(get(paste0("mydf",i))$DIM[1]))*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:a2){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j]-get(paste0("mydf",i))$DIM[j-1])*(get(paste0("mydf",i))[[x]][j]+get(paste0("mydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                 assign(paste0("milk",a2+1),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2+1), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2+1),get))))
                       }else if(get(paste0("mydf",i))$DIM[2]<40){
                             milk2=0.5*(1.001-0.00042*(get(paste0("mydf",i))$DIM[2])+0.0109*((get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])+1)/(get(paste0("mydf",i))$DIM[2]))*(get(paste0("mydf",i))$DIM[3]-get(paste0("mydf",i))$DIM[2])*(get(paste0("mydf",i))[[x]][3]+get(paste0("mydf",i))[[x]][2])+0.5*(get(paste0("mydf",i))$DIM[2]-get(paste0("mydf",i))$DIM[1])*(get(paste0("mydf",i))[[x]][2]+get(paste0("mydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:a2-1){
                                  assign(paste0("milk",j),0.5*(get(paste0("mydf",i))$DIM[j+1]-get(paste0("mydf",i))$DIM[j])*(get(paste0("mydf",i))[[x]][j+1]+get(paste0("mydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("mydf",i))$DIM[a2+1])){
                                 assign(paste0("milk",a2),(get(paste0("mydf",i))[[x]][a2]+0.5*(get(paste0("mydf",i))[[x]][a2+1]-get(paste0("mydf",i))[[x]][a2])/(get(paste0("mydf",i))$DIM[a2+1]-get(paste0("mydf",i))$DIM[a2])*(306-get(paste0("mydf",i))$DIM[a2]))*(305-get(paste0("mydf",i))$DIM[a2]))       ##第三个矫正时间段
                             }else{assign(paste0("milk",a2), get(paste0("mydf",i))[[x]][a2]*(305-get(paste0("mydf",i))$DIM[a2])*(1-0.5*0.144*(305-get(paste0("mydf",i))$DIM[a2])/(71-0.144*get(paste0("mydf",i))$DIM[a2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("my_305all",i),sum(unlist(lapply(paste0("milk",1:a2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("mydf",i))[[model$rnd_gen]][a2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("mydf",i))[[model$fix_fct[k]]][a2])}
               assign(paste0("Parity",i),get(paste0("mydf",i))[["Parity"]][a2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("mydf",i))[[model$fix_cov[k]]][a2])}
           }
           phe1_my305=data.table(unlist(lapply(paste0("my_305all",1:length(ID_Parity_all)),get)))
           phe1_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe1_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe1_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe1_fixfct=cbind(phe1_fixfct,phe1_fixfct1)}
           phe1_fixfct=phe1_fixfct[,-1]
           phe1_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe1_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe1_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe1_fixcov=cbind(phe1_fixcov,phe1_fixcov1)}
           phe1_fixcov=phe1_fixcov[,-1]
           phe11=cbind(phe1_rndgen,phe1_fixfct,phe1_Parity,phe1_fixcov,phe1_my305)
        }
        if(str_detect(x,fixed("fatyield",ignore_case = T))){
           phe2=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,x))}
           phe2=phe2[order(phe2[[1]],phe2$Parity,phe2$DIM),]
           phe2$ID_Parity=paste(setDF(phe2)[[1]],setDF(phe2)$Parity,sep="")
           ID_Parity_all=unique(phe2$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("fydf",i),phe2[phe2$ID_Parity==ID_Parity_all[i],])
               b2=sum(get(paste0("fydf",i))$DIM<305)
               if(get(paste0("fydf",i))$Parity[1]==1){
                   fat1=get(paste0("fydf",i))$DIM[1]*get(paste0("fydf",i))[[x]][1]*(0.235+0.239*sqrt(get(paste0("fydf",i))$DIM[1])-0.0225*get(paste0("fydf",i))$DIM[1]+0.000069*(get(paste0("fydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("fydf",i))$DIM[1]>40){
                             for(j in 2:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                 }else if(get(paste0("fydf",i))$DIM[1]<40){
                       if(get(paste0("fydf",i))$DIM[2]>40){
                             fat2=0.5*(0.998+1.23/((get(paste0("fydf",i))$DIM[1])^2)+0.0113*((get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])+1)/(get(paste0("fydf",i))$DIM[1]))*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                       }else if(get(paste0("fydf",i))$DIM[2]<40){
                             fat2=0.5*(0.998+1.23/((get(paste0("fydf",i))$DIM[2])^2)+0.0113*((get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])+1)/(get(paste0("fydf",i))$DIM[2]))*(get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])*(get(paste0("fydf",i))[[x]][3]+get(paste0("fydf",i))[[x]][2])+0.5*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:b2-1){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j+1]-get(paste0("fydf",i))$DIM[j])*(get(paste0("fydf",i))[[x]][j+1]+get(paste0("fydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0025*(305-get(paste0("fydf",i))$DIM[b2])/(2.03-0.0025*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2),get))))
                       }
                   }
                 }
               }else if(get(paste0("fydf",i))$Parity[1]>1){
                   fat1=get(paste0("fydf",i))$DIM[1]*get(paste0("fydf",i))[[x]][1]*(0.476+0.146*sqrt(get(paste0("fydf",i))$DIM[1])-0.0115*get(paste0("fydf",i))$DIM[1]+0.000038*(get(paste0("fydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("fydf",i))$DIM[1]>40){
                             for(j in 2:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                  assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                 }else if(get(paste0("fydf",i))$DIM[1]<40){
                       if(get(paste0("fydf",i))$DIM[2]>40){
                             fat2=0.5*(1.001-0.00042*(get(paste0("fydf",i))$DIM[1])+0.0109*((get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])+1)/(get(paste0("fydf",i))$DIM[1]))*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:b2){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j]-get(paste0("fydf",i))$DIM[j-1])*(get(paste0("fydf",i))[[x]][j]+get(paste0("fydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                 assign(paste0("fat",b2+1),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2+1), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2+1),get))))
                       }else if(get(paste0("fydf",i))$DIM[2]<40){
                             fat2=0.5*(1.001-0.00042*(get(paste0("fydf",i))$DIM[2])+0.0109*((get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])+1)/(get(paste0("fydf",i))$DIM[2]))*(get(paste0("fydf",i))$DIM[3]-get(paste0("fydf",i))$DIM[2])*(get(paste0("fydf",i))[[x]][3]+get(paste0("fydf",i))[[x]][2])+0.5*(get(paste0("fydf",i))$DIM[2]-get(paste0("fydf",i))$DIM[1])*(get(paste0("fydf",i))[[x]][2]+get(paste0("fydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:b2-1){
                                  assign(paste0("fat",j),0.5*(get(paste0("fydf",i))$DIM[j+1]-get(paste0("fydf",i))$DIM[j])*(get(paste0("fydf",i))[[x]][j+1]+get(paste0("fydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("fydf",i))$DIM[b2+1])){
                                 assign(paste0("fat",b2),(get(paste0("fydf",i))[[x]][b2]+0.5*(get(paste0("fydf",i))[[x]][b2+1]-get(paste0("fydf",i))[[x]][b2])/(get(paste0("fydf",i))$DIM[b2+1]-get(paste0("fydf",i))$DIM[b2])*(306-get(paste0("fydf",i))$DIM[b2]))*(305-get(paste0("fydf",i))$DIM[b2]))       ##第三个矫正时间段
                             }else{assign(paste0("fat",b2), get(paste0("fydf",i))[[x]][b2]*(305-get(paste0("fydf",i))$DIM[b2])*(1-0.5*0.0052*(305-get(paste0("fydf",i))$DIM[b2])/(2.78-0.0052*get(paste0("fydf",i))$DIM[b2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("fy_305all",i),sum(unlist(lapply(paste0("fat",1:b2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("fydf",i))[[model$rnd_gen]][b2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("fydf",i))[[model$fix_fct[k]]][b2])}
               assign(paste0("Parity",i),get(paste0("fydf",i))[["Parity"]][b2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("fydf",i))[[model$fix_cov[k]]][b2])}
           }
           phe2_fy305=data.table(unlist(lapply(paste0("fy_305all",1:length(ID_Parity_all)),get)))
           phe2_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe2_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe2_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe2_fixfct=cbind(phe2_fixfct,phe2_fixfct1)}
           phe2_fixfct=phe2_fixfct[,-1]
           phe2_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe2_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe2_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe2_fixcov=cbind(phe2_fixcov,phe2_fixcov1)}
           phe2_fixcov=phe2_fixcov[,-1]
           phe22=cbind(phe2_rndgen,phe2_fixfct,phe2_Parity,phe2_fixcov,phe2_fy305)
        }
        if(str_detect(x,fixed("proyield",ignore_case = T))){
           phe3=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,x))}
           phe3=phe3[order(phe3[[1]],phe3$Parity,phe3$DIM),]
           phe3$ID_Parity=paste(setDF(phe3)[[1]],setDF(phe3)$Parity,sep="")
           ID_Parity_all=unique(phe3$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("pydf",i),phe3[phe3$ID_Parity==ID_Parity_all[i],])
               c2=sum(get(paste0("pydf",i))$DIM<305)
               if(get(paste0("pydf",i))$Parity[1]==1){
                   pro1=get(paste0("pydf",i))$DIM[1]*get(paste0("pydf",i))[[x]][1]*(0.136+0.316*sqrt(get(paste0("pydf",i))$DIM[1])-0.0351*get(paste0("pydf",i))$DIM[1]+0.00013*(get(paste0("pydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("pydf",i))$DIM[1]>40){
                             for(j in 2:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                 }else if(get(paste0("pydf",i))$DIM[1]<40){
                       if(get(paste0("pydf",i))$DIM[2]>40){
                             pro2=0.5*(0.998+1.23/((get(paste0("pydf",i))$DIM[1])^2)+0.0113*((get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])+1)/(get(paste0("pydf",i))$DIM[1]))*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                       }else if(get(paste0("pydf",i))$DIM[2]<40){
                             pro2=0.5*(0.998+1.23/((get(paste0("pydf",i))$DIM[2])^2)+0.0113*((get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])+1)/(get(paste0("pydf",i))$DIM[2]))*(get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])*(get(paste0("pydf",i))[[x]][3]+get(paste0("pydf",i))[[x]][2])+0.5*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:c2-1){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j+1]-get(paste0("pydf",i))$DIM[j])*(get(paste0("pydf",i))[[x]][j+1]+get(paste0("pydf",i))[[x]][j]))}  
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0019*(305-get(paste0("pydf",i))$DIM[c2])/(1.67-0.0019*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2),get))))
                       }
                   }
                 }
               }else if(get(paste0("pydf",i))$Parity[1]>1){
                   pro1=get(paste0("pydf",i))$DIM[1]*get(paste0("pydf",i))[[x]][1]*(0.177+0.324*sqrt(get(paste0("pydf",i))$DIM[1])-0.0366*get(paste0("pydf",i))$DIM[1]+0.000141*(get(paste0("pydf",i))$DIM[1])^2)         ##第一个矫正时间段
                 if(get(paste0("pydf",i))$DIM[1]>40){
                             for(j in 2:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                  assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                  assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                 }else if(get(paste0("pydf",i))$DIM[1]<40){
                       if(get(paste0("pydf",i))$DIM[2]>40){
                             pro2=0.5*(1.001-0.00042*(get(paste0("pydf",i))$DIM[1])+0.0109*((get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])+1)/(get(paste0("pydf",i))$DIM[1]))*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])       ##第二个矫正时间段
                             for(j in 3:c2){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j]-get(paste0("pydf",i))$DIM[j-1])*(get(paste0("pydf",i))[[x]][j]+get(paste0("pydf",i))[[x]][j-1]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                 assign(paste0("pro",c2+1),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2+1), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2+1),get))))
                       }else if(get(paste0("pydf",i))$DIM[2]<40){
                             pro2=0.5*(1.001-0.00042*(get(paste0("pydf",i))$DIM[2])+0.0109*((get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])+1)/(get(paste0("pydf",i))$DIM[2]))*(get(paste0("pydf",i))$DIM[3]-get(paste0("pydf",i))$DIM[2])*(get(paste0("pydf",i))[[x]][3]+get(paste0("pydf",i))[[x]][2])+0.5*(get(paste0("pydf",i))$DIM[2]-get(paste0("pydf",i))$DIM[1])*(get(paste0("pydf",i))[[x]][2]+get(paste0("pydf",i))[[x]][1])    ##第二个矫正时间段
                             for(j in 3:c2-1){
                                  assign(paste0("pro",j),0.5*(get(paste0("pydf",i))$DIM[j+1]-get(paste0("pydf",i))$DIM[j])*(get(paste0("pydf",i))[[x]][j+1]+get(paste0("pydf",i))[[x]][j]))}     ##中间无需矫正的时间段
                             if(!is.na(get(paste0("pydf",i))$DIM[c2+1])){
                                 assign(paste0("pro",c2),(get(paste0("pydf",i))[[x]][c2]+0.5*(get(paste0("pydf",i))[[x]][c2+1]-get(paste0("pydf",i))[[x]][c2])/(get(paste0("pydf",i))$DIM[c2+1]-get(paste0("pydf",i))$DIM[c2])*(306-get(paste0("pydf",i))$DIM[c2]))*(305-get(paste0("pydf",i))$DIM[c2]))       ##第三个矫正时间段
                             }else{assign(paste0("pro",c2), get(paste0("pydf",i))[[x]][c2]*(305-get(paste0("pydf",i))$DIM[c2])*(1-0.5*0.0045*(305-get(paste0("pydf",i))$DIM[c2])/(2.39-0.0045*get(paste0("pydf",i))$DIM[c2])))}       ##305d后测定日不存在，第四个矫正时间段
                                 assign(paste0("py_305all",i),sum(unlist(lapply(paste0("pro",1:c2),get))))
                       }
                   }
               }
               assign(paste0(model$rnd_gen,i),get(paste0("pydf",i))[[model$rnd_gen]][c2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("pydf",i))[[model$fix_fct[k]]][c2])}
               assign(paste0("Parity",i),get(paste0("pydf",i))[["Parity"]][c2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("pydf",i))[[model$fix_cov[k]]][c2])}
           }
           phe3_py305=data.table(unlist(lapply(paste0("py_305all",1:length(ID_Parity_all)),get)))
           phe3_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe3_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe3_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe3_fixfct=cbind(phe3_fixfct,phe3_fixfct1)}
           phe3_fixfct=phe3_fixfct[,-1]
           phe3_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe3_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe3_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe3_fixcov=cbind(phe3_fixcov,phe3_fixcov1)}
           phe3_fixcov=phe3_fixcov[,-1]
           phe33=cbind(phe3_rndgen,phe3_fixfct,phe3_Parity,phe3_fixcov,phe3_py305)
        }
        if(str_detect(x,fixed("SCS",ignore_case = T))){
           phe4=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,x))}
           phe4=phe4[order(phe4[[1]],phe4$Parity,phe4$DIM),]
           phe4$ID_Parity=paste(setDF(phe4)[[1]],setDF(phe4)$Parity,sep="")
           ID_Parity_all=unique(phe4$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("scdf",i),phe4[phe4$ID_Parity==ID_Parity_all[i],])
               d2=sum(get(paste0("scdf",i))$DIM<305)
               assign(paste0("scc",i),2^(get(paste0("scdf",i))[[x]]-3)*100)
               assign(paste0("scc_305",i),mean(get(paste0("scc",i))))              ##SCS与SCC直接取平均值
               assign(paste0("scs_305",i),log2(get(paste0("scc_305",i))/100)+3)
               assign(paste0(model$rnd_gen,i),get(paste0("scdf",i))[[model$rnd_gen]][d2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("scdf",i))[[model$fix_fct[k]]][d2])}
               assign(paste0("Parity",i),get(paste0("scdf",i))[["Parity"]][d2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("scdf",i))[[model$fix_cov[k]]][d2])}
           }
           phe4_scc305=data.table(unlist(lapply(paste0("scc_305",1:length(ID_Parity_all)),get)))
           phe4_scs305=data.table(unlist(lapply(paste0("scs_305",1:length(ID_Parity_all)),get)))
           phe4_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe4_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe4_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe4_fixfct=cbind(phe4_fixfct,phe4_fixfct1)}
           phe4_fixfct=phe4_fixfct[,-1]
           phe4_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe4_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe4_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe4_fixcov=cbind(phe4_fixcov,phe4_fixcov1)}
           phe4_fixcov=phe4_fixcov[,-1]
           phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_scc305,phe4_scs305)
        }
        if(str_detect(x,fixed("SCC",ignore_case = T))){
           phe4=select(phe,c(model$rnd_gen,model$fix_fct,"DIM","Parity",model$fix_cov,x))}
           phe4=phe4[order(phe4[[1]],phe4$Parity,phe4$DIM),]
           phe4$ID_Parity=paste(setDF(phe4)[[1]],setDF(phe4)$Parity,sep="")
           ID_Parity_all=unique(phe4$ID_Parity)
           for(i in 1:length(ID_Parity_all)){
               assign(paste0("scdf",i),phe4[phe4$ID_Parity==ID_Parity_all[i],])
               d2=sum(get(paste0("scdf",i))$DIM<305)
               assign(paste0("scc_305",i),mean(get(paste0("scdf",i))[[x]]))              ##SCS与SCC直接取平均值
               assign(paste0("scs_305",i),log2(get(paste0("scc_305",i))/100)+3)
               assign(paste0(model$rnd_gen,i),get(paste0("scdf",i))[[model$rnd_gen]][d2])
               for(k in 1:length(model$fix_fct)){
                    assign(paste0(model$fix_fct[k],i),get(paste0("scdf",i))[[model$fix_fct[k]]][d2])}
               assign(paste0("Parity",i),get(paste0("scdf",i))[["Parity"]][d2])
               for(k in 1:length(model$fix_cov)){
                     assign(paste0(model$fix_cov[k],i),get(paste0("scdf",i))[[model$fix_cov[k]]][d2])}
           }
           phe4_scc305=data.table(unlist(lapply(paste0("scc_305",1:length(ID_Parity_all)),get)))
           phe4_scs305=data.table(unlist(lapply(paste0("scs_305",1:length(ID_Parity_all)),get)))
           phe4_rndgen=data.table(unlist(lapply(paste0(model$rnd_gen,1:length(ID_Parity_all)),get)))
           phe4_fixfct=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_fct)){
               phe4_fixfct1=data.table(unlist(lapply(paste0(model$fix_fct[i],1:length(ID_Parity_all)),get)))
               phe4_fixfct=cbind(phe4_fixfct,phe4_fixfct1)}
           phe4_fixfct=phe4_fixfct[,-1]
           phe4_Parity=data.table(unlist(lapply(paste0("Parity",1:length(ID_Parity_all)),get)))
           phe4_fixcov=data.table(rep(NA,length(ID_Parity_all)))
           for(i in 1:length(model$fix_cov)){
               phe4_fixcov1=data.table(unlist(lapply(paste0(model$fix_cov[i],1:length(ID_Parity_all)),get)))
               phe4_fixcov=cbind(phe4_fixcov,phe4_fixcov1)}
           phe4_fixcov=phe4_fixcov[,-1]
           phe44=cbind(phe4_rndgen,phe4_fixfct,phe4_Parity,phe4_fixcov,phe4_scc305,phe4_scs305)
        }
   )}
   if(exists("phe11") & exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   phe_final$propercent=phe_final[[ncol(phe11)+2]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","proyield","SCC","SCS","fatpercent","propercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","proyield","SCC","SCS","fatpercent","propercent")
   }else if(exists("phe11") & exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe33[[ncol(phe33)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   phe_final$propercent=phe_final[[ncol(phe11)+2]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","proyield","fatpercent","propercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","proyield","fatpercent","propercent")
   }else if(exists("phe11") & exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","fatpercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","fatpercent")
   }else if(exists("phe11") & !exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe11,phe33[[ncol(phe33)]])
   phe_final$propercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","proyield","propercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","proyiled","propercent")
   }else if(exists("phe11") & exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe22[[ncol(phe22)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$fatpercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","fatyiled","SCC","SCS","fatpercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","fatyiled","SCC","SCS","fatpercent")
   }else if(exists("phe11") & !exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   phe_final$propercent=phe_final[[ncol(phe11)+1]]/phe_final[[ncol(phe11)]]*100
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","proyiled","SCC","SCS","propercent")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","proyiled","SCC","SCS","propercent")
   }else if(exists("phe11") & !exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe11,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled","SCC","SCS")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("milkyiled","SCC","SCS")
   }else if(exists("phe11") & !exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=phe11
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"milkyiled")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits="milkyiled"
   }else if(!exists("phe11") & exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=cbind(phe22,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","SCC","SCS")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","SCC","SCS")
   }else if(!exists("phe11") & exist("phe22") & !exist("phe33") & !exists("phe44")){
   phe_final=phe22
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled")
   }else if(!exists("phe11") & !exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe33,phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"proyiled","SCC","SCS")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("proyiled","SCC","SCS")
   }else if(!exists("phe11") & !exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=phe33
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"proyiled")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("proyiled")
   }else if(!exists("phe11") & exist("phe22") & exist("phe33") & exists("phe44")){
   phe_final=cbind(phe22,phe33[[ncol(phe33)]],phe44[[ncol(phe44)-1]],phe44[[ncol(phe44)]])
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","proyield","SCC","SCS")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","proyield","SCC","SCS")
   }else if(!exists("phe11") & exist("phe22") & exist("phe33") & !exists("phe44")){
   phe_final=cbind(phe22,phe33[[ncol(phe33)]])
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"fatyiled","proyield")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("fatyiled","proyield")
   }else if(!exists("phe11") & !exist("phe22") & !exist("phe33") & exists("phe44")){
   phe_final=phe44
   names(phe_final)=c(model$rnd_gen,model$fix_fct,"Parity",model$fix_cov,"SCC","SCS")
   write.table(phe_final,file="phe_after_adjust.txt",row.names = FALSE,sep="\t",quote = FALSE)
   model$traits=c("SCC","SCS")
   }
}


  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec<-c(model$fix_fct,model$rnd_non)
  fct<-map_dfc(fct_vec,function(x)as.numeric(as.factor(phe_final[[x]])))
  fct_name<-map(fct_vec,function(x){
    name<-levels(as.factor(phe_final[[x]]))
    data.frame(Name=name,Code=1:length(name))
  })
  names(fct_name)<-fct_vec
  names(fct)<-fct_vec
  
  ## 提取实数变量——协变量+性状表型值
  rl<-select(phe_final,c(model$fix_cov,model$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe_final[[x]],ped[[1]])]
    })
    rnd_gen_name<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }else{
    rnd_gen<-map_dfc(model$rnd_gen,function(x){
      ord[match(phe_final[[x]],fam_new[[1]])]
    })
    rnd_gen_name<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen)<-model$rnd_gen
    
  }
  
  rnd_gen[is.na(rnd_gen)]<--9999
  fct[is.na(fct)]<--9999
  rl[is.na(rl)]<--9999
  
  # 开始进行计算
  int_name<-c(model$rnd_gen,model$fix_fct,model$rnd_non)
  int_code<-c(match(model$fix_fct,int_name),
              match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  rnd_code<-c(match(model$rnd_gen,int_name),
              match(model$rnd_non,int_name))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name[match(fam$V1,rnd_gen_name$Name),'Code']
    grm2<-data.frame(fam_id[grm$V1],fam_id[grm$V2],grm$V4)
    return(grm2)
  }
  
  create_ginv<-function(infile="final.grm.gz",fam=fam_new){
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
  
  ### BLUP
  if(method=="BLUP"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

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
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
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
  # BLUP_VAR
  }else if(method=="BLUP_VAR"){
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$VAR_STR 1 PED 2 ASCII ../DMU_PED

$PRIOR VAR_INPUT

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
    
    ### ssGBLUP_VAR
  }else if(method=="ssGBLUP_VAR"){
    grm<-get_grm()
    write.table(grm,file="GMAT",row.names = F,col.names = F,quote=F)
    write.table(unique(grm[[1]]),file="GENO_ID",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
    
    ### GBLUP_VAR
  }else if(method=="GBLUP_VAR"){
    ginv<-create_ginv()
    write.table(ginv,file="GINV",row.names = F,col.names = F,quote=F)
    
    walk(model$traits,function(x){
      rl_tmp<-c(model$fix_cov,x)
      rl_dat<-select(rl,all_of(rl_tmp))
      dmu_dat<-cbind(rnd_gen,fct,rl_dat)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat<-dmu_dat[order(dmu_dat[[1]]),]
      write.table(dmu_dat,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      var_comps<-var_comp_input[which(model$traits==x),]
      var_output<-data.frame(V1=1:ncol(var_comp_input),V2=1,V3=1,V4=map_dbl(var_comps,~.x))
      write.table(var_output,file=paste0(outdir,"/VAR_INPUT"),row.names = F,
                  col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 11 9 0 0

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

$PRIOR VAR_INPUT

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
  }else{
    
    converge_msg<-map_dfr(model$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
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
        if(need_pedigree){
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
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model$rnd_non)>0){
          if(model$rnd_non=="ID"){
           non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          }else{
            non_genetic_sol<-map_dfr(model$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model$rnd_gen)+which(model$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name[[x]]$Name,Code=fct_name[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          }
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
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass){
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         var_summary_msg)
    
  }else{
    oneLine<-str_replace(oneLine, "var_data_summary_is_here", 
                         "**[ERROR: ] Something is wrong, no variance component data summary message is shown. Please check your data first!**")
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no converge check message is shown. Please check your data first!\\*\\*", 
                         paste0("* ",paste(converge_msg$MSG,collapse = "#' * ")))
    
    
  }
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass & Success!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }
  writeLines(oneLine,con_out)
}

close(con)
close(con_out)

system(paste0("cp ",report_path,"/Rmd.css ."))
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


smtp <- server(host = "smtp.qq.com",
               port = 465,
               username = from_email,
               password = "cwdqmmxjoimhbbac")

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
    message("[WARNING:] Emailing the results caused a wanning!\n")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  }
)
}









