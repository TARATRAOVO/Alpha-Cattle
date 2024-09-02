library(data.table)
library(pedigree)
library(tidyverse)
library(emayili)

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
model1<-res$model1
model2<-res$model2
dir<-res$dir
phe_file<-res$phe_file
ped_file<-res$ped_file
gen_file<-res$gen_file
DMU1<-"/srv/shiny-server/dmu/dmu1"
DMUAI<-"/srv/shiny-server/dmu/dmuai"
DMU4<-"/srv/shiny-server/dmu/dmu4"
plink<-"/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink"
gcta<-"/public/home/chenjc/miniconda3/pkgs/gcta-1.93.2beta-h9ee0642_1/bin/gcta64"
from_email<-"852322127@qq.com"
to_email<-res$email
task<-res$task

# 定义一个确定是否需要基因型或者系谱数据的变量
need_genotype<-TRUE
if(method=="BLUP"){
  need_genotype<-FALSE
}

need_pedigree<-TRUE
if(method=="GBLUP"){
  need_pedigree<-FALSE
}

# 定义一个是否需要方差组分文件的变量
need_var_comp<-TRUE
if(!str_detect(method,"VAR$")){
  need_var_comp<-FALSE
}

# 检查model是否正常，我们的model必须包含固定效应（至少有一个mean）、随机效应
model_check_pass<-TRUE

if(length(model1$fix_fct)==0 | length(model1$rnd_gen)==0){
  model_check_pass<-FALSE
}else{
  model1$rnd_gen<-model1$rnd_gen[1]
} ###这里如果model check不通过，是不是会没有报错信息。

if(length(model2$fix_fct)!=0 & length(model2$rnd_gen)!=0){
  model2$rnd_gen<-model2$rnd_gen[1]
}


# 数据读入、质控与编码
phe<-fread(phe_file,fill=T)
phe<-unique(phe)

if(need_pedigree){
  ped<-fread(ped_file,fill=T)
  ped<-unique(ped)
}

if(need_genotype){
  fam<-fread(res$gen_file[str_detect(res$gen_file,"fam$")][1],h=F)
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
  if(method=="BLUP"){
    phe_not_in_ped<-phe[[1]][!phe[[1]]%in%ped[[1]]]
    res$phe_not_in_ped<-phe_not_in_ped
    if(length(phe_not_in_ped)>0){
      msg<-sprintf("%d individual(s) with phenotype not in pedigree were removed. You can find removed ID in `phe_id_not_in_ped.txt` file.\n", length(phe_not_in_ped))
    }else{
      msg<-"all individuals with phenotype are in pedigree.\n"
    }
  }
  
  if(method=="GBLUP"){
    phe_not_in_fam<-phe[[1]][!phe[[1]]%in%fam$V1]
    res$phe_not_in_fam<-phe_not_in_fam
    if(length(phe_not_in_fam)>0){
      msg<-sprintf("%d individual(s) with phenotype not in genotype data were removed. You can find removed ID in `phe_id_not_in_fam.txt` file.\n", length(phe_not_in_fam))
    }else{
      msg<-"all individuals with phenotype are in genotype data.\n"
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
  
  
  phe_summary_num1<-map_dbl(model1$traits,function(x)sum(!is.na(phe[[x]])))
  phe_summary_mean1<-map_dbl(model1$traits,function(x)mean(phe[[x]],na.rm=T))
  phe_summary_sd1<-map_dbl(model1$traits,function(x)sd(phe[[x]],na.rm=T))

  phe_summary1<-data.frame(
    Trait=model1$traits,
    Size=phe_summary_num1,
    Mean=phe_summary_mean1,
    SD=phe_summary_sd1
  )
  
   phe_summary<-phe_summary1

if(!is.null(model2$traits)){
  phe_summary_num2<-map_dbl(model2$traits,function(x)sum(!is.na(phe[[x]])))
  phe_summary_mean2<-map_dbl(model2$traits,function(x)mean(phe[[x]],na.rm=T))
  phe_summary_sd2<-map_dbl(model2$traits,function(x)sd(phe[[x]],na.rm=T))
  
    phe_summary2<-data.frame(
      Trait=model2$traits,
      Size=phe_summary_num2,
      Mean=phe_summary_mean2,
      SD=phe_summary_sd2
  )
    phe_summary<-rbind(phe_summary1,phe_summary2)
}


  if(nrow(phe)<=0 | sum(phe_summary_num1)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }

if(!is.null(model2$traits)){
  if(nrow(phe)<=0 | sum(phe_summary_num1)==0 | sum(phe_summary_num2)==0){
    phenotype_check_pass<-FALSE
    phenotype_summary_msg<-"**[ERROR:] After data filter, we found no phenotype record is left in the data, so the whole analysis is stopped.**"
  }
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
          system(paste0(gcta, " --bfile ", geno_names[i], ".cln --maf 0 --make-grm --out ", geno_names[i], " --make-grm-alg 1"))
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
    dmu_ped[is.na(dmu_ped)]<-0
    
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
Success1<-0
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass){
  
  ## 编码固定效应中的因子与非遗传随机效应
  fct_vec1<-c(model1$fix_fct,model1$rnd_non)
  fct1<-map_dfc(fct_vec1,function(x)as.numeric(as.factor(phe[[x]])))
  fct_name1<-map(fct_vec1,function(x){
    name1<-levels(as.factor(phe[[x]]))
    data.frame(Name=name1,Code=1:length(name1))
  })
  names(fct_name1)<-fct_vec1
  names(fct1)<-fct_vec1
  
  ## 提取实数变量——协变量+性状表型值
  rl1<-select(phe,c(model1$fix_cov,model1$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen1<-map_dfc(model1$rnd_gen,function(x){
      ord[match(phe[[x]],ped[[1]])]
    })
    rnd_gen_name1<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen1)<-model1$rnd_gen
    
  }else{
    rnd_gen1<-map_dfc(model1$rnd_gen,function(x){
      ord[match(phe[[x]],fam_new[[1]])]
    })
    rnd_gen_name1<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen1)<-model1$rnd_gen
    
  }
  
  rnd_gen1[is.na(rnd_gen1)]<--999999
  fct1[is.na(fct1)]<--999999
  rl1[is.na(rl1)]<--999999
  
  # 开始进行计算
  int_name1<-c(model1$rnd_gen,model1$fix_fct,model1$rnd_non)
  int_code1<-c(match(model1$fix_fct,int_name1),
              match(model1$rnd_gen,int_name1),
              match(model1$rnd_non,int_name1))
  rnd_code1<-c(match(model1$rnd_gen,int_name1),
              match(model1$rnd_non,int_name1))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm1<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name1[match(fam$V1,rnd_gen_name1$Name),'Code']
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
    fam_id<-rnd_gen_name1[match(fam$V1,rnd_gen_name1$Name),'Code']
    #fam_id<-1:nrow(fam)
    ginv_dat<-data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
    return(ginv_dat)
  }
  
  ### BLUP
  if(method=="BLUP"){
    walk(model1$traits,function(x){
      rl_tmp1<-c(model1$fix_cov,x)
      rl_dat1<-select(rl1,all_of(rl_tmp1))
      dmu_dat1<-cbind(rnd_gen1,fct1,rl_dat1)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat1,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 1 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1
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
",method,x,length(int_name1),length(rl_tmp1),
                  paste(int_name1,collapse = " "),
                  paste(rl_tmp1,collapse = " "),
                  length(rl_tmp1),length(int_code1),
                  paste(int_code1,collapse = " "),
                  length(rnd_code1),paste(1:length(rnd_code1),collapse = " "),
                  length(model1$fix_cov), paste(1:length(model1$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm1()
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
1
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
    
    walk(model1$traits,function(x){
      rl_tmp1<-c(model1$fix_cov,x)
      rl_dat1<-select(rl1,all_of(rl_tmp1))
      dmu_dat1<-cbind(rnd_gen1,fct1,rl_dat1)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat1<-dmu_dat1[order(dmu_dat1[[1]]),]
      write.table(dmu_dat1,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 31 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1
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
",method,x,length(int_name1),length(rl_tmp1),
                  paste(int_name1,collapse = " "),
                  paste(rl_tmp1,collapse = " "),
                  length(rl_tmp1),length(int_code1),
                  paste(int_code1,collapse = " "),
                  length(rnd_code1),paste(1:length(rnd_code1),collapse = " "),
                  length(model1$fix_cov), paste(1:length(model1$fix_cov),collapse = " ")),
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
1
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
    grm<-get_grm1()
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
1
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
1
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
  
  dmu_log1<-map(model1$traits,function(x){
    outdir<-toupper(x)
    setwd(outdir)
    dmu1_log1<-system(paste0(DMU1," <DMU.DIR>DMU.lst"), intern=T)
    if(!str_detect(method,"VAR$")){
      dmuai_log1<-system(paste0(DMUAI," >>DMU.lst"), intern=T)
    }else{
      dmuai_log1<-system(paste0(DMU4," >>DMU.lst"), intern=T)
    }
    
    setwd("../")
    return(sum(attr(dmu1_log1,'status'),attr(dmuai_log1,'status')))
  })
  names(dmu_log1)<-model1$traits
  
  if(any(dmu_log1>0)){
    Success1=1
    error_trait1<-names(dmu_log1)[dmu_log1>0]
  }
  
  if(all(dmu_log1>0)){
    Success1=2
  }
  
  
  ## 先检验SOL文件是否生成，再检验模型是否收敛
  
  if(!str_detect(method,"VAR$")){
    converge_msg1<-map_dfr(model1$traits, function(x){
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
      converge_msg1<-""
      if(sol_exists & converged){
        converge_msg1<-paste0("Variance components estimation for ",x," has successfully converged.\n")
      }
      if(sol_exists & !converged){
        converge_msg1<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & converged){
        converge_msg1<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & !converged){
        converge_msg1<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg1,SOL=sol_exists, CON=converged)
      return(res)
    })
  }else{
    
    converge_msg1<-map_dfr(model1$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg1<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg1<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg1<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg1,SOL=sol_exists, CON=converged)
      return(res)
    })
  }
  
  if(Success1!=2){
    ## 生成几个结果的表格
    ### Variance component
    
    if(!str_detect(method,"VAR$")){
      var_comp1<-as.data.frame(matrix(NA,length(model1$traits),length(rnd_code1)+2))
      colnames(var_comp1)<-c("Trait",int_name1[rnd_code1],"Residual")
      var_comp1$Trait<-model1$traits
      
      walk(model1$traits,function(x){
        outdir<-toupper(x)
        if(converge_msg1$SOL[converge_msg1$TRT==x] & converge_msg1$CON[converge_msg1$TRT==x]){
          var_com1<-read.table(paste0(outdir,"/PAROUT"),h=F)
          var_comp1[var_comp1$Trait==x,2:ncol(var_comp1)]<<-var_com1$V4
        }
      })
    }else{
      var_comp1<-cbind(model1$traits,var_comp_input)
      colnames(var_comp1)<-c("Trait",int_name1[rnd_code1],"Residual")
    }
    write.table(var_comp1,file="Variance_component1.txt",row.names = F,quote=F,sep="\t")
    
    get_estimate1<-function(trait){
      res<-list()
      res$fix_eff<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      res$genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      if(length(model1$rnd_non)>0){res$non_genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)}
      outdir<-toupper(trait)
      if(converge_msg1$SOL[converge_msg1$TRT==trait] & converge_msg1$CON[converge_msg1$TRT==trait]){
        ### 读取solution
        sol<-fread(paste0(outdir,"/SOL"))
        fix_sol<-sol[sol$V1==2,]
        gen_sol<-sol[sol$V1==4,]
        nongen_sol<-sol[sol$V1==3,]
        if(method=="GBLUP"){
        gen_sol<-sol[sol$V1==3,]
        nongen_sol<-sol[sol$V1==4,]}
        if(length(model1$fix_cov)>0){cov_sol<-sol[sol$V1==1,]}
        
        ### fixed factor effect
        fixed_fct_name<-int_name1[-rnd_code1]
        
        fix_fct_sol<-map_dfr(fixed_fct_name,function(x){
          code<-fct_name1[[x]]
          sol_dat<-fix_sol[fix_sol$V4==which(fixed_fct_name==x),]
          sol_dat2<-sol_dat[match(code$Code,sol_dat$V5),]
          sol_res<-data.frame(Name=x,Value=code$Name,Code=code$Code,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
          return(sol_res)
        })
        
        ### Covariate effect
        if(length(model1$fix_cov)>0){
          cov_name<-model1$fix_cov
          fix_cov_sol<-data.frame(Name=cov_name,Value=1,Code=1:length(cov_name),Estimate=cov_sol$V8,SE=cov_sol$V9)
        }else{
          fix_cov_sol<-NULL
        }
        
        ### Genetic effect
        if(need_pedigree){
          genetic_sol<-map_dfr(model1$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model1$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model1$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model1$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model1$rnd_non)>0){
          non_genetic_sol<-map_dfr(model1$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model1$rnd_gen)+which(model1$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name1[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name1[[x]]$Name,Code=fct_name1[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          non_genetic_sol<-non_genetic_sol[complete.cases(non_genetic_sol),]
          
        }
        
        
        
        res$fix_eff<-rbind(fix_fct_sol,fix_cov_sol)
        res$genetic_sol<-genetic_sol
        if(length(model1$rnd_non)>0){res$non_genetic_sol<-non_genetic_sol}
      }
      return(res)
    }
    
    all_estimate1<-map(model1$traits,get_estimate1)
    names(all_estimate1)<-model1$traits
    
    fix_eff1<-rbindlist(map(model1$traits,function(x)all_estimate1[[x]]$fix_eff),idcol="Trait")
    fix_eff1$Trait<-model1$traits[fix_eff1$Trait]
    
    gen_eff1<-rbindlist(map(model1$traits,function(x)all_estimate1[[x]]$genetic_sol),idcol="Trait")
    gen_eff1$Trait<-model1$traits[gen_eff1$Trait]
    
    if(length(model1$rnd_non)>0){
      nongen_eff1<-rbindlist(map(model1$traits,function(x)all_estimate1[[x]]$non_genetic_sol),idcol="Trait")
      nongen_eff1$Trait<-model1$traits[nongen_eff1$Trait]
    }
    
    fwrite(fix_eff1, file="fixed_effect_estimate1.txt",row.names=F,quote=F,sep="\t")
    fwrite(gen_eff1, file="genetic_effect_estimate1.txt",row.names=F,quote=F,sep="\t")
    if(length(model1$rnd_non)>0){fwrite(nongen_eff1, file="nongenetic_random_effect_estimate1.txt",row.names=F,quote=F,sep="\t")}
  }
}else{
  Success1<-2
}

##以下为体型性状计算流程
Success2<-0
if(!is.null(model2$traits)){
if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass ){

  ## 编码+固定效应中的因子与非遗传随机效应
  fct_vec2<-c(model2$fix_fct,model2$rnd_non)
  fct2<-map_dfc(fct_vec2,function(x)as.numeric(as.factor(phe[[x]])))
  fct_name2<-map(fct_vec2,function(x){
    name2<-levels(as.factor(phe[[x]]))
    data.frame(Name=name2,Code=1:length(name2))
  })
  names(fct_name2)<-fct_vec2
  names(fct2)<-fct_vec2
  
  ## 提取实数变量——协变量+性状表型值
  rl2<-select(phe,c(model2$fix_cov,model2$traits))
  
  ## 编码遗传随机效应
  if(need_pedigree){
    rnd_gen2<-map_dfc(model2$rnd_gen,function(x){
      ord[match(phe[[x]],ped[[1]])]
    })
    rnd_gen_name2<-data.frame(Name=ped[[1]],Code=ord)
    
    names(rnd_gen2)<-model2$rnd_gen
    
  }else{
    rnd_gen2<-map_dfc(model2$rnd_gen,function(x){
      ord[match(phe[[x]],fam_new[[1]])]
    })
    rnd_gen_name2<-data.frame(Name=fam_new[[1]],Code=ord)
    
    names(rnd_gen2)<-model2$rnd_gen
    
  }
  
  rnd_gen2[is.na(rnd_gen2)]<-0
  fct2[is.na(fct2)]<-0
  rl2[is.na(rl2)]<--999999
  
  # 开始进行计算
  int_name2<-c(model2$rnd_gen,model2$fix_fct,model2$rnd_non)
  int_code2<-c(match(model2$fix_fct,int_name2),
              match(model2$rnd_gen,int_name2),
              match(model2$rnd_non,int_name2))
  rnd_code2<-c(match(model2$rnd_gen,int_name2),
              match(model2$rnd_non,int_name2))
  
  ## 生成数据文件与DIR文件
  
  ### Define functions
  
  get_grm2<-function(grm_file="final.grm.gz",fam=fam_new){
    grm<-fread(grm_file)
    fam_id<-rnd_gen_name2[match(fam$V1,rnd_gen_name2$Name),'Code']
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
    fam_id<-rnd_gen_name2[match(fam$V1,rnd_gen_name2$Name),'Code']
    #fam_id<-1:nrow(fam)
    ginv_dat<-data.frame(c(0,fam_id[grm$V1]),c(0,fam_id[grm$V2]),
                         c(determinant(ginv)$modulus,ginv[upper.tri(ginv,diag = T)]))
    return(ginv_dat)
  }
  
  ### BLUP
  if(method=="BLUP"){
    walk(model2$traits,function(x){
      rl_tmp2<-c(model2$fix_cov,x)
      rl_dat2<-select(rl2,all_of(rl_tmp2))
      dmu_dat2<-cbind(rnd_gen2,fct2,rl_dat2)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      write.table(dmu_dat2,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 1 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1
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
",method,x,length(int_name2),length(rl_tmp2),
                  paste(int_name2,collapse = " "),
                  paste(rl_tmp2,collapse = " "),
                  length(rl_tmp2),length(int_code2),
                  paste(int_code2,collapse = " "),
                  length(rnd_code2),paste(1:length(rnd_code2),collapse = " "),
                  length(model2$fix_cov), paste(1:length(model2$fix_cov),collapse = " ")),
          file=paste0(outdir,"/DMU.DIR"))
    })
    
    ### ssGBLUP
  }else if(method=="ssGBLUP"){
    grm<-get_grm2()
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
1
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
    
    walk(model2$traits,function(x){
      rl_tmp2<-c(model2$fix_cov,x)
      rl_dat2<-select(rl2,all_of(rl_tmp2))
      dmu_dat2<-cbind(rnd_gen2,fct2,rl_dat2)
      outdir<-toupper(x)
      dir.create(outdir, showWarnings=F)
      dmu_dat2<-dmu_dat2[order(dmu_dat2[[1]]),]
      write.table(dmu_dat2,file=paste0(outdir,"/DAT"),row.names = F,col.names = F,quote = F)
      
      cat(sprintf("$COMMENT
%s for %s trait

$ANALYSE 1 31 0 0

$DATA  ASCII (%d,%d,-9999) DAT

$VARIABLE
%s
%s

$MODEL
1
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
",method,x,length(int_name2),length(rl_tmp2),
                  paste(int_name2,collapse = " "),
                  paste(rl_tmp2,collapse = " "),
                  length(rl_tmp2),length(int_code2),
                  paste(int_code2,collapse = " "),
                  length(rnd_code2),paste(1:length(rnd_code2),collapse = " "),
                  length(model2$fix_cov), paste(1:length(model2$fix_cov),collapse = " ")),
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
1
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
    grm<-get_grm2()
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
1
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
1
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
  
  dmu_log2<-map(model2$traits,function(x){
    outdir<-toupper(x)
    setwd(outdir)
    dmu1_log2<-system(paste0(DMU1," <DMU.DIR>DMU.lst"), intern=T)
    if(!str_detect(method,"VAR$")){
      dmuai_log2<-system(paste0(DMUAI," >>DMU.lst"), intern=T)
    }else{
      dmuai_log2<-system(paste0(DMU4," >>DMU.lst"), intern=T)
    }
    
    setwd("../")
    return(sum(attr(dmu1_log2,'status'),attr(dmuai_log2,'status')))
  })
  names(dmu_log2)<-model2$traits
  
  if(any(dmu_log2>0)){
    Success2=1
    error_trait2<-names(dmu_log2)[dmu_log2>0]
  }
  
  if(all(dmu_log2>0)){
    Success2=2
  }
  
  
  ## 先检验SOL文件是否生成，再检验模型是否收敛
  
  if(!str_detect(method,"VAR$")){
    converge_msg2<-map_dfr(model2$traits, function(x){
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
      converge_msg2<-""
      if(sol_exists & converged){
        converge_msg2<-paste0("Variance components estimation for ",x," has successfully converged.\n")
      }
      if(sol_exists & !converged){
        converge_msg2<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & converged){
        converge_msg2<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      if(!sol_exists & !converged){
        converge_msg2<-paste0("**[ERROR:]** Variance components estimation for ",x," does not converge. Please check your model, data and `",x,".dmu_run.lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg2,SOL=sol_exists, CON=converged)
      return(res)
    })
  }else{
    
    converge_msg2<-map_dfr(model2$traits, function(x){
      outdir<-toupper(x)
      sol_exists<-file.exists(paste0(outdir,"/SOL"))
      converge_msg2<-""
      converged<-FALSE
      if(sol_exists){
        converge_msg2<-paste0("Calculation for ",x," has successfully completed.\n")
        converged<-TRUE
      }else{
        converge_msg2<-paste0("**[ERROR:]** Solution for ",x," is not successfully created Please check your model, data and `",trait,".lst` file.\n")
      }
      res<-data.frame(TRT=x,MSG=converge_msg2,SOL=sol_exists, CON=converged)
      return(res)
    })
  }
  
  if(Success2!=2){
    ## 生成几个结果的表格
    ### Variance component
    
    if(!str_detect(method,"VAR$")){
      var_comp2<-as.data.frame(matrix(NA,length(model2$traits),length(rnd_code2)+2))
      colnames(var_comp2)<-c("Trait",int_name2[rnd_code2],"Residual")
      var_comp2$Trait<-model2$traits
      
      walk(model2$traits,function(x){
        outdir<-toupper(x)
        if(converge_msg2$SOL[converge_msg2$TRT==x] & converge_msg2$CON[converge_msg2$TRT==x]){
          var_com2<-read.table(paste0(outdir,"/PAROUT"),h=F)
          var_comp2[var_comp2$Trait==x,2:ncol(var_comp2)]<<-var_com2$V4
        }
      })
    }else{
      var_comp2<-cbind(model2$traits,var_comp_input)
      colnames(var_comp2)<-c("Trait",int_name2[rnd_code2],"Residual")
    }
    write.table(var_comp2,file="Variance_component2.txt",row.names = F,quote=F,sep="\t")
    
    get_estimate2<-function(trait){
      res<-list()
      res$fix_eff<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      res$genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)
      if(length(model2$rnd_non)>0){res$non_genetic_sol<-data.frame(Name=NA,Value=NA,Code=NA,Estimate=NA,SE=NA)}
      outdir<-toupper(trait)
      if(converge_msg2$SOL[converge_msg2$TRT==trait] & converge_msg2$CON[converge_msg2$TRT==trait]){
        ### 读取solution
        sol<-fread(paste0(outdir,"/SOL"))
        fix_sol<-sol[sol$V1==2,]
        gen_sol<-sol[sol$V1==4,]
        nongen_sol<-sol[sol$V1==3,]
        if(method=="GBLUP"){
        gen_sol<-sol[sol$V1==3,]
        nongen_sol<-sol[sol$V1==4,]}
        if(length(model2$fix_cov)>0){cov_sol<-sol[sol$V1==1,]}
        
        ### fixed factor effect
        fixed_fct_name<-int_name2[-rnd_code2]
        
        fix_fct_sol<-map_dfr(fixed_fct_name,function(x){
          code<-fct_name2[[x]]
          sol_dat<-fix_sol[fix_sol$V4==which(fixed_fct_name==x),]
          sol_dat2<-sol_dat[match(code$Code,sol_dat$V5),]
          sol_res<-data.frame(Name=x,Value=code$Name,Code=code$Code,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
          return(sol_res)
        })
        
        ### Covariate effect
        if(length(model2$fix_cov)>0){
          cov_name<-model2$fix_cov
          fix_cov_sol<-data.frame(Name=cov_name,Value=1,Code=1:length(cov_name),Estimate=cov_sol$V8,SE=cov_sol$V9)
        }else{
          fix_cov_sol<-NULL
        }
        
        ### Genetic effect
        if(need_pedigree){
          genetic_sol<-map_dfr(model2$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model2$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=ped[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }else{
          genetic_sol<-map_dfr(model2$rnd_gen,function(x){
            
            sol_dat<-gen_sol[gen_sol$V4==which(model2$rnd_gen==x),]
            sol_dat2<-sol_dat[match(ord,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fam_new[[1]],Code=ord,Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
        }
        
        ### Non-Genetic random effect
        if(length(model2$rnd_non)>0){
          non_genetic_sol<-map_dfr(model2$rnd_non,function(x){
            sol_dat<-nongen_sol[nongen_sol$V4==(length(model2$rnd_gen)+which(model2$rnd_non==x)),]
            sol_dat2<-sol_dat[match(fct_name2[[x]]$Code,sol_dat$V5),]
            sol_res<-data.frame(Name=x,Value=fct_name2[[x]]$Name,Code=fct_name2[[x]]$Code,
                                Estimate=sol_dat2$V8,SE=sol_dat2$V9)
            return(sol_res)
          })
          non_genetic_sol<-non_genetic_sol[complete.cases(non_genetic_sol),]
          
        }
        
        
        
        res$fix_eff<-rbind(fix_fct_sol,fix_cov_sol)
        res$genetic_sol<-genetic_sol
        if(length(model2$rnd_non)>0){res$non_genetic_sol<-non_genetic_sol}
      }
      return(res)
    }
    
    all_estimate2<-map(model2$traits,get_estimate2)
    names(all_estimate2)<-model2$traits
    
    fix_eff2<-rbindlist(map(model2$traits,function(x)all_estimate2[[x]]$fix_eff),idcol="Trait")
    fix_eff2$Trait<-model2$traits[fix_eff2$Trait]
    
    gen_eff2<-rbindlist(map(model2$traits,function(x)all_estimate2[[x]]$genetic_sol),idcol="Trait")
    gen_eff2$Trait<-model2$traits[gen_eff2$Trait]
    
    if(length(model2$rnd_non)>0){
      nongen_eff2<-rbindlist(map(model2$traits,function(x)all_estimate2[[x]]$non_genetic_sol),idcol="Trait")
      nongen_eff2$Trait<-model2$traits[nongen_eff2$Trait]
    }
    
    fwrite(fix_eff2, file="fixed_effect_estimate2.txt",row.names=F,quote=F,sep="\t")
    fwrite(gen_eff2, file="genetic_effect_estimate2.txt",row.names=F,quote=F,sep="\t")
    if(length(model2$rnd_non)>0){fwrite(nongen_eff2, file="nongenetic_random_effect_estimate2.txt",row.names=F,quote=F,sep="\t")}
  }
}else{
  Success2<-2
}
}

if(!is.null(model2$traits) & Success1==0 & Success2==0){
var_comp<-rbind(var_comp1,var_comp2)
converge_msg<-rbind(converge_msg1,converge_msg2)
fix_eff<-rbind(fix_eff1,fix_eff2)
gen_eff<-rbind(gen_eff1,gen_eff2)
}else{
var_comp<-var_comp1
converge_msg<-converge_msg1
fix_eff<-fix_eff1
gen_eff<-gen_eff1
}

if(!is.null(model2$traits) & Success1==0 & Success2==0 & length(model1$rnd_non)>0 & length(model2$rnd_non)>0){
nongen_eff<-rbind(nongen_eff1,nongen_eff2)
}else if(is.null(model2$traits) & Success1==0 & Success2==0 & length(model1$rnd_non)>0){
nongen_eff<-nongen_eff1
}

CPI_trait_check<-TRUE
trait<-unique(gen_eff[[1]])
if(is.null(model2$traits) & length(trait)==4 & sum(str_detect(trait,fixed("milk",ignore_case = T)))==1 & sum(str_detect(trait,fixed("fatpct",ignore_case = T)))==1 &sum(str_detect(trait,fixed("propct",ignore_case = T)))==1 & sum(str_detect(trait,fixed("scs",ignore_case = T)))==1){
ebv_all<-paste0(trait,"_EBV")
Value<-as.character(unique(gen_eff[[3]]))
Code<-as.numeric(unique(gen_eff[[4]]))
cpi_sol<-data.frame(ID=Value,Code=Code)
for(i in 1:length(trait)){
cpi_sol[ebv_all[i]]<-assign(ebv_all[i],gen_eff[[5]][gen_eff[[1]]==trait[i]])
}
cpi_sol["CPI2"]<-20*(cpi_sol[,str_detect(names(cpi_sol),fixed("milk",ignore_case = T))]*30/459+cpi_sol[,str_detect(names(cpi_sol),fixed("fatpct",ignore_case = T))]*15/0.16+cpi_sol[,str_detect(names(cpi_sol),fixed("propct",ignore_case = T))]*25/0.08-(cpi_sol[,str_detect(names(cpi_sol),fixed("scs",ignore_case = T))]-3)*10/0.16)
}else if(!is.null(model2$traits) & method=="BLUP" & length(trait)==7 & sum(str_detect(trait,fixed("milk",ignore_case = T)))==1 & sum(str_detect(trait,fixed("fatpct",ignore_case = T)))==1 &sum(str_detect(trait,fixed("propct",ignore_case = T)))==1 & sum(str_detect(trait,fixed("scs",ignore_case = T)))==1 & sum(str_detect(trait,fixed("type",ignore_case = T)))==1 & sum(str_detect(trait,fixed("ms",ignore_case = T)))==1 & sum(str_detect(trait,fixed("fl",ignore_case = T)))==1){
ebv_all<-paste0(trait,"_EBV")
Value<-as.character(unique(gen_eff[[3]]))
Code<-as.numeric(unique(gen_eff[[4]]))
cpi_sol<-data.frame(ID=Value,Code=Code)
for(i in 1:length(trait)){
cpi_sol[ebv_all[i]]<-assign(ebv_all[i],gen_eff[[5]][gen_eff[[1]]==trait[i]])
}
cpi_sol["CPI1"]<-20*(cpi_sol[,str_detect(names(cpi_sol),fixed("milk",ignore_case = T))]*30/459+cpi_sol[,str_detect(names(cpi_sol),fixed("fatpct",ignore_case = T))]*15/0.16+cpi_sol[,str_detect(names(cpi_sol),fixed("propct",ignore_case = T))]*25/0.08+cpi_sol[,str_detect(names(cpi_sol),fixed("type",ignore_case = T))]*5/5+cpi_sol[,str_detect(names(cpi_sol),fixed("ms",ignore_case = T))]*10/5+cpi_sol[,str_detect(names(cpi_sol),fixed("fl",ignore_case = T))]*5/5-(cpi_sol[,str_detect(names(cpi_sol),fixed("scs",ignore_case = T))]-3)*10/0.16)
cpi_sol["CPI3"]<-20*(cpi_sol[,str_detect(names(cpi_sol),fixed("milk",ignore_case = T))]*30/800+cpi_sol[,str_detect(names(cpi_sol),fixed("fatpct",ignore_case = T))]*10/0.3+cpi_sol[,str_detect(names(cpi_sol),fixed("propct",ignore_case = T))]*20/0.12+cpi_sol[,str_detect(names(cpi_sol),fixed("type",ignore_case = T))]*5/5+cpi_sol[,str_detect(names(cpi_sol),fixed("ms",ignore_case = T))]*15/5+cpi_sol[,str_detect(names(cpi_sol),fixed("fl",ignore_case = T))]*10/5-(cpi_sol[,str_detect(names(cpi_sol),fixed("scs",ignore_case = T))]-3)*10/0.46)
}else if(!is.null(model2$traits) & method=="GBLUP" & length(trait)==7 & sum(str_detect(trait,fixed("milk",ignore_case = T)))==1 & sum(str_detect(trait,fixed("fatpct",ignore_case = T)))==1 &sum(str_detect(trait,fixed("propct",ignore_case = T)))==1 & sum(str_detect(trait,fixed("scs",ignore_case = T)))==1 & sum(str_detect(trait,fixed("type",ignore_case = T)))==1 & sum(str_detect(trait,fixed("ms",ignore_case = T)))==1 & sum(str_detect(trait,fixed("fl",ignore_case = T)))==1){
ebv_all<-paste0(trait,"_GEBV")
Value<-as.character(unique(gen_eff[[3]]))
Code<-as.numeric(unique(gen_eff[[4]]))
cpi_sol<-data.frame(ID=Value,Code=Code)
for(i in 1:length(trait)){
cpi_sol[ebv_all[i]]<-assign(ebv_all[i],gen_eff[[5]][gen_eff[[1]]==trait[i]])
}
cpi_sol["GCPI"]<-20*(cpi_sol[,str_detect(names(cpi_sol),fixed("milk",ignore_case = T))]*30/800+cpi_sol[,str_detect(names(cpi_sol),fixed("fatpct",ignore_case = T))]*15/0.3+cpi_sol[,str_detect(names(cpi_sol),fixed("propct",ignore_case = T))]*25/0.12+cpi_sol[,str_detect(names(cpi_sol),fixed("type",ignore_case = T))]*5/5+cpi_sol[,str_detect(names(cpi_sol),fixed("ms",ignore_case = T))]*10/5+cpi_sol[,str_detect(names(cpi_sol),fixed("fl",ignore_case = T))]*5/5-(cpi_sol[,str_detect(names(cpi_sol),fixed("scs",ignore_case = T))]-3)*10/0.46)+80
}else{
CPI_trait_check<-FALSE
}

if(CPI_trait_check==TRUE){
write.table(cpi_sol,file="CPI_calculation.txt",row.names = F,quote=F,sep="\t")
}


# 输出html格式报告

con<-file(paste0(report_path,"/report_template1.R"), open = "r")
con_out<-file("report.R", open = "w")
model_formula1<-paste(unlist(lapply(model1$traits,function(x)paste0("\\\\text{",x,"}"))),
                     paste0("=",paste0("\\\\text{",paste(c(model1$fix_fct,model1$xp,model1$rnd_gen,model1$rnd_non),
                                                         collapse = "} + \\\\text{"),"}")),collapse = "\\\\\\\\")
if(!is.null(model2$traits)){
model_formula2<-paste(unlist(lapply(model2$traits,function(x)paste0("\\\\text{",x,"}"))),
                     paste0("=",paste0("\\\\text{",paste(c(model2$fix_fct,model2$xp,model2$rnd_gen,model2$rnd_non),
                                                         collapse = "} + \\\\text{"),"}")),collapse = "\\\\\\\\")
model_formula<-paste0("The model is as follows:\n#' $$\n#' ",model_formula1,model_formula2,"\n#' $$")
}else{
model_formula<-paste0("The model is as follows:\n#' $$\n#' ",model_formula1,"\n#' $$")
}

data_summary_entry<-"The basic information of your input data is as follows:\n#+ echo=FALSE\nknitr::kable(data_summary)"

if(length(model1$rnd_non)>0){
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

if(CPI_trait_check==TRUE){
  cpi_entry<-"**Note that you can find all CPI results in supplementary files.**
#' 
#' ### CPI estimates of first 50 records
#' 
#+ echo=FALSE
cpi_sol %>% 
  do(head(.,50)) %>%
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
  
  if(data_summary_check & phenotype_check_pass & genotype_check_pass & method!="GBLUP"){
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
  
  if(model_check_pass & data_summary_check & phenotype_check_pass & genotype_check_pass & pedigree_check_pass & var_check_pass & Success1!=2 & Success2!=2){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, the calculation is not run and no effect estimate is shown. Please check your data first!\\*\\*", 
                         model_res_entry)
  }

  if(CPI_trait_check<-TRUE){
    oneLine<-str_replace(oneLine, "\\*\\*\\[ERROR: \\] Something is wrong, maybe the trait names in Phenotype file are inconformity with naming rules. Please check the trait number and their symbol first!\\*\\*", 
                         cpi_entry)
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
res_files<-map_chr(model1$traits,function(x){
  system(paste0("cp ",toupper(x),"/DMU.lst", " ./",x,".dmu_run.lst"))
  return(paste0(x,".dmu_run.lst"))
})

if(!is.null(model2$traits)){
res_files2<-map_chr(model2$traits,function(x){
  system(paste0("cp ",toupper(x),"/DMU.lst", " ./",x,".dmu_run.lst"))
  return(paste0(x,".dmu_run.lst"))
})
res_files<-c(res_files,res_files2)
}

if(Success1!=2 & Success2!=2 & is.null(model2$traits) & CPI_trait_check==FALSE){
  if(length(model1$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt")
  }else{
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
             "Variance_component1.txt")
  }
}else if(Success1!=2 & Success2!=2 & !is.null(model2$traits) & CPI_trait_check==FALSE){
  if(length(model1$rnd_non)>0 & length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt",
                 "genetic_effect_estimate2.txt","fixed_effect_estimate2.txt",
                 "Variance_component2.txt","nongenetic_random_effect_estimate2.txt")
  }else if(length(model1$rnd_non)>0 & !length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt",
                  "genetic_effect_estimate2.txt","fixed_effect_estimate2.txt",
                 "Variance_component2.txt")
  }else if(!length(model1$rnd_non)>0 & length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","genetic_effect_estimate2.txt",
                 "fixed_effect_estimate2.txt",
                 "Variance_component2.txt","nongenetic_random_effect_estimate2.txt")
  }else if(!length(model1$rnd_non)>0 & !length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","genetic_effect_estimate2.txt",
                 "fixed_effect_estimate2.txt","Variance_component2.txt")
  }
}else if(Success1!=2 & Success2!=2 & is.null(model2$traits) & CPI_trait_check==TRUE){
  if(length(model1$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt","CPI_calculation.txt")
  }else{
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
             "Variance_component1.txt","CPI_calculation.txt")
  }
}else if(Success1!=2 & Success2!=2 & !is.null(model2$traits) & CPI_trait_check==TRUE){
  if(length(model1$rnd_non)>0 & length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt",
                 "genetic_effect_estimate2.txt","fixed_effect_estimate2.txt",
                 "Variance_component2.txt","nongenetic_random_effect_estimate2.txt","CPI_calculation.txt")
  }else if(length(model1$rnd_non)>0 & !length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","nongenetic_random_effect_estimate1.txt",
                  "genetic_effect_estimate2.txt","fixed_effect_estimate2.txt",
                 "Variance_component2.txt","CPI_calculation.txt")
  }else if(!length(model1$rnd_non)>0 & length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","genetic_effect_estimate2.txt",
                 "fixed_effect_estimate2.txt",
                 "Variance_component2.txt","nongenetic_random_effect_estimate2.txt","CPI_calculation.txt")
  }else if(!length(model1$rnd_non)>0 & !length(model2$rnd_non)>0){
    res_files<-c(res_files,"genetic_effect_estimate1.txt","fixed_effect_estimate1.txt",
                 "Variance_component1.txt","genetic_effect_estimate2.txt",
                 "fixed_effect_estimate2.txt","Variance_component2.txt","CPI_calculation.txt")
  }
}


system(paste0("zip hegs_res.zip ", paste(res_files,collapse = " ")))

# 发送邮件
Success<-Success1+Success2
Success<-as.character(Success)
mail_text<-switch(Success,
                  "0"="Dear user,\n\nYour run of HEGS has sucessfully finished. Please find data analysis report and results in the attachment. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "1"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "2"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "3"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team",
                  "4"="Dear user,\n\nYour run of HEGS has run into some errors. Please find data analysis report in the attachment and check your data and model once more. If you have any question, please do not hesitate to contact zhe_zhang@zju.edu.cn.\n\nBest\n\nHEGS Team")

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












