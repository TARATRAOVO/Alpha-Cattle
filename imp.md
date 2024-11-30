以下是对`Impute_job_submit.R`文件中基因型填补步骤的详细解释：

```Impute_job_submit.R
system(paste0("./plink/plink --bfile ",
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
system(paste0("./plink/plink --vcf ",geno_names[i],".beagle.WGS.vcf.gz --make-bed --out ",geno_names[i],".beagle.WGS"))
setwd(path)
```

### 步骤详解

1. **将PLINK二进制文件转换为VCF格式**

    ```R
    system(paste0("./plink/plink --bfile ",
                  geno_names[i],
                  " --export vcf  --out ",geno_names[i]))
    ```

    - 使用PLINK工具将基因型数据（`.bed`, `.bim`, `.fam`文件）转换为VCF格式，便于后续处理。
    - 参数解释：
        - `--bfile`: 输入的PLINK二进制文件前缀。
        - `--export vcf`: 指定导出格式为VCF。
        - `--out`: 输出文件的前缀。

2. **压缩VCF文件**

    ```R
    system(paste0("bgzip ",geno_names[i],".vcf"))
    ```

    - 使用`bgzip`工具对生成的VCF文件进行压缩，生成`.vcf.gz`文件，节省存储空间并加快后续处理速度。

3. **为压缩后的VCF文件创建索引**

    ```R
    system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".vcf.gz"))
    ```

    - 使用`bcftools`为压缩的VCF文件创建索引文件（`.csi`或`.tbi`），以便快速随机访问。

4. **提取染色体名称**

    ```R
    ref_names[i]=str_extract(geno_names[i], "chr\\d")
    ```

    - 使用正则表达式从`geno_names[i]`中提取染色体名称（如`chr1`, `chr2`等），用于后续与参考基因型数据对齐。

5. **对齐基因型数据与参考基因型数据**

    ```R
    system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform.150K"))
    ```

    - 使用`conform-gt`工具将样本的VCF文件与参考基因型数据对齐。
    - 参数解释：
        - `ref`: 参考基因型VCF文件的路径。
        - `gt`: 样本的基因型VCF文件。
        - `strict=false`: 非严格模式，允许部分不匹配。
        - `match=POS`: 按位置匹配SNP。
        - `chrom`: 染色体编号。
        - `out`: 输出对齐后的VCF文件前缀。

6. **为对齐后的VCF文件创建索引**

    ```R
    system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.150K.vcf.gz"))
    ```

    - 同样使用`bcftools`为对齐后的VCF文件创建索引文件，便于后续快速访问。

7. **进行基因型填补（BEAGLE）**

    ```R
    system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=../",breed,"_ref/",ref_names[i],".vcf.gz gt=",geno_names[i],".conform.150K.vcf.gz ne=100 chrom=",i," out=",geno_names[i],".beagle.150K"))
    ```

    - 使用BEAGLE工具进行基因型填补（Imputation）。
    - 参数解释：
        - `ref`: 参考基因型VCF文件。
        - `gt`: 对齐后的样本基因型VCF文件。
        - `ne=100`: 突变重组事件的数目，影响填补精度和计算时间。
        - `chrom`: 染色体编号。
        - `out`: 输出填补后的VCF文件前缀。

8. **为填补后的VCF文件创建索引**

    ```R
    system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".beagle.150K.vcf.gz"))
    ```

    - 使用`bcftools`为填补后的VCF文件创建索引。

9. **再次对齐填补后的基因型数据与更大范围的参考基因型数据（Whole Genome）**

    ```R
    system(paste0("java -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/conform-gt.24May16.cee.jar ref=/public/home/chenjc/ltw/nyjh/fill/imputed1000/genotype_beagle_",ref_names[i],".vcf gt=",geno_names[i],".beagle.150K.vcf.gz strict=false match=POS chrom=",i," out=",geno_names[i],".conform.WGS"))
    ```

    - 使用`conform-gt`工具将已经填补的基因型数据与一个更大范围的参考基因型数据（如Whole Genome Sequence）对齐，以进一步提高填补的准确性和覆盖范围。

10. **为再次对齐后的VCF文件创建索引**

    ```R
    system(paste0("/srv/bcftools-1.14/bcftools index ",geno_names[i],".conform.WGS.vcf.gz"))
    ```

    - 使用`bcftools`为再次对齐后的VCF文件创建索引。

11. **再次进行基因型填补（BEAGLE）**

    ```R
    system(paste0("java -Xmx60g -jar /public/home/chenjc/ltw/nyjh/fill/softwareandImputeacc/beagle.18May20.d20.jar ref=/public/home/chenjc/ltw/nyjh/fill/imputed1000/genotype_beagle_",ref_names[i],".vcf gt=",geno_names[i],".conform.WGS.vcf.gz ne=100 chrom=",i," out=",geno_names[i],".beagle.WGS"))
    ```

    - 使用BEAGLE工具基于更大范围的参考基因型数据进行第二轮基因型填补。
    - **注意**：此步骤分配了`-Xmx60g`，即允许Java虚拟机使用60GB内存，以处理大规模的数据。

12. **将最终填补后的VCF文件转换回PLINK二进制格式**

    ```R
    system(paste0("./plink/plink --vcf ",geno_names[i],".beagle.WGS.vcf.gz --make-bed --out ",geno_names[i],".beagle.WGS"))
    ```

    - 使用PLINK将最终填补后的VCF文件转换回PLINK二进制格式（`.bed`, `.bim`, `.fam`），便于后续的基因型分析和关联研究。

13. **恢复工作目录**

    ```R
    setwd(path)
    ```

    - 将当前工作目录设置回初始路径，以确保后续操作的路径正确。

### 总结

整个填补过程主要包括以下几个关键步骤：

1. **格式转换**：将PLINK格式转换为VCF格式，方便与参考基因型数据对齐和填补工具的使用。
2. **数据压缩与索引**：使用`bgzip`和`bcftools`对VCF文件进行压缩和索引，提高数据处理效率。
3. **基因型对齐**：使用`conform-gt`工具确保样本基因型数据与参考基因型数据在结构和位点上的一致性。
4. **基因型填补**：使用BEAGLE工具基于参考数据对缺失或未测位点进行填补，提高基因型数据的完整性和准确性。
5. **格式还原**：将填补后的VCF文件转换回PLINK二进制格式，便于后续分析。
6. **循环处理**：该过程在基因型文件列表`geno_names`中的每个文件上循环执行，以批量处理多个样本或染色体的数据。

通过上述步骤，脚本实现了对基因型数据的高效填补，确保数据的完整性和质量，为后续的遗传分析提供可靠的数据基础。