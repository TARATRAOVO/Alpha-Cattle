import os
import glob
import subprocess
from multiprocessing import Pool
from pathlib import Path
from tqdm import tqdm

def create_index(vcf_file):
    """为单个VCF文件创建索引"""
    try:
        # 检查是否已经存在索引文件
        if os.path.exists(f"{vcf_file}.tbi"):
            print(f"索引已存在: {vcf_file}.tbi")
            return True
            
        # 执行tabix命令创建索引
        cmd = f"tabix -p vcf {vcf_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"成功创建索引: {vcf_file}")
            return True
        else:
            print(f"创建索引失败: {vcf_file}")
            print(f"错误信息: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"处理文件时出错 {vcf_file}: {str(e)}")
        return False

def main():
    # 设置VCF文件所在目录
    vcf_dir = "/srv/shiny-server/data_base/genomic_data"
    
    # 获取所有.vcf.gz文件
    vcf_files = glob.glob(os.path.join(vcf_dir, "*.vcf.gz"))
    
    if not vcf_files:
        print(f"在 {vcf_dir} 中没有找到.vcf.gz文件")
        return
        
    print(f"找到 {len(vcf_files)} 个VCF文件")
    
    # 创建进程池，使用CPU核心数的进程
    with Pool(processes=10) as pool:
        results = list(tqdm(
            pool.imap(create_index, vcf_files),
            total=len(vcf_files),
            desc="创建索引"
        ))
    
    # 统计结果
    success = sum(1 for r in results if r)
    failed = len(results) - success
    
    print("\n索引创建完成:")
    print(f"成功: {success}")
    print(f"失败: {failed}")

if __name__ == "__main__":
    main() 