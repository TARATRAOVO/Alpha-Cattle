import pandas as pd

# 读取CSV文件
df = pd.read_csv('liftover/variant_coordinates.csv')

# 筛选条件:
# 1. Reference Sequence 为 ARS-UCD1.2
# 2. mutation列为形如 X>Y 的格式(表示单点突变)
filtered_df = df[
    (df['Reference Sequence'] == 'ARS-UCD1.2') & 
    (df['mutation'].str.contains('>',na=False))
]

# 选择需要的列
result_df = filtered_df

# 保存结果到新的CSV文件
result_df.to_csv('single_nucleotide_variants_1.2.csv', index=False)

# 打印找到的记录数
print(f"找到 {len(result_df)} 条符合条件的记录") 