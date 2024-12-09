import pandas as pd
import re

# 读取CSV文件
df = pd.read_csv('variants.csv')

def extract_positions(coord):
    if pd.isna(coord):
        return pd.NA, pd.NA, pd.NA
    
    coord = str(coord)
    
    # 处理包含参考序列号的格式，如 NC_037328.1:g.81082187C>T
    if ':g.' in coord:
        coord = coord.split(':g.')[1]
    
    # 处理单点替换，如 g.123A>G 或 81082187C>T
    if '>' in coord:
        pos = re.search(r'(\d+)([ATCG])>([ATCG])', coord)
        if pos:
            start = int(pos.group(1))
            mutation = f"{pos.group(2)}>{pos.group(3)}"
            return start, start + 1, mutation
    
    # 处理删除，如 g.123_456del 或 g.123del
    elif 'del' in coord:
        positions = re.search(r'(\d+)(?:_(\d+))?del([ATCG]*)', coord)
        if positions:
            start = int(positions.group(1))
            end = int(positions.group(2)) + 1 if positions.group(2) else start + 1
            deleted_seq = positions.group(3) if positions.group(3) else 'del'
            return start, end, deleted_seq
    
    # 处理插入，如 g.123_124insA 或 g.123insA
    elif 'ins' in coord:
        positions = re.search(r'(\d+)(?:_(\d+))?ins([ATCG]*)', coord)
        if positions:
            start = int(positions.group(1))
            end = int(positions.group(2)) if positions.group(2) else start + 1
            inserted_seq = positions.group(3) if positions.group(3) else 'ins'
            return start, end, inserted_seq
    
    # 处理重复，如 g.123_456dup
    elif 'dup' in coord:
        positions = re.search(r'(\d+)(?:_(\d+))?dup', coord)
        if positions:
            start = int(positions.group(1))
            end = int(positions.group(2)) + 1 if positions.group(2) else start + 1
            return start, end, 'dup'
            
    return pd.NA, pd.NA, pd.NA

# 添加变异类型说明
variant_types = {
    'missense': '单个核苷酸改变导致氨基酸改变',
    'nonsense': '单个核苷酸改变导致提前终止密码子',
    'splicing': '影响RNA剪接位点的变异',
    'deletion': '序列缺失',
    'insertion': '序列插入',
    'duplication': '序列重复',
    'delins': '序列缺失并插入',
}

# 提取所需信息
variant_info = df[[
    'Variant Phenotype',     # 变异表型
    'g. or m.',              # 基因组坐标
    'Type of Variant',       # 变异类型
    'Gene',                  # 基因名称
    'Chr.',                   # 染色体位置
    'Reference Sequence'      # 参考基因组
]].copy()

# 添加变异类型说明列
variant_info['variant_description'] = variant_info['Type of Variant'].map(variant_types)

# 提取开始和结束位置以及变异形式
variant_info[['start', 'end', 'mutation']] = variant_info['g. or m.'].apply(lambda x: pd.Series(extract_positions(x)))

# 去除空值行
variant_info = variant_info.dropna(subset=['start', 'end'])

# 将位置转换为整数
variant_info['start'] = variant_info['start'].astype(int)
variant_info['end'] = variant_info['end'].astype(int)

# 按染色体位置和起始位置排序
variant_info = variant_info.sort_values(['Chr.', 'start'])

# 导出为CSV文件
variant_info.to_csv('variant_coordinates.csv', index=False)