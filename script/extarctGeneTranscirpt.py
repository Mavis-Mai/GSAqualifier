##希望通过输入的注释文件，然后生成一个transcript-GeneID这样的映射关系，并保存到文件内。
##我需要一个python流程
#!/usr/bin/env python3
import sys
import re
import pandas as pd
from collections import defaultdict
import argparse

def extract_mappings(gxf_file, output_prefix, transcript_pattern=None, gene_pattern=None):
    """
    从GFF/GTF文件中提取transcript-to-GeneID的映射关系
    
    参数:
        gxf_file: 输入GFF3/GTF文件路径
        output_prefix: 输出文件前缀
        transcript_pattern: 转录本ID提取模式(可选)
        gene_pattern: 基因ID提取模式(可选)
        
    返回:
        DataFrame: 包含映射关系的数据框
        dict: 统计信息
    """
    # 设置默认的正则表达式模式
    if transcript_pattern is None:
        transcript_pattern = re.compile(r'transcript_id[=\s]*["\']?([A-Za-z0-9_.\-:]+)["\']?[;\s]?')
    if gene_pattern is None:
        gene_pattern = re.compile(r'gene_id[=\s]*["\']?([A-Za-z0-9_.\-:]+)["\']?[;\s]?')
    
    mappings = defaultdict(dict)
    feature_counts = defaultdict(int)
    
    with open(gxf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) < 9:
                continue
                
            feature_type = fields[2].lower()
            attrs = fields[8]
            
            # 对转录本和基因特征特别注意
            if feature_type in ['transcript', 'mrna', 'gene']:
                try:
                    # 尝试提取转录本ID和基因ID
                    transcript_id = transcript_pattern.search(attrs).group(1) if transcript_pattern else None
                    gene_id = gene_pattern.search(attrs).group(1) if gene_pattern else None
                    
                    if transcript_id and gene_id:
                        mappings[transcript_id]['GeneID'] = gene_id
                        mappings[transcript_id]['Chr'] = fields[0]
                        mappings[transcript_id]['Strand'] = fields[6]
                        mappings[transcript_id]['FeatureType'] = feature_type
                    elif gene_id and not transcript_id:
                        # 基因级别的条目
                        pass
                        
                    feature_counts[feature_type] += 1
                except (AttributeError, IndexError) as e:
                    continue
    
    # 创建数据框
    df = pd.DataFrame.from_dict(mappings, orient='index')
    df.index.name = 'TranscriptID'
    df = df.reset_index()

    
    # 保存结果
    output_file = f"{output_prefix}_transcript_gene_mapping.tsv"
    df.to_csv(output_file, sep='\t', index=False)
    
    return mappings


def main():
    """命令行界面"""
    parser = argparse.ArgumentParser(
        description='从GFF/GTF文件提取transcript-to-GeneID映射关系',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('input', help='输入GFF/GTF文件')
    parser.add_argument('-o', '--output', default='transcript_gene',
                       help='输出文件前缀(不包含扩展名)')
    parser.add_argument('--transcript_pattern', default=None,
                       help='自定义转录本ID提取模式(正则表达式)')
    parser.add_argument('--gene_pattern', default=None,
                       help='自定义基因ID提取模式(正则表达式)')
    
    args = parser.parse_args()
    
    print(f"\n从 {args.input} 提取转录本-基因映射关系...")
    
    # 编译正则表达式
    transcript_re = re.compile(args.transcript_pattern) if args.transcript_pattern else None
    gene_re = re.compile(args.gene_pattern) if args.gene_pattern else None
    
    # 执行提取
    df, stats = extract_mappings(
        args.input,
        args.output,
        transcript_pattern=transcript_re,
        gene_pattern=gene_re
    )
    
    # 打印摘要信息
    print("\n处理完成! 结果摘要:")
    print(f"- 共找到 {stats['total_transcripts']} 个转录本")
    print(f"- 其中 {stats['transcripts_with_gene']} 个有基因ID映射")
    print(f"- 输出文件已保存为: {args.output}_transcript_gene_mapping.tsv")
    print(f"- 统计信息保存在: {args.output}_mapping_stats.txt")
    
    # 显示前几行
    print("\n输出文件前5行:")
    print(df.head())


if __name__ == "__main__":
    main()
