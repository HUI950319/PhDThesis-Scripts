import argparse
import subprocess
import sys
import os
from pathlib import Path
import pandas as pd
from pyscenic.cli.utils import load_signatures


def get_cistarget_files(species):
    """
    根据物种获取对应的cisTarget文件路径
    
    Args:
        species (str): 物种名称 ('human' 或 'mouse')
        
    Returns:
        dict: 包含各种文件路径的字典
    """
    base_dir = Path("cistarget")
    
    if species.lower() == "human":
        return {
            "tfs": base_dir / "hsa_hgnc_tfs.motifs-v10.txt",
            "motif_path": base_dir / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
            "db_500bp": base_dir / "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "db_10kb": base_dir / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        }
    elif species.lower() == "mouse":
        return {
            "tfs": base_dir / "mmu_mgi_tfs.motifs-v10.txt",
            "motif_path": base_dir / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
            "db_500bp": base_dir / "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "db_10kb": base_dir / "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
        }
    else:
        raise ValueError(f"不支持的物种: {species}。支持: human, mouse")


def validate_files(cistarget_files):
    """
    验证cisTarget文件是否存在
    
    Args:
        cistarget_files (dict): cisTarget文件路径字典
        
    Raises:
        FileNotFoundError: 当文件不存在时抛出异常
    """
    for file_type, file_path in cistarget_files.items():
        if not file_path.exists():
            raise FileNotFoundError(f"cisTarget文件不存在: {file_path}")


def get_motif_logo(regulon):
    """
    获取motif logo的URL
    
    Args:
        regulon: regulon对象
        
    Returns:
        str: motif logo的URL
    """
    base_url = "http://motifcollections.aertslab.org/v10nr_clust/logos/"
    for elem in regulon.context:
        if elem.endswith('.png'):
            return base_url + elem
    return ""


def regulon_to_gmt(regulon_file, project_prefix, min_regulon_size=10):
    """
    将regulon文件转换为GMT格式
    
    Args:
        regulon_file (str): regulon文件路径
        project_prefix (str): 输出文件前缀
        min_regulon_size (int): 最小regulon基因数，用于过滤regulon
        
    Returns:
        tuple: (gmt_file_path, txt_file_path) 输出文件路径
    """
    regulons = load_signatures(regulon_file)
    select_cols = [i.name for i in regulons if len(i.genes) >= min_regulon_size]
    
    gmt_file = f"{project_prefix}.gmt"
    txt_file = f"{project_prefix}.txt"
    
    with open(gmt_file, 'w') as fo1, open(txt_file, 'w') as fo2:
        for i in regulons:
            if i.name in select_cols:
                motif = get_motif_logo(i)
                genes = "\t".join(i.genes)
                tf = f"{i.transcription_factor}({len(i.genes)}g)"
                fo1.write(f"{tf}\t{motif}\t{genes}\n")
                fo2.write(f"{tf}\t{motif}\t{genes.replace(chr(9), ',')}\n")
    
    print(f"GMT文件已生成: {gmt_file}")
    print(f"TXT文件已生成: {txt_file}")
    
    return gmt_file, txt_file


def main():
    parser = argparse.ArgumentParser(description='Run pySCENIC pipeline with external parameters.')
    parser.add_argument('--input_csv', required=True, help='Input expression matrix (CSV file)')
    parser.add_argument('--output_path', required=True, help='Output directory path')
    parser.add_argument('--species', required=True, choices=['human', 'mouse'], 
                       help='Species (human or mouse)')
    parser.add_argument('--num_workers', type=int, default=10, help='Number of workers (default: 10)')
    parser.add_argument('--seed', type=int, default=777, help='Random seed (default: 777)')
    parser.add_argument('--min_regulon_size', type=int, default=10, help='Minimum regulon size for filtering (default: 10)')
    args = parser.parse_args()

    # 获取cisTarget文件路径
    try:
        cistarget_files = get_cistarget_files(args.species)
        validate_files(cistarget_files)
        print(f"使用 {args.species} 物种的cisTarget文件:")
        for file_type, file_path in cistarget_files.items():
            print(f"  {file_type}: {file_path}")
    except Exception as e:
        print(f"cisTarget文件错误: {e}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    output_dir = Path(args.output_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 定义输出文件名
    grn_output = output_dir / "grn.tsv"
    ctx_output = output_dir / "ctx.tsv"

    print(f"输出文件:")
    print(f"  GRN: {grn_output}")
    print(f"  CTX: {ctx_output}")

    # Step 1: Build GRN
    if grn_output.exists():
        print(f"GRN文件已存在，跳过GRN构建步骤: {grn_output}")
    else:
        grn_cmd = [
            'arboreto_with_multiprocessing.py',
            args.input_csv,
            str(cistarget_files["tfs"]),
            '--method', 'grnboost2',
            '--output', str(grn_output),
            '--num_workers', str(args.num_workers),
            '--seed', str(args.seed)
        ]
        print('Running:', ' '.join(grn_cmd))
        ret = subprocess.run(grn_cmd, stdout=sys.stdout, stderr=sys.stderr, bufsize=1, universal_newlines=True)
        if ret.returncode != 0:
            print('Error in arboreto_with_multiprocessing.py step', file=sys.stderr)
            sys.exit(1)

    # Step 2: cisTarget
    if ctx_output.exists():
        print(f"CTX文件已存在，跳过cisTarget步骤: {ctx_output}")
    else:
        ctx_cmd = [
            'pyscenic', 'ctx',
            str(grn_output),
            str(cistarget_files["db_500bp"]), str(cistarget_files["db_10kb"]),
            '--annotations_fname', str(cistarget_files["motif_path"]),
            '--expression_mtx_fname', args.input_csv,
            '--output', str(ctx_output),
            '--num_workers', str(args.num_workers)
        ]
        print('Running:', ' '.join(ctx_cmd))
        ret = subprocess.run(ctx_cmd, stdout=sys.stdout, stderr=sys.stderr, bufsize=1, universal_newlines=True)
        if ret.returncode != 0:
            print('Error in pyscenic ctx step', file=sys.stderr)
            sys.exit(1)

    # Step 3: Convert regulons to GMT format
    try:
        gmt_file, txt_file = regulon_to_gmt(
            str(ctx_output), 
            str(output_dir / "regulons"), 
            args.min_regulon_size
        )
        print(f"Regulon转换完成!")
        print(f"GMT文件: {gmt_file}")
        print(f"TXT文件: {txt_file}")
    except Exception as e:
        print(f"Regulon转换错误: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"pySCENIC分析完成!")
    print(f"结果保存在: {output_dir}")

if __name__ == '__main__':
    main()
