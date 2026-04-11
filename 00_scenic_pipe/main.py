#!/usr/bin/env python3
"""
SCENIC Pipeline Main Script
实现GRN (Gene Regulatory Network) 分析的主控制脚本
"""

import os
import sys
import subprocess
import configparser
import argparse
from pathlib import Path


def read_config(config_file):
    """
    读取配置文件
    
    Args:
        config_file (str): 配置文件路径
        
    Returns:
        configparser.ConfigParser: 配置对象
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"配置文件不存在: {config_file}")
    
    config = configparser.ConfigParser()
    config.read(config_file, encoding='utf-8')
    return config


def validate_config(config):
    """
    验证配置文件的完整性
    
    Args:
        config (configparser.ConfigParser): 配置对象
        
    Raises:
        ValueError: 当配置不完整时抛出异常
    """
    required_sections = ['global', 'preprocess', 'pyscenic', 'aucell']
    for section in required_sections:
        if section not in config.sections():
            raise ValueError(f"配置文件缺少必需的节: [{section}]")
    
    # 验证global节
    global_required = ['species', 'output_path']
    for param in global_required:
        if not config.has_option('global', param):
            raise ValueError(f"global节缺少必需参数: {param}")
    
    # 验证物种设置
    species = config.get('global', 'species')
    if species not in ['human', 'mouse']:
        raise ValueError(f"不支持的物种: {species}。支持: human, mouse")
    
    # 验证preprocess节
    preprocess_required = ['input_file']
    for param in preprocess_required:
        if not config.has_option('preprocess', param):
            raise ValueError(f"preprocess节缺少必需参数: {param}")
    
    # 检查Rscript路径
    rscript_path = config.get('preprocess', 'rscript_path', fallback='Rscript')
    if not os.path.exists(rscript_path):
        print(f"警告: Rscript路径不存在: {rscript_path}")
        print("将尝试使用系统默认的Rscript")
        rscript_path = 'Rscript'
    
    # 验证Rscript是否可用
    try:
        subprocess.run([rscript_path, "--version"], check=True, capture_output=True)
        print(f"✓ Rscript可用: {rscript_path}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"✗ Rscript不可用: {rscript_path}")
        print("请检查Rscript路径是否正确")
        sys.exit(1)
    
    # 验证pyscenic节 - 现在只需要可选参数
    # 所有必需参数都由脚本自动处理
    
    # 验证aucell节 - 验证可选参数
    method = config.get('aucell', 'method', fallback='AUCell')
    if method not in ['AUCell', 'UCell']:
        raise ValueError(f"不支持的AUCell方法: {method}。支持: AUCell, UCell")


def run_preprocess(config):
    """
    运行预处理脚本 (preprocess.R)
    
    Args:
        config (configparser.ConfigParser): 配置对象
        
    Returns:
        bool: 成功返回True，失败返回False
    """
    print("=" * 60)
    print("步骤 1: 运行预处理脚本 (preprocess.R)")
    print("=" * 60)
    
    # 获取Rscript路径
    rscript_path = config.get('preprocess', 'rscript_path', fallback='Rscript')
    
    # 构建R脚本命令
    r_script = "preprocess/preprocess.R"
    cmd = [
        rscript_path, r_script,
        "--input_file", config.get('preprocess', 'input_file'),
        "--output_path", config.get('global', 'output_path'),
        "--species", config.get('global', 'species'),
        "--n_cells", config.get('preprocess', 'n_cells'),
        "--threads", config.get('preprocess', 'threads'),
        "--k", config.get('preprocess', 'k'),
        "--min_cells_per_gene", config.get('preprocess', 'min_cells_per_gene')
    ]
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr, bufsize=1, universal_newlines=True)
        print("预处理脚本执行成功!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"预处理脚本执行失败!")
        print(f"错误代码: {e.returncode}")
        return False


def run_pyscenic(config):
    """
    运行pySCENIC脚本 (run_pyscenic.py)
    
    Args:
        config (configparser.ConfigParser): 配置对象
        
    Returns:
        bool: 成功返回True，失败返回False
    """
    print("=" * 60)
    print("步骤 2: 运行pySCENIC脚本 (run_pyscenic.py)")
    print("=" * 60)
    
    # 构建输入文件路径
    input_csv = os.path.join(config.get('global', 'output_path'), 'imputed.mat.csv')
    
    # 检查输入文件是否存在
    if not os.path.exists(input_csv):
        print(f"错误: 输入文件不存在: {input_csv}")
        print("请先运行预处理步骤或确保文件存在")
        return False
    
    # 构建Python脚本命令
    py_script = "preprocess/run_pyscenic.py"
    cmd = [
        sys.executable, py_script,
        "--input_csv", input_csv,
        "--output_path", config.get('global', 'output_path'),
        "--species", config.get('global', 'species'),
        "--num_workers", config.get('pyscenic', 'num_workers'),
        "--seed", config.get('pyscenic', 'seed')
    ]
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr, bufsize=1, universal_newlines=True)
        print("pySCENIC脚本执行成功!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"pySCENIC脚本执行失败!")
        print(f"错误代码: {e.returncode}")
        print(f"请检查上面的错误信息，通常包括:")
        print(f"  - 输入文件格式是否正确")
        print(f"  - cisTarget文件是否存在")
        print(f"  - 内存是否足够")
        print(f"  - 依赖包是否正确安装")
        return False


def run_aucell(config):
    """
    运行AUCell脚本 (run_aucell.R)
    
    Args:
        config (configparser.ConfigParser): 配置对象
        
    Returns:
        bool: 成功返回True，失败返回False
    """
    print("=" * 60)
    print("步骤 3: 运行AUCell脚本 (run_aucell.R)")
    print("=" * 60)
    
    # 获取Rscript路径
    rscript_path = config.get('preprocess', 'rscript_path', fallback='Rscript')
    
    # 构建输入文件路径
    input_file = config.get('preprocess', 'input_file')
    regulon_file = os.path.join(config.get('global', 'output_path'), 'regulons.gmt')
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在: {input_file}")
        print("请确保Seurat对象文件存在")
        return False
    
    if not os.path.exists(regulon_file):
        print(f"错误: regulon文件不存在: {regulon_file}")
        print("请先运行pySCENIC步骤或确保regulons.gmt文件存在")
        return False
    
    # 构建R脚本命令
    r_script = "preprocess/run_aucell.R"
    cmd = [
        rscript_path, r_script,
        "--input_file", input_file,
        "--regulon_file", regulon_file,
        "--output_path", config.get('global', 'output_path'),
        "--method", config.get('aucell', 'method'),
        "--min_size", config.get('aucell', 'min_size'),
        "--batch_size", config.get('aucell', 'batch_size'),
        "--cores", config.get('aucell', 'cores')
    ]
    
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, stdout=sys.stdout, stderr=sys.stderr, bufsize=1, universal_newlines=True)
        print("AUCell脚本执行成功!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"AUCell脚本执行失败!")
        print(f"错误代码: {e.returncode}")
        return False


def check_dependencies():
    """
    检查依赖项是否可用
    
    Returns:
        bool: 所有依赖项可用返回True，否则返回False
    """
    print("检查依赖项...")
    
    # 检查R和Rscript - 这里只检查默认路径，具体路径在运行时检查
    try:
        subprocess.run(["Rscript", "--version"], check=True, capture_output=True)
        print("✓ 系统默认Rscript 可用")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ 系统默认Rscript 不可用")
        print("  请在配置文件中指定正确的Rscript路径")
        return False
    
    # 检查Python包
    try:
        import pandas
        print("✓ pandas 可用")
    except ImportError:
        print("✗ pandas 不可用，请安装: pip install pandas")
        return False
    
    # 检查pySCENIC相关命令
    try:
        subprocess.run(["arboreto_with_multiprocessing.py", "--help"], 
                      check=True, capture_output=True)
        print("✓ arboreto_with_multiprocessing.py 可用")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ arboreto_with_multiprocessing.py 不可用，请确保pySCENIC已安装")
        return False
    
    try:
        subprocess.run(["pyscenic", "--help"], check=True, capture_output=True)
        print("✓ pyscenic 可用")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("✗ pyscenic 不可用，请确保pySCENIC已安装")
        return False
    
    return True


def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(description="SCENIC Pipeline - GRN分析主控制脚本")
    parser.add_argument("--config", default="configs.txt", 
                       help="配置文件路径 (默认: configs.txt)")
    parser.add_argument("--skip-preprocess", action="store_true",
                       help="跳过预处理步骤")
    parser.add_argument("--skip-pyscenic", action="store_true", 
                       help="跳过pySCENIC步骤")
    parser.add_argument("--skip-aucell", action="store_true",
                       help="跳过AUCell步骤")
    parser.add_argument("--check-deps", action="store_true",
                       help="仅检查依赖项")
    parser.add_argument("--continue-on-error", action="store_true",
                       help="遇到错误时继续执行后续步骤")
    
    args = parser.parse_args()
    
    print("\n\nSCENIC Pipeline - GRN分析")
    print("=" * 60)
    
    # 检查依赖项
    if args.check_deps:
        if check_dependencies():
            print("\n所有依赖项检查通过!")
        else:
            print("\n依赖项检查失败，请安装缺失的依赖项")
            sys.exit(1)
        return
    
    if not check_dependencies():
        print("依赖项检查失败，请安装缺失的依赖项")
        sys.exit(1)
    
    # 读取配置文件
    try:
        print(f"读取配置文件: {args.config}")
        config = read_config(args.config)
        validate_config(config)
        print("配置文件验证通过!")
    except Exception as e:
        print(f"配置文件错误: {e}")
        sys.exit(1)
    
    # 检查imputed.mat.csv是否存在
    imputed_csv_path = os.path.join(config.get('global', 'output_path'), 'imputed.mat.csv')
    
    # 运行预处理
    if not args.skip_preprocess:
        if os.path.exists(imputed_csv_path):
            print(f"发现已存在的文件: {imputed_csv_path}")
            print("跳过预处理步骤，直接进行pySCENIC分析")
        else:
            if not run_preprocess(config):
                print("预处理步骤失败，终止执行")
                sys.exit(1)
    else:
        print("跳过预处理步骤")
    
    # 运行pySCENIC
    if not args.skip_pyscenic:
        if not run_pyscenic(config):
            print("pySCENIC步骤失败")
            if not args.continue_on_error:
                print("终止执行")
                sys.exit(1)
            else:
                print("继续执行后续步骤...")
        else:
            print("pySCENIC步骤成功完成")
    else:
        print("跳过pySCENIC步骤")
    
    # 运行AUCell
    if not args.skip_aucell:
        if not run_aucell(config):
            print("AUCell步骤失败")
            if not args.continue_on_error:
                print("终止执行")
                sys.exit(1)
            else:
                print("继续执行后续步骤...")
        else:
            print("AUCell步骤成功完成")
    else:
        print("跳过AUCell步骤")
    
    print("=" * 60)
    print("SCENIC Pipeline 执行完成!")
    print("=" * 60)


if __name__ == "__main__":
    main()
