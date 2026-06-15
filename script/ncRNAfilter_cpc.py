#!/usr/bin/env python3
import sys
import os
from pathlib import Path
import pandas as pd 
import importlib.util
import subprocess

def run_tool(args_list):
    # 获取 CPC2 脚本的绝对路径
    cpc2_dir = Path(__file__).parent / "CPC2_standalone-1.0.1" / "bin"
    script_path = cpc2_dir / "CPC2.py"  # 假设主脚本是 CPC2.py
    sys.path.insert(0, str(cpc2_dir))
    print(f"Script path: {script_path}")
    
    if not script_path.exists():
        raise FileNotFoundError(f"CPC2 script not found at {script_path}")
    
    print("Executing tool: CPC2")

    # 动态导入
    module_name = "CPC2"  # 不含 .py
    spec = importlib.util.spec_from_file_location(module_name, str(script_path))
    if spec is None:
        raise ImportError(f"Failed to create spec from {script_path}")
    
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    print("Module loaded successfully")
    module.main(args_list)  # 调用 CPC2 的 main 函数

def runncRNAfilter(transcript_fasta,input_gtf,output_prefix):
    cpc2_output=f"{output_prefix}_cpc"
    argslist=["-i",transcript_fasta,"-o",cpc2_output,"-r","--ORF"]
    run_tool(argslist)
    condingdata=pd.read_csv(f"{cpc2_output}.txt",sep="\t",header=0)
    noncodingID=condingdata.loc[condingdata["label"]=="noncoding"]["#ID"].tolist()
    
    output_gtf=f"{output_prefix}_ncRNAfilter.gtf"
    with open(input_gtf, "r+") as f_in, open(output_gtf, "w+") as f_out:
        for line in f_in:
            # Skip comment lines (optional)
            if line.startswith("#"):
                f_out.write(line)
                continue
            
            # Extract ID (assumes format like `gene_id "XYZ";`)
            if "gene_id" in line:  # Adjust field name if needed (e.g., "transcript_id")
                line_id = line.split('gene_id "')[1].split('"')[0].strip()
                
                # Only write if ID is NOT in the noncoding set
                if line_id not in noncodingID:
                    f_out.write(line)
        
    return output_gtf
