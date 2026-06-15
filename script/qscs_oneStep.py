import sys
from pathlib import Path
import importlib.util
import argparse
import os 
import subprocess
import pandas as pd
import shutil
import glob

sys.path.insert(0, str(Path(__file__).parent))
try:
    from build_junction_lib_v5 import *
    from loaddata import *
    from qscs_calculationv2 import *
    from generatereport import copypict,update_csq,safe_copy_dir
except ImportError as e:
    print(f"Error: 必需的模块 overlap_resolver 未找到。请确保它在脚本所在目录的父目录中。")
    print(f"错误详情: {e}")
    sys.exit(1)

def busco_table(busco_txt:"full_table.tsv"):
    busco_info=pd.read_csv(busco_txt,sep="\t",comment="#",names=["busco_id","status","gene_id","score","length"],usecols=[0, 1, 2, 3, 4])
    genelist=busco_info[["gene_id","busco_id"]]
    genelist=genelist.fillna("0")
    return genelist

def build_junciton_lib(fasta,gtf,type,output,args):
    # 解析GTF文件
    print("Parsing GTF file...", file=sys.stderr)
    outputgtf=f"{gtf}.gtf"
    gff3_to_gtf(gtf,outputgtf)
    transcripts, exons, cdses = parse_gtf(outputgtf, args)
    #print(cdses)    
    if not transcripts:
        print("Error: No transcripts found in GTF/GFF file!", file=sys.stderr)
        sys.exit(1)
    
    if type=="exon":
        merged_exons = merge_exons_cds(exons, cdses)
    else:
        merged_exons = merge_exons_cds({}, cdses)
    
    # 加载基因组
    print("Loading genome sequence...", file=sys.stderr)
    genome = load_genome(fasta)
    print(f"Loaded {len(genome)} chromosomes", file=sys.stderr)

    # 处理每个transcript
    outputfile=f"{output}_lib.tsv"
    results=process_transcript(transcripts,merged_exons,genome)
    print(f"Writing results to {outputfile}...", file=sys.stderr)
    write_tsv(outputfile, results)
    print(f"Done! Processed {len(results)} transcripts.", file=sys.stderr)
    print(f"Sample output (first 3):", file=sys.stderr)
    for tr_id in sorted(results.keys())[:3]:
        print(results[tr_id], file=sys.stderr)
    return outputfile

def qscs_calculation(a_data,busco_data,genelist,output):
    #处理基因列表 最好为geneid\tbusco
    results = []
    missinglist=[]
    ouput_withoutIntron=open(f"gene_withoutIntron.txt",'w+')
    buscokeylist=list(busco_data.keys())
    print(len(buscokeylist))
    for i in range(len(genelist)):
        gene1,busco=genelist.iloc[i]
        if gene1!="0":
            if gene1 not in a_data:
                #ouput_withoutIntron.write(f'{gene1}\n')
                continue
            try:
                busco_data_tmp=busco_data[busco]
                cds_a=a_data[gene1]
                for gene_b in busco_data_tmp:           
                    score= calc_qscs(cds_a,busco_data_tmp[gene_b])
                    results.append({
                    'BUSCO': busco,
                    'GeneA': gene1,
                    'GeneB': gene_b,
                     **score
                    })
                # 将当前 BUSCO 标记为已处理
                if busco in buscokeylist:
                    buscokeylist.remove(busco)

            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"csq_calculation error in : {e.stderr.decode().strip()}")
        
    missinglist=buscokeylist
        

    ouputfile=f'{output}.csq.tsv'
    pd.DataFrame(results).to_csv(ouputfile, sep='\t', index=False)
    #print(f'busco:{len(busco_list)}')
    #print(f'adata:{len(busco_adata)}')
    print(f"Results have save in {ouputfile}.")  
    return ouputfile,missinglist

def qscs_calculation_genepairs(a_data,busco_data,genepair,output):
    results = []
    busco_list = list(busco_data.keys())
    busco_adata=[]

    ouput_withoutIntron=open(f"gene_withoutIntron.txt",'w+')
    for i in range(len(genepair)):
        gene1,busco_id,gene2 = genepair.iloc[i]
        if gene1 not in a_data:
            ouput_withoutIntron.write(f'{gene1}\n')
            continue
        try:
            busco_data_tmp=busco_data[busco_id][gene2]
            cds_a=a_data[gene1]
            score= calc_qscs(cds_a,busco_data_tmp)
            busco_adata.append(busco_id)
            results.append({
                    'BUSCO': busco_id,
                    'GeneA': gene1,
                    'GeneB': gene2,
                     **score
                    })
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"csq_calculation error in : {e.stderr.decode().strip()}")

    print(f'busco:{len(busco_list)}')
    print(f'adata:{len(busco_adata)}')
    missinglist= [x for x in busco_list if x not in busco_adata]
    ouputfile=f'{output}.csq.tsv'
    pd.DataFrame(results).to_csv(ouputfile, sep='\t', index=False)
    print(f"Results have save in {ouputfile}.")          
    return ouputfile,missinglist

def summarize(tsvfile,threshold_table,output,missinglist):
    module_name="qSCS_summarize.py"
    script_dir=os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir,module_name )      
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"scirpt {module_name} cannot found in {script_dir} ")
                
    # 加载模块
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
                
            # 确保有main函数
    if not hasattr(module, 'main'):
        raise AttributeError(f"Tool module '{module_name}' has no 'main' function")
                
    # 执行工具的主函数并传递参数
    print(f"\nExecuting tool: {module_name}")
    print("-"*60)
    import subprocess
    args_list=['-i',tsvfile,'-o',output,'-t',threshold_table,'--miss',missinglist]
    module.main(args_list)
    
def diamond(seq_q,seq_r,threads="8"):  
    diamond_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"diamond")
    try:
        print(f"running diamond")
        subprocess.run(
            [diamond_path, "makedb", "--in", os.path.abspath(seq_r), 
            "--db",os.path.abspath(seq_r)+".dmnd", "--threads", threads],
            check=True)
        subprocess.run(
            [diamond_path, "blastp", "--db", os.path.abspath(seq_r)+".dmnd", 
            "--query",  os.path.abspath(seq_q), "--threads", threads,"-o","seq2BUSCO.blast" ,"-k","2","--evalue","1e-5"],
            check=True)

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"DIAMOND failed with exit code {e.returncode}\nCommand: {e.cmd}\nError: {e.stderr}")
        sys.exit(0)
    
    blastdata=pd.read_csv("seq2BUSCO.blast",usecols=[0,1,11],names=['gene_id','busco_id','score'],sep="\t")
    #取第二个最优匹配的
    result = blastdata.loc[blastdata.groupby('gene_id')['score'].apply(lambda x: x.nlargest(1).idxmin())] 
    result = result.reset_index(drop=True)
    
    split_cols = result['busco_id'].str.split('-', n=1, expand=True)
    result['busco_id_'] = split_cols[0] 
    result['busco_gene'] = split_cols[1] 
    genepairs=result[["gene_id","busco_id_","busco_gene"]]
    genepairs.colname=["gene_id","busco_id","busco_gene"]
    print(genepairs)
    return genepairs

def _translate(sequence, table=1, stop_symbol="*"):
    """Helper function for frame translation"""
    from Bio.Seq import translate
    try:
        return translate(
            sequence,
            table=table,
            to_stop=False,
            stop_symbol=stop_symbol,
            cds=False
        )
    except Exception as e:
        raise ValueError(f"Translation failed: {str(e)}")

def getProt(libdata,output,isbusco=True):
    #如果是有busco列，可能还是要增加进去。
    fastafile=f'{output}.fa'

    with open(fastafile,'w+') as f:
        if isbusco:
            for busco in libdata:
                for geneid in libdata[busco]:
                    #print(libdata[busco])
                    f.write(f">{busco}-{geneid}\n")
                    seq="".join(libdata[busco][geneid]["seq"])
                    aaseq=_translate(seq)
                    f.write(f"{aaseq}\n")                   
        else:
            for geneid in libdata:
                f.write(f">{geneid}\n")
                seq="".join(libdata[geneid]["seq"])
                aaseq=_translate(seq)
                f.write(f"{aaseq}\n")
    return fastafile

def main(args_list=None):
    parser = argparse.ArgumentParser(
        description='Process GTF and genome files to extract splice site information.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 创建一个互斥组(--showdb 和其他参数互斥)
    group = parser.add_mutually_exclusive_group(required=True)
    
    # 选项1: 必须提供 -f/-g/-b 进行正常分析
    group.add_argument(
        '-f', '--fasta', 
        help='Input genome FASTA file (required for analysis mode)',
        metavar='FILE'
    )
    parser.add_argument(
        '-g', '--gtf', 
        help='Input GTF annotation file (required for analysis mode)',
        metavar='FILE'
    )
    parser.add_argument(
        '-b', '--library', 
        help='The highly conserved BUSCO gene library (eudicots/monocot/embryophyta)',
        metavar='LIB'
    )
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    
    # 选项2: 只显示数据库(--showdb)
    group.add_argument(
        '--showdb', 
        action='store_true',
        help='List all included built-in databases and exit'
    )
    
    # 其他可选参数
    parser.add_argument(
        '-l', '--busco-table', 
        default=None,
        help='Detailed BUSCO protein table (e.g., full_table.tsv)',
        metavar='FILE'
    )
    parser.add_argument(
        '--max-target-seqs', 
        type=int,  # 修正好你的参数名
        default=2, 
        help='Maximum number of Diamond hits'
    )
    parser.add_argument(
        '-o', '--output', 
        default='output', 
        help='Output prefix'
    )
    parser.add_argument(
        '-d', '--type', 
        default="CDS", 
        choices=['CDS', 'exon'],
        help='Analysis mode: CDS or exon'
    )
    parser.add_argument(
        '--report',
        action="store_true", 
        default=False,
        help='Whether to generate report(.md)'
    )
    args = parser.parse_args(args_list) if args_list else parser.parse_args()
    if args.showdb:
        print("The built-in database includes:")
        print("{:<12} {}".format("Plants:", "eudicots, monocots, embryophyta, eudicots_odb10, monocot_odb10, embryophyta_odb10"))
        print("{:<12} {}".format("Animals:", "aves, mammalia, sauropsida, actinopterygii, vertebrata"))
        sys.exit(0)
    else:
        if not all([args.fasta, args.gtf, args.library]):
            parser.error("Analysis mode requires -f/--fasta, -g/--gtf and -b/--library")
    

    if args.showdb:
        print(f'The built-in database includes:')
        print(f'plant database: eudicots, monocot, embryophyta,eudicots_odb10,monocot_odb10 and embryophyta_odb10')
        print(f'animal database: aves, mammalia, sauropsida,actinopterygii,vertebrata')
        sys.exit(0)
    
    #读取数据
    outputfile=build_junciton_lib(args.fasta,args.gtf,args.type,args.output,args)

    dirpath = Path(__file__).parent.parent
    datasetpath=os.path.join(dirpath,"Dataset")
    if args.library.lower() == "eudicots_odb10":
        busco_lib=os.path.join(datasetpath,"archive","eudicots_lib.txt")
        threshold_table=os.path.join(datasetpath,"eudicots_threshold.txt")
    elif  args.library.lower() == "eudicots":
        busco_lib=os.path.join(datasetpath,"eudicots_v2_lib.txt")
        threshold_table=os.path.join(datasetpath,"eudicots_v2_threshold.txt")
    elif  args.library.lower() == "monocot_odb10":
        busco_lib=os.path.join(datasetpath,"archive","liliopsida_lib.txt")
        threshold_table=os.path.join(datasetpath,"liliopsida_threshold.txt")
    elif  args.library.lower() == "monocot":
        busco_lib=os.path.join(datasetpath,"liliopsida_v3_lib.txt")
        threshold_table=os.path.join(datasetpath,"liliopsida_v3_threshold.txt")
    elif args.library.lower() == "embryophyta":
        busco_lib=os.path.join(datasetpath,"embryophyta_lib_v2.txt")
        threshold_table=os.path.join(datasetpath,"embryophyta_threshold_v2.txt")
    elif args.library.lower() == "embryophyta_odb10":
        busco_lib=os.path.join(datasetpath,"embryophyta_lib.txt")
        threshold_table=os.path.join(datasetpath,"archive","embryophyta_threshold.txt")
        #######动物####
    elif args.library.lower() == "aves":
        busco_lib=os.path.join(datasetpath,"Aves_libv2.txt")
        threshold_table=os.path.join(datasetpath,"Aves_thresholdv2.txt")
    elif args.library.lower() == "mammalia":
        busco_lib=os.path.join(datasetpath,"Mammalia_lib_v2.txt")
        threshold_table=os.path.join(datasetpath,"Mammalia_thresholdv2.txt")
    elif args.library.lower() == "sauropsida":
        busco_lib=os.path.join(datasetpath,"Sauropsida_lib.txt")
        threshold_table=os.path.join(datasetpath,"Sauropsida_threshold.txt")
    elif args.library.lower() == "vertebrata":
        busco_lib=os.path.join(datasetpath,"vertebrata_lib.txt")
        threshold_table=os.path.join(datasetpath,"vertebrata_threshold.txt")        
    elif args.library.lower() == "actinopterygii":
        busco_lib=os.path.join(datasetpath,"actinopterygii_v2_lib.txt")
        threshold_table=os.path.join(datasetpath,"actinopterygii_v2_threshold.txt")  
    else:
        busco_lib=f"{args.library}_lib.txt"
        threshold_table=f"{args.library}_threshold.txt"
        if os.path.exists(busco_lib) and os.path.exists(threshold_table):
            pass
        else:
            print("The built-in database includes:")
            print("{:<12} {}".format("Plants:", "eudicots, monocots, embryophyta, eudicots_odb10, monocot_odb10, embryophyta_odb10"))
            print("{:<12} {}".format("Animals:", "aves, mammalia, sauropsida, actinopterygii, vertebrata"))
            sys.exit(0)


    busco_data=load_busco_junction_data(busco_lib)
    a_data = load_junction_data(outputfile)
    print("The data has been loaded.")
    #形成对应列表
    
    if args.busco_table == None:
        seq1file=getProt(a_data,"seq1",False)
        buscofile=getProt(busco_data,"busco",True)
        genepair=diamond(seq1file,buscofile)
        print(genepair)
        qCSStsvfile, missinglist =qscs_calculation_genepairs(a_data,busco_data,genepair,args.output)
    
    else:
        genelist=busco_table(args.busco_table)
        print(genelist)
        qCSStsvfile, missinglist=qscs_calculation(a_data,busco_data,genelist,args.output)

    with open("Missing_busco.id",'w+')as output:
        for id in missinglist:
            output.write(f"{id}\n")
    summarize(qCSStsvfile,threshold_table,args.output,"Missing_busco.id")

    if args.report:
        templatePath=os.path.join(dirpath,"Template/Template.md")
        target_template=os.path.join(os.getcwd(),"Template.md")
        safe_copy_dir(templatePath, target_template)
        pattern = os.path.join(os.getcwd(), f"{args.output}*pie.png")
        copypict(pattern,"CSQ.png","Template.md")
        
        reportfile=os.path.join(os.getcwd(),"Template.md/CSQ-template.md")
        pattern = os.path.join(os.getcwd(), "*_stat.txt")
        CSQstat=glob.glob(pattern)[0]
        CSQdf = pd.read_csv(CSQstat,header=0,sep="\t")
        #print(CSQdf)
        busco_data = list(zip(CSQdf['label'], CSQdf['count'], CSQdf['frequency(%)']))
        update_csq( template_path=reportfile,new_table_data=busco_data)
        #除了CSQ-template.md其他报告模板都删了
        dir_path = Path("Template.md")
        exclude_file = "CSQ-template.md"
        for filename in os.listdir(dir_path):
            filepath = os.path.join(dir_path, filename)
            if os.path.isfile(filepath) and filename != exclude_file:
                os.remove(filepath)
        
        #除了RNAseqvenn.png其他图片都删了  
        dir_path = Path("Template.md/_v_images")
        exclude_file = "CSQ.png"
        for filename in os.listdir(dir_path):
            filepath = os.path.join(dir_path, filename)
            if os.path.isfile(filepath) and filename != exclude_file:
                os.remove(filepath)            
        print("The report has been generated.")
if __name__ == '__main__':
    main()

