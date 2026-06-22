> **GSAqualifier: **A novel algorithm/toolkit assessing genome annotation accuracy through splice junction pattern analysis
>
> citation:
>

## <font style="color:rgb(64, 64, 64);background-color:rgb(252, 252, 252);">Installation</font>
```plain
wget https://github.com/Mavis-Mai/GSAqualifier/releases/download/1.0/GSAqualifier_v1.0.zip

unzip GSAqualifier_v1.0.zip && cd GSAqualifier_v1.0
mamba  create -n GSAqualifier -c conda-forge  python=3.10.5 
conda activate GSAqualifier

mamba  install  -c bioconda  stringtie  mash  regtools gffread
pip install -r requirements.txt
```

## <font style="color:rgb(31, 35, 40);">Contents</font>
<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/jpeg/29299085/1782093480095-27a12e0a-fb48-4278-893a-bff0eb717594.jpeg)

```plain
$GSAqualifier
```

```plain
Available tools:

Overview:
  csq                            The integration process of CSQ, including build, calculate and summarize.

Gobal assessment:
  csq2sp                         The integration process of CSQ comparing in two species
  rna                               The integration assessment process of RNA-seq, including the detection of gene number, splice junction, TES/TSS from RNA-seq, and reporter generation.

Other Tools:
  checkgene                      Check the gene (id/CDS-seq) by RNA-seq raw data

Usage: python GSAQualifier.py <tool_name>
Usage: GSAQualifier <tool_name>
```



### CSQ-mode
```plain
GSAqualifier csq --help
```

```plain
usage: GSAqualifier [-h] [-f FILE] [-g FILE] [-b LIB] [--debug] [--showdb] [-l FILE]
                    [--max-target-seqs MAX_TARGET_SEQS] [-o OUTPUT] [-d {CDS,exon}] [--report]

Process GTF and genome files to extract splice site information.

options:
  -h, --help            show this help message and exit
  -f FILE, --fasta FILE
                        Input genome FASTA file (required for analysis mode) (default: None)
  -g FILE, --gtf FILE   Input GTF annotation file (required for analysis mode) (default: None)
  -b LIB, --library LIB
                        The highly conserved BUSCO gene library (eudicots/monocot/embryophyta) (default: None)
  --debug               Enable debug mode (default: False)
  --showdb              List all included built-in databases and exit (default: False)
  -l FILE, --busco-table FILE
                        Detailed BUSCO protein table (e.g., full_table.tsv) (default: None)
  --max-target-seqs MAX_TARGET_SEQS
                        Maximum number of Diamond hits (default: 2)
  -o OUTPUT, --output OUTPUT
                        Output prefix (default: output)
  -d {CDS,exon}, --type {CDS,exon}
                        Analysis mode: CDS or exon (default: CDS)
  --report              Whether to generate report(.md) (default: False)
```

**<font style="background-color:#D8DAD9;">The options of library included:  eudicots, monocot, embryophta,  actinopterygii, aves, vertebrate, mammalia. </font>**



**Example: **

```plain
# Provide the genome file and annotation file, select the library 
GSAqualifier csq -f Arabidopsis_thaliana.genome.fna -g TAIR10_GFF3_genes.gff -b eudicots -o Ath

# If you have a corresponding library list, such as the gene mapping file for BUSCO (odb12), you can specify it using the -l option.
GSAqualifier csq -f Arabidopsis_thaliana.genome.fna -g TAIR10_GFF3_genes.gff -b eudicots -l  run_eudicotyledons_odb12/full_table.tsv -o Ath
```

**Results:**

+ <font style="background-color:#D8DAD9;">Ath_pie.svg</font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/svg/29299085/1782096124160-ee9de913-31ba-4738-a536-c390fbce142e.svg)

+ <font style="background-color:#D8DAD9;">Ath.csq.tsv : Detailed metrics table for all homologous gene comparisons</font>
+ <font style="background-color:#D8DAD9;">Ath_label_CSQ.txt:  Detailed metrics table for all homologous gene comparisons with Labels.</font>
+ **<font style="background-color:#D8DAD9;">Ath_max_CSQ.txt:</font>**<font style="background-color:#D8DAD9;"> Filtered metrics table retaining only reference-library genes with maximum CSQ values in homologous comparisons (Main Results)</font>

****

<font style="color:rgb(0, 0, 0);">The correspondence between the splice reference library and the BUSCO library is as follows： </font>

| CSQ library | BUSCO library |
| --- | --- |
| eudicots | eudicotyledons_odb12 |
| monocot | liliopsida_odb12 |
| embryophyta | embryophyta_odb12 |
| vertebrata | vertebrata_odb12 |
| mammalia | mammalia_odb12 |
| aves | aves_odb12 |
| actinopterygii | actinopterygii_odb12 |
| sauropsida | sauropsida_odb12 |


### CSQ2sp-mode
```plain
GSAqualifier csq2sp
```

```plain
usage: GSAqualifier [-h] -f QUERY_FASTA -g QUERY_GTF -F REF_FASTA -G REF_GTF [-l GENEPAIRS] [-k MAX_TARGET_SEQS]
                    [-d TYPE] [-t THREADS] [--report] [-o OUTPUT] [--fig_format FIG_FORMAT] [--debug]

Process GTF and genome files to extract splice site information.
options:
  -h, --help            show this help message and exit
  -f QUERY_FASTA, --query_fasta QUERY_FASTA
                        Input query genome FASTA file (required) (default: None)
  -g QUERY_GTF, --query_gtf QUERY_GTF
                        Input query GTF annotation file (required) (default: None)
  -F REF_FASTA, --ref_fasta REF_FASTA
                        Input reference genome FASTA file (required) (default: None)
  -G REF_GTF, --ref_gtf REF_GTF
                        Input reference GTF annotation file (required) (default: None)
  -l GENEPAIRS, --genepairs GENEPAIRS
                        Provide high-confidence gene pair information; If this file is not provided, reciprocal
                        alignment of diamond will be used to map gene pairs (default: None)
  -k MAX_TARGET_SEQS, ---max-target-seqs MAX_TARGET_SEQS
                        hit the number of Diamonds (default: 2)
  -d TYPE, --type TYPE  mode:CDS or exon. (default: CDS)
  -t THREADS, --threads THREADS
                        cpu number (default: 8)
  --report              Whether to generate report(.md) (default: False)
  -o OUTPUT, --output OUTPUT
                        output prefix (default: output)
  --fig_format FIG_FORMAT
                        The selection of graphic format. such as:svg, png, pdf (default: svg)
  --debug               Enable debug mode (default: False)

```



**Example: **

```plain
# Provide the target genome file and annotation file, and the high quality reference genome file and annotation file
GSAqualifier csq -F Arabidopsis_thaliana.genome.fna -G TAIR10_GFF3_genes.gff -f Oryza_sativa_NIP.genomic.fa -g Oryza_sativa_NIP.genomic.gff -o Ath-Osa
```

**Results:**

+ <font style="background-color:#D8DAD9;">CSQ-density.png&CSQ-density.svg: The disribution of CSQ in all homologous genes. Q1 can represent the overall similarity of splicing patterns between homologous genes across two species.</font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/png/29299085/1782097749883-a5a13a5f-0cb2-48b4-b1b5-326f0cd9fa92.png)

+ **<font style="background-color:#D8DAD9;">Ath_Osa.csq.tsv: </font>**<font style="background-color:#D8DAD9;">Detailed metrics table for all homologous gene in two species. (Main Results)</font>
+ <font style="background-color:#D8DAD9;">reciprocal_pairs.tsv: the reciprocal pairs in two species.</font>



### RNA-mode
```plain
GSAqualifier rna
```

```plain
usage: GSAqualifier [-h] -g GXF -b BAM_FLODER_FILE [BAM_FLODER_FILE ...] -f FASTA [-o OUTPUT] [--force]
                    [--ncRNA_filter] [-m MIN_OVERLAP] [-c MIN_SUPPORT] [-d MAX_GOOD] [-D MAX_MODERATE] [-t THREADS]
                    [-r RECORDTYEP] [--idItem IDITEM] [-w WINDOW_SIZE] [--item ITEM] [--report]

Compare the gene regions, junction, TES/TSS between RNA-seq files and annotation

optional arguments:
  -h, --help            show this help message and exit
  -g GXF, --gxf GXF     gene anotation (gff/gtf). (default: None)
  -b BAM_FLODER_FILE [BAM_FLODER_FILE ...], --bam_floder_file BAM_FLODER_FILE [BAM_FLODER_FILE ...]
                        A directory with RNA-seq files(bam) or a bam file. (default: None)
  -f FASTA, --fasta FASTA
                        gene assembly (fasta). (default: None)
  -o OUTPUT, --output OUTPUT
                        output prefiex (default: output)
  --force               Force regenerate the result (overwrite the existing file) (default: False)
  -m MIN_OVERLAP, --min_overlap MIN_OVERLAP
                        Minimum overlap ratio (0-1) (default: 0.5)
  -c MIN_SUPPORT, --min_support MIN_SUPPORT
                        The minimum number threshold of supported read segments (default: 5)
  -d MAX_GOOD, --max_good MAX_GOOD
                        maximum distance threshold for "1-Good" (default: 200 bp)
  -D MAX_MODERATE, --max_moderate MAX_MODERATE
                        maximum distance threshold for "2-Moderate" (default: 500 bp)
  -t THREADS, --threads THREADS
                        Number of parallel threads (default:8 CPUs) (default: 8)
  -r RECORDTYEP, --recordTyep RECORDTYEP
                        extract junction from feature type (default: exon)
  --idItem IDITEM       extract gene name from item. with --gxfform (default: )
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        size(bp) of sliding window (default: 100)
  --item ITEM           extract exon information from item. with --gxfform (default: )
  --report              Whether to generate report(.md) (default: False)

```



**Example: **

```plain
# Provide the target genome file and annotation file, and the alignment files/folds(*.bam) of RNA-seq data 
GSAqualifier rna -f Arabidopsis_thaliana.genome.fna -g TAIR10_GFF3_genes.gff -b bamfile/  -o Ath -t 10
```

**Results:**

**<font style="background-color:#D8DAD9;">Running_dir</font>**<font style="background-color:#D8DAD9;">:</font>

<font style="background-color:#D8DAD9;">-Ath_riskscore_distribution.svg:</font>

<font style="background-color:#D8DAD9;">The disribution of risk score (RS) in all expression genes. Top 95% value of RS can represent the overall splicing quality of genes in target specie.</font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/png/29299085/1782110807607-bf7557ff-fbff-44c8-b48e-629156785e00.png =600x)

**<font style="background-color:#D8DAD9;">Splice_junctions</font>**<font style="background-color:#D8DAD9;">:</font>

<font style="background-color:#D8DAD9;">-SJ_output_pieplot.svg</font>

    - **<font style="color:rgb(0, 0, 0);">LEVEL A</font>**<font style="color:rgb(0, 0, 0);">: SJs fully supported by transcriptome data.</font>
    - **<font style="color:rgb(0, 0, 0);">LEVEL B</font>**<font style="color:rgb(0, 0, 0);">: SJs partially conflicting with transcriptome data.</font>
    - **<font style="color:rgb(0, 0, 0);">NA</font>**<font style="color:rgb(0, 0, 0);">: No transcriptomic support, possibly not expressed.</font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/png/29299085/1782110906597-8b4dc9da-0670-43c6-96c3-9c40dc172a2a.png)

**<font style="background-color:#D8DAD9;">Gene_Region</font>**

<font style="background-color:#D8DAD9;">-Gene_output_venn.png </font>

<font style="background-color:#D8DAD9;">Left: 87.9% of </font>_<font style="background-color:#D8DAD9;">de novo</font>_<font style="background-color:#D8DAD9;"> RNA-seq loci match known gene annotations, with some novel loci remaining.</font>

<font style="background-color:#D8DAD9;">Right: 85.0% of gene annotations align with </font>_<font style="background-color:#D8DAD9;">de novo</font>_<font style="background-color:#D8DAD9;"> assembled transcripts.</font>

<font style="background-color:#D8DAD9;"></font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/png/29299085/1782111232999-88b12fea-1e55-43e7-bd1b-425bb52ffdd1.png)

### BuildDB
<font style="color:rgb(0, 0, 0);">For species not covered by standard references, </font>**<font style="color:rgb(0, 0, 0);">GSAqualifier</font>**<font style="color:rgb(0, 0, 0);"> supports user-defined library construction.</font>

<!-- 这是一张图片，ocr 内容为： -->
![](https://cdn.nlark.com/yuque/0/2026/jpeg/29299085/1782111957883-cb8e984b-1ae3-4f4d-b36c-40478f6733c2.jpeg)

```plain
python GSAqualifier/script/builddatabase_one_step.py
```

```plain
usage: builddatabase_one_step.py [-h] -f GENOME_DIR -l BUSCO_TSV_DIR [-o OUTPUT] [-n CPU] [--debug]

Self-build conservative database.

options:
  -h, --help            show this help message and exit
  -f GENOME_DIR, --genome_dir GENOME_DIR
                        Input directory including genome(fa) and annotation(gxf) of different species (default: None)
  -l BUSCO_TSV_DIR, --busco_tsv_dir BUSCO_TSV_DIR
                        Input directory including busco gene pairs list(tsv) of different species (default: None)
  -o OUTPUT, --output OUTPUT
                        output prefix of database directory (default: output)
  -n CPU, --cpu CPU     cpu number used in parallel jobs (default: 10)
  --debug               Enable debug mode (default: False)

```

**Example: **

**-genome: **Directory for storing the genomic files and annotation files of species that are constructed <font style="color:rgb(0, 0, 0);">user-defined library.</font>

```plain
|____genome
| |____sp3
| | |____GCA_978657495.1_TAIR12_genomic_clean.fna
| | |____GCA_978657495.1_TAIR12_genomic.gff
| |____sp2
| | |____GWHERGL00000000.genome.fasta
| | |____GWHERGL00000000.gff
| |____sp1
| | |____Mperfoliatum_583_v2.0.fa
| | |____Mperfoliatum_583_v2.1.gene.gff3

```

**-busco: The mapping files stored all gene-to-ortholog relationships (e.g., OG1\tComplete\tAJChr06G01814). These files are preserved in a standardized format, such as BUSCO output files, ensuring consistent orthologous mappings across all species studied.**

```plain
|____busco
| |____Brapa_full_table.tsv
| |____Mper_full_table.tsv
| |____conserved_buscos.id
| |____conserved_buscos.txt
| |____TAIR12_full_table.tsv
```

```plain
python GSAqualifier/script/builddatabase_one_step.py -f genome/ -l busco/ -n 10 -o test
```

**Results:**

<font style="background-color:#D8DAD9;">-test_threshold.txt</font>

<font style="background-color:#D8DAD9;">-test_lib.txt</font>

<font style="background-color:#D8DAD9;">The threshold file and the library file generated and used by CSQ mode.</font>

```plain
GSAqualifier csq -f Mperfoliatum_583_v2.0.fa -g  Mperfoliatum_583_v2.1.gene.gff3 -b <path>/test
```

