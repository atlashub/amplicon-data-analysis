# 扩增子数据分析之： dix-seq 生物信息数据分析工作流解析


`约定`: 本材料适合安装 `biostack` 虚拟机或者 `biostack` 一体机的机器, 材料中涉及到的工具和数据库都默认已经安装。


### 0. 准备工作：


使用 `MobaXterm` 登录 `biostack` 生物信息一体机终端，

添加辅助程序路径至环境变量：

```bash
export PATH=$PATH:/biostack/tools/protocols/dix-seq-0.0.2/utils
export PATH=$PATH:/biostack/tools/protocols/dix-seq-0.0.2/binaries
biostack
```

准备好原始数据，样品信息表、以及引物信息。
软件需要安装：`R, seqtk, usearch, tsv-utils, fastx-utils, atlas-utils, biom以及fasttree、 Trimal, mafft, KronaTools 以及 phylommand`

`primers`文件

    $ cat primers.txt  | head -n 4


    >304F
    CCTAYGGGRBGCASCAG
    >515R
    GGACTACNNGGGTATCTAAT


未拆分数据文件：

    $ ls  raw_data


    illumina-1.fastq  illumina-2.fastq


`barcodes.txt` 文件：

    $ cat barcodes.txt | head -n 6


    A-1     CGACTAC GTATA
    A-2     GATGA   CGATATC
    A-3     AACAT   TACTATG
    B-1     CTTGTA  GCCAAT
    B-2     GCGAAGT CGAGG
    B-3     TCATT   ATGAGTC


`mapping_file.txt` 文件 

    $ cat mapping_file.txt | head


    #SampleID       Description
    A-1     A
    A-2     A
    A-3     A
    B-1     B
    B-2     B
    B-3     B


数据演示当前目录设定为： `/project/dix-demo`；



### 1. 数据拆分

根据`PCR`扩增实验使用的`PCR primer`序列和`barcode`序列，使用 `atlas-utils` [demultiplex](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/demultiplex/demultiplex.md) 进行数据拆分。

```bash
atlas-utils  demultiplex -b barcodes.txt  -1  raw_data/illumina-1.fastq  -2 raw_data/illumina-2.fastq  -l 5  -d raw_data
ls raw_data
```

    A-1_1.fastq  A-1_2.fastq  A-2_1.fastq  A-2_2.fastq  A-3_1.fastq  A-3_2.fastq  B-1_1.fastq  B-1_2.fastq  B-2_1.fastq  B-2_2.fastq  B-3_1.fastq  B-3_2.fastq  illumina-1.fastq  illumina-2.fastq

`注意：`分拆的数据如果想保留barcode序列可以使用 atlas-utils [split](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/split/split.md) 程序：


```bash
atlas-utils split -b barcodes.txt  -1  raw_data/illumina-1.fastq  -2 raw_data/illumina-2.fastq  -l 5  -d raw_data
```

得到的分拆后序列可以下一轮数据预处理。



### 2.质量控制，去除接头和低质量序列

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/trimming/reads/
mkdir -p /project/dix-demo/data_analysis/trimming/report/
```

使用 [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) 去除拆分序列中可能可能存在的接头和切掉低质量的碱基， 针对一些小片段扩增产物，测序数据可能包含接头序列。

```bash
java -jar /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/Trimmomatic/trimmomatic.jar PE \
    -threads 12 \
    -phred33    \
    /project/dix-demo/raw_data/A-1_1.fastq /project/dix-demo/raw_data/A-1_2.fastq \
    /project/dix-demo/data_analysis/trimming/reads/A-1_1.fastq \
    /project/dix-demo/data_analysis/trimming/reads/A-1_1.singleton.fastq \
    /project/dix-demo/data_analysis/trimming/reads/A-1_2.fastq \
    /project/dix-demo/data_analysis/trimming/reads/A-1_2.singleton.fastq \
    ILLUMINACLIP:/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:15 \
    LEADING:3 \
    TRAILING:3 \
    MINLEN:36
```

针对多样本可使用 `bash` 进行循环迭代：

```bash
cut -f1  /project/dix-demo/mapping_file.txt  | \
    grep -v "#" | \
    perl -ane 'print qq{java -jar /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/Trimmomatic/trimmomatic.jar PE -threads 12 -phred33 /project/dix-demo/raw_data/$F[0]\_1.fastq /project/dix-demo/raw_data/$F[0]\_2.fastq /project/dix-demo/data_analysis/trimming/reads/$F[0]_1.fastq /project/dix-demo/data_analysis/trimming/reads/$F[0]_1.singleton.fastq /project/dix-demo/data_analysis/trimming/reads/$F[0]_2.fastq /project/dix-demo/data_analysis/trimming/reads/$F[0]_2.singleton.fastq ILLUMINACLIP:/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 ;\n}' | \
    bash
```

**参数解析 :**


    命令模式:  SE指定单端数据 ，PE指定双端数据，对于PE模式: 2个输入文件（正向和反向reads）和4个输出文件（正向配对、正向未配对、反向配对和反向未配对）
    threads： 设定线程数
    phred<qual>： 设定碱基质量编码模式，两种模式 `-phred33` 或 `-phred64`， 如果未指定质量编码，则将自动识别, `fastq` 质量值编码可以参考：[https://wiki.bits.vib.be/index.php/FASTQ](https://wiki.bits.vib.be/index.php/FASTQ)
    
    ILUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>： 该步骤用于寻找并去除Illumina接头    
    
        <fastaWithAdaptersEtc1>:     指定包含所有接头、PCR序列等的fasta文件的路径;
        <seed mismatches>：          指定允许执行完全匹配的最大不匹配数；
        <palindrome clip threshold>：指定两个成对接头reads之间的匹配对于双端回文read比对的精度；
        <Simple clip threshold>：    指定任何接头等序列与read之间的匹配精度阈值；
        <minAdapterLength>:          设定接头的最小长度。如果未指定，默认为8个bp;
        <keepBothReads>：            在回文模式检测到read测穿并删除接头序列后，反向读取包含与正向读取相同的序列信息。因此，默认行为是完全删除反向读取。通过为该参数指定true，反向读取也将被保留
    
    SLIDINGWINDOW:<windowSize>:<requiredQuality>: 滑窗修剪, 它从5'端开始扫描，当窗口内的平均质量低于阈值时，剔除该窗口内的所有碱基;
    
        <windowSize>:       设定窗口大小, 覆盖碱基数量；
        <requiredQuality>:  设定平均质量。
    
    HEADCROP:<length> ：切除read起始端低于阈值的碱基
    
        <length> 从`read` 起始端开始要切除的长度
    
    TRAILING:<quality>：从末端移除低质量的碱基。只要碱基的质量值低于阈值，则切除该碱基，并调查下一个碱基（因为Trimmomatic从3'primeend开始，将是位于刚切除碱基之前的碱基）。 此方法可用于去除Illumina低质量段区域（质量分数标记为2）, 
     
        <quality>: 指定保留碱基所需的最低质量
     
    MINLEN:<length>: 设置保留reads的最小长度。
        <length>:可被保留的最短 read 长度。


`trimmomatic.jar`命令参考：[http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


**统计序列数目 :**


切除`adaptor`和低质量的碱基后的序列，然后使用`atlas-utils` [fqchk](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/fqchk/fqchk.md) 进行序列质量统计

```bash
cat /project/dix-demo/data_analysis/trimming/reads/A-1_1.fastq   \
    /project/dix-demo/data_analysis/trimming/reads/A-1_2.fastq | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils fqchk -p -q 33 -l A-1 -  \
    >/project/dix-demo/data_analysis/trimming/reads/A-1.txt
```

可使用 `bash` 进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
grep -v "#" | \
perl -ane 'print qq{ cat /project/dix-demo/data_analysis/trimming/reads/$F[0]_1.fastq /project/dix-demo/data_analysis/trimming/reads/$F[0]_2.fastq | /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils fqchk -p -q 33 -l $F[0] - > /project/dix-demo/data_analysis/trimming/reads/$F[0].txt ;\n}' | \
bash
```

**统计信息 :**

将所有样本的碱基序列质量统合为`trimming.stats.txt`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/usearch-utils trimming           \
     /project/dix-demo/mapping_file.txt /project/dix-demo/data_analysis/trimming/reads \
     >/project/dix-demo/data_analysis/trimming/report/trimming.stats.txt
```



### 3. 合并双端序列

**创建目录 :**

```bash
mkdir -p    /project/dix-demo/data_analysis/mergepairs/reads
mkdir -p    /project/dix-demo/data_analysis/mergepairs/report
```

**合并双端序列 :**


使用`usearch` [-fastq_mergepairs](https://www.drive5.com/usearch/manual/merge_options.html) 命令合并`A1-1`的双端序列，并输出结果到`data_analysis/mergepairs/reads/`目录下:

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -fastq_mergepairs /project/dix-demo/data_analysis/trimming/reads/A-1_1.fastq \
    -reverse /project/dix-demo/data_analysis/trimming/reads/A-1_2.fastq \
    -fastqout /project/dix-demo/data_analysis/mergepairs/reads/illumina.fastq \
    -fastqout_notmerged_fwd  /project/dix-demo/data_analysis/mergepairs/reads/A-1.notmerged_fwd.fastq \
    -fastqout_notmerged_rev  /project/dix-demo/data_analysis/mergepairs/reads/A-1.notmerged_rev.fastq \
    -log  /project/dix-demo/data_analysis/mergepairs/reads/A-1.log \
    -report /project/dix-demo/data_analysis/mergepairs/reads/A-1.report \
    -threads 1 \
    -fastq_minmergelen 250 \
    -fastq_maxmergelen 500 \
    -fastq_maxdiffs 10  \
    -fastq_pctid 90  \
    -fastq_minovlen 16   \
    -fastq_trunctail 2 \
    -fastq_minlen 64
```

可使用 `bash` 进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
grep -v "#" |                                 \
perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch -fastq_mergepairs /project/dix-demo/data_analysis/trimming/reads/$F[0]_1.fastq -reverse /project/dix-demo/data_analysis/trimming/reads/$F[0]_2.fastq -fastqout /project/dix-demo/data_analysis/mergepairs/reads/$F[0].fastq -fastqout_notmerged_fwd  /project/dix-demo/data_analysis/mergepairs/reads/$F[0].notmerged_fwd.fastq -fastqout_notmerged_rev  /project/dix-demo/data_analysis/mergepairs/reads/$F[0].notmerged_rev.fastq -log  /project/dix-demo/data_analysis/mergepairs/reads/$F[0].log -report /project/dix-demo/data_analysis/mergepairs/reads/$F[0].report -threads 1 -fastq_minmergelen 250 -fastq_maxmergelen 500 -fastq_maxdiffs 10  -fastq_pctid 90  -fastq_minovlen 16   -fastq_trunctail 2 -fastq_minlen 64 ;\n}' | \
bash
```

`usearch -fastq_mergepairs`参数解析

    -fastq_mergepairs       正向FASTQ文件名（输入）
    -reverse                反向FASTQ文件名（输入）
    -fastqout               合并后的FASTQ文件名（输出）
    -fastqout_notmerged_fwd 未合并的正向序列FASTQ文件名（输出）
    -fastqout_notmerged_rev 未合并的反向序列FASTQ文件名（输出）
    -log                    log文件名（输出）
    -report                 总结报告文件名（输出）
    -threads                线程数
    -fastq_minmergelen      合并序列的最小长度
    -fastq_maxmergelen      合并序列的最大长度  
    -fastq_maxdiffs         最大错配数，默认为5
    -fastq_pctid            最小序列质量，默认为90%
    -fastq_minovlen         合并后的序列如果低于指定值则过滤，默认为16
    -fastq_trunctail        在第一个小于指定Q值处截断序列，默认为
    -fastq_minlen           如果-fastq_trunctail截断后的序列小于指定值，则过滤掉序列

`usearch -fastq_mergepairs` 命令参考：[http://www.drive5.com/usearch/manual/merge_options.html](http://www.drive5.com/usearch/manual/merge_options.html)

** 统计序列信息 :**

使用`atlas-utils`  [fqchk](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/fqchk/fqchk.md)  统计序列质量

```bash
cat /project/dix-demo/data_analysis/mergepairs/reads/A-1.fastq |                            /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils fqchk -q 33 -l A-1 - \
    > /project/dix-demo/data_analysis/mergepairs/reads/A-1.txt ;

```

可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | grep -v "#" | perl -ane 'print qq{cat /project/dix-demo/data_analysis/mergepairs/reads/$F[0].fastq | /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils fqchk -q 33 -l $F[0] - > /project/dix-demo/data_analysis/mergepairs/reads/$F[0].txt ;\n}' |bash

```

**统计信息 :**

合并所有`mergepairs/reads/`目录下后缀为`.txt`的文件并传递给`tsv-utils` [view](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/view/view.md) 去除中间的注释行, 使用`usearch-utils mergepairs`统计`log`文件

```bash
cat /project/dix-demo/data_analysis/mergepairs/reads/*.txt |                  \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils view - \
    >/project/dix-demo/data_analysis/mergepairs/report/sequencing.stats.txt
```

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/usearch-utils \mergepairs             /project/dix-demo/mapping_file.txt \
    /project/dix-demo/data_analysis/mergepairs/reads  \
    >/project/dix-demo/data_analysis/mergepairs/report/mergepairs.stats.txt
```

### 4. 引物识别和切除引物


**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/primer_strip/reads
mkdir -p /project/dix-demo/data_analysis/primer_strip/report
```

**切除引物 :**


使用`usearch` [search_pcr2](http://www.drive5.com/usearch/manual/cmd_search_pcr2.html) 切除引物序列并输出结果到`primer_strip/reads/`目录

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -search_pcr2 /project/dix-demo/data_analysis/mergepairs/reads/illumina.fastq \
    -fwdprimer CCTAYGGGRBGCASCAG     \
    -revprimer GGACTACNNGGGTATCTAAT  \
    -minamp 250                      \
    -maxamp 480                      \
    -strand both                     \
    -maxdiffs 2                      \
    -tabbedout /project/dix-demo/data_analysis/primer_strip/reads/illumina.txt   \
    -fastqout  /project/dix-demo/data_analysis/primer_strip/reads/illumina.fastq \
    -log /project/dix-demo/data_analysis/primer_strip/reads/illumina.log
```

`usearch -search_pcr2`参数解析:

    -fwdprimer      指定正向引物序列
    -revprimer      指定反向引物序列
    -minamp         指定引物之间的序列最小长度
    -maxamp         指定引物之间的序列最大长度
    -strand         指定plus或者both
    -maxdiffs       指定引物的最大允许错配属
    -tabbedout      输出结果列表（输出）
    -fastqout       输出引物之间片段的FASTQ文件（输出）
    -log            日志（输出）


可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
   grep -v "#" |                              \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch -search_pcr2 /project/dix-demo/data_analysis/mergepairs/reads/$F[0].fastq -fwdprimer CCTAYGGGRBGCASCAG -revprimer GGACTACNNGGGTATCTAAT -minamp 250 -maxamp 480 -strand both -maxdiffs 2 -tabbedout /project/dix-demo/data_analysis/primer_strip/reads/$F[0].txt -fastqout  /project/dix-demo/data_analysis/primer_strip/reads/$F[0].fastq -log /project/dix-demo/data_analysis/primer_strip/reads/$F[0].log ;\n}' | \
   bash
```


`usearch -search_pcr2`命令参考：[http://www.drive5.com/usearch/manual/cmd_search_pcr2.html](http://www.drive5.com/usearch/manual/cmd_search_pcr2.html)


**统计信息 :**


使用`usearch-utils pcrsearch`统计`primer_strip/reads`中的文件

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/usearch-utils pcrsearch               /project/dix-demo/mapping_file.txt \
    /project/dix-demo/data_analysis/primer_strip/reads  \
   >/project/dix-demo/data_analysis/primer_strip/report/pcrsearch.stats.txt
```


### 5. 添加序列注释信息和合并样本

**创建目录:**

```bash
mkdir -p    /project/dix-demo/data_analysis/zotu/labels
mkdir -p    /project/dix-demo/data_analysis/zotu/reads
mkdir -p    /project/dix-demo/data_analysis/zotu/filter
```

使用`fastx-utils` [rename]() 修改`primer_strip/reads/A-1.fastq`序列的名字前缀为`A-1`，然后使用`fastx-utils` [label]() 在序列标签后缀添加`;sample=A-1;`并导出`zotu/labels/A-1.fastq`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils rename \
    /project/dix-demo/data_analysis/primer_strip/reads/A-1.fastq A-1 | \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils label -  ";sample=A-1;" \
>/project/dix-demo/data_analysis/zotu/labels/A-1.fastq
```

此时序列标识符为: 

    @illumina_1;sample=A-1;


可使用 `bash` 进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  |  \
   grep -v "#" |                               \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils rename /project/dix-demo/data_analysis/primer_strip/reads/$F[0].fastq $F[0] | /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils label -  ";sample=$F[0];" > /project/dix-demo/data_analysis/zotu/labels/$F[0].fastq ;\n}' | \
   bash
```


### 6. 质量控制, 去除低质量序列


使用`usearch` [fastq_filter](http://drive5.com/usearch/manual/cmd_fastq_filter.html) 过滤`zotu/labels/A-1.fastq`中的序列并导出 `zotu/filter/A-1.fasta`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -fastq_filter /project/dix-demo/data_analysis/zotu/labels/A-1.fastq \
    -relabel A-1_Filt \
    -fastq_maxee 1 \
    -fastaout  /project/dix-demo/data_analysis/zotu/filter/A-1.fasta \
    -threads  72 \
    -log /project/dix-demo/data_analysis/zotu/filter/A-1.txt
```

`usearch -fastq_filter`参数解析：

    -relabel        生成新的前缀
    -fastq_maxee    过滤掉大于指定E值的序列
    -fastaout       输出的FASTA文件（输出）
    -threads        线程数
    -log            日志文件（输出）


`usearch -fastq_filter` 命令参考：[http://www.drive5.com/usearch/manual/cmd_fastq_filter.html](http://www.drive5.com/usearch/manual/cmd_fastq_filter.html)


可使用 `bash` 进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
   grep -v "#" |                              \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch -fastq_filter /project/dix-demo/data_analysis/zotu/labels/$F[0].fastq -relabel $F[0]_Filt -fastq_maxee 1 -fastaout  /project/dix-demo/data_analysis/zotu/filter/$F[0].fasta -threads  72 -log /project/dix-demo/data_analysis/zotu/filter/$F[0].txt ;\n}' | \
bash
```

**注意事项:** 质量控制后的序列仅仅用于构建`ZOTU`表； 丰度定量使用合并后的`strip.fastq`文件。



### 7. 去除非冗余序列，添加大小注释


合并`zotu/filter/`中所有的`fasta`序列并导出`zotu/reads/filtered.fasta`

```bash
cat /project/dix-demo/data_analysis/zotu/filter/*.fasta \
    >/project/dix-demo/data_analysis/zotu/reads/filtered.fasta
```

使用`atlas-utils` [uniques](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/uniques/uniques.md) 获得所有非冗余序列，添加大小注释，导出`zotu/reads/derep.fasta`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils uniques -a  \
    /project/dix-demo/data_analysis/zotu/reads/filtered.fasta   \
    >/project/dix-demo/data_analysis/zotu/reads/derep.fasta
```


### 8. unoise3 去噪

**创建目录 :**

```bash
mkdir -p  /project/dix-demo/data_analysis/zotu/unoise/
```

**去噪 :**

使用`usearch` [unoise3](http://www.drive5.com/usearch/manual/cmd_unoise3.html) 进行去噪，形成代表序列。

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -unoise3 /project/dix-demo/data_analysis/zotu/reads/derep.fasta \
    -minsize 8  \
    -unoise_alpha 2 \
    -zotus /project/dix-demo/data_analysis/zotu/unoise/denoise.fasta \
    -tabbedout /project/dix-demo/data_analysis/zotu/unoise/unoise.txt \
    -log /project/dix-demo/data_analysis/zotu/unoise/unoise.log
```

`usearch -unoise3`参数解析：

    -minsize        指定最小丰度，如果序列低于指定丰度则过滤，默认为8
    -unoise_alpha   指定alpha参数，默认为2
    -zotus          降噪后的FASTA文件（输出）
    -tabbedout      序列的处理结果列表（输出）
    -log            日志文件（输出）


`usearch -unoise3` 命令参考：[http://www.drive5.com/usearch/manual/cmd_unoise3.html](http://www.drive5.com/usearch/manual/cmd_unoise3.html)



### 9. 删除非特异性扩增，格式化ZOTU代表序列


**删除非特异性扩增 :**

使用`usearch` [usearch_global](http://www.drive5.com/usearch/manual/cmd_usearch_global.html)进行与数据库进行全局比对, 这一步设置的相似度比较低，可以过滤掉非特异性扩增序列，这一步可以不做。

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch                \
    -usearch_global /project/dix-demo/data_analysis/zotu/unoise/denoise.fasta   \
    -db /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/rdp_16s_v16_sp.fasta \
    -id 0.4                                                                     \
    -maxhits 1                                                                  \
    -blast6out /project/dix-demo/data_analysis/zotu/unoise/align.txt            \
    -strand both                                                                \
    -log /project/dix-demo/data_analysis/zotu/unoise/align.log
```

`usearch -usearch_global`参数解析：

    -db         指定数据库文件
    -id         最小的序列一致性
    -maxhits    最大比对数
    -blast6out  blast格式的输出列表
    -strand     指定plus或者both
    -log        日志文件

`usearch -usearch_global`命令参考： [http://www.drive5.com/usearch/manual/cmd_usearch_global.html](http://www.drive5.com/usearch/manual/cmd_usearch_global.html)


**格式化序列 :**

```bash
mkdir -p    /project/dix-demo/data_analysis/zotu/report/
```

**格式化结果文件 :**

- 使用`tsv-utils` [cut](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/cut/cut.md)获取全局比对后的序列ID
- 使用`fastx-utils` [subseq](https://github.com/atlashub/biostack-suits-docs/blob/master/fastx-utils/subseq/subseq.md)获取全局比对后的序列
- 使用`fastx-utils` [rename](https://github.com/atlashub/biostack-suits-docs/blob/master/fastx-utils/rename/rename.md)修改序列标签


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils cut -f1 \
	/project/dix-demo/data_analysis/zotu/unoise/align.txt | \
	/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils subseq /project/dix-demo/data_analysis/zotu/unoise/denoise.fasta - | \
	/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils rename - ZOTU \
	>/project/dix-demo/data_analysis/zotu/report/zotus.fasta
```


### 10. 构建ZOTU表

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/zotu/otutab/
```

**ZOTU丰度 :**

使用`usearch` [otutab](http://www.drive5.com/usearch/manual/cmd_otutab.html) 构建`OTU`表

```bash
cat /project/dix-demo/data_analysis/zotu/labels/*.fastq                        \
   > /project/dix-demo/data_analysis/zotu/reads/striped.fastq;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch               \
    -otutab /project/dix-demo/data_analysis/zotu/reads/striped.fastq           \
    -zotus  /project/dix-demo/data_analysis/zotu/report/zotus.fasta            \
    -strand plus                                                               \
    -id 0.97                                                                   \
    -otutabout /project/dix-demo/data_analysis/zotu/otutab/zotu_table.txt      \
    -mapout /project/dix-demo/data_analysis/zotu/otutab/A-1.map.txt            \
    -threads 72; 
```

`usearch -otutab`参数解析：

    -zotus          指定数据库文件
    -strand         指定plus或者both
    -id             最小的序列一致性
    -otutabout      生成OTU表
    -mapout         生成序列标签和OTU标签表
    -threads        线程数

`usearch -otutab`命令参考：[http://www.drive5.com/usearch/manual/cmd_otutab.html](http://www.drive5.com/usearch/manual/cmd_otutab.html)


如果因为`USEARCH 32bit` 受到限制,可以使用 `atlas-utils` [rare](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/rare/rare.md)  或者 分拆样本模式:


**样本分拆 :**

使用`usearch` [otutab](http://www.drive5.com/usearch/manual/cmd_otutab.html) 构建`OTU`表

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch              \
    -otutab  /project/dix-demo/data_analysis/zotu/labels/A-1.fastq            \
    -zotus   /project/dix-demo/data_analysis/zotu/report/zotus.fasta          \
    -strand  plus                                                             \
    -id 0.97                                                                  \
    -otutabout /project/dix-demo/data_analysis/zotu/otutab/A-1.zotu_table.txt \
    -mapout /project/dix-demo/data_analysis/zotu/otutab/A-1.map.txt           \
    -threads 72; 
```

可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
   grep -v "#" | \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch -otutab  /project/dix-demo/data_analysis/zotu/labels/$F[0].fastq -zotus /project/dix-demo/data_analysis/zotu/report/zotus.fasta -strand plus -id 0.97 -otutabout /project/dix-demo/data_analysis/zotu/otutab/$F[0].zotu_table.txt -mapout /project/dix-demo/data_analysis/zotu/otutab/$F[0].map.txt -threads 72;\n}' | \
   bash
```

**合并ZOTU表 :**

使用 `tsv-utils` [join](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/join/join.md) 合并各个样本匹配的`zotu_table.txt`文件，将表头的`catalog`替换为`OTU ID`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils join -p 0 /project/dix-demo/data_analysis/zotu/otutab/A-1.zotu_table.txt \
/project/dix-demo/data_analysis/zotu/otutab/A-2.zotu_table.txt \
/project/dix-demo/data_analysis/zotu/otutab/A-3.zotu_table.txt \
/project/dix-demo/data_analysis/zotu/otutab/B-1.zotu_table.txt \
/project/dix-demo/data_analysis/zotu/otutab/B-2.zotu_table.txt \
/project/dix-demo/data_analysis/zotu/otutab/B-3.zotu_table.txt | \
sed 's/catalog/OTU ID/'>/project/dix-demo/data_analysis/zotu/unoise/zotu_table.txt
```

**ZOTU表重新排序 :**

使用`fastx-utils` [view](https://github.com/atlashub/biostack-suits-docs/blob/master/fastx-utils/view/view.md) 获取 `zotu` 的顺序并使用 `tsv-utils` [reorder](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/reorder/reorder.md) 将`zotu_table.txt`进行重排

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils view \
/project/dix-demo/data_analysis/zotu/report/zotus.fasta |          \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils reorder \
/project/dix-demo/data_analysis/zotu/unoise/zotu_table.txt -       \
>/project/dix-demo/data_analysis/zotu/report/zotu_table.txt 
```

**转换相对丰度 :**

使用`atlas-utils` [counts2freqs](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/counts2freqs/counts2freqs.md) 转换`zotu`表的计数转为相对丰度;

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils counts2freqs /project/dix-demo/data_analysis/zotu/report/zotu_table.txt \
>/project/dix-demo/data_analysis/zotu/report/zotu_table_freqs.txt
```

**统计ZOTU表信息 :**

使用`usearch` [otutab_stats](https://www.drive5.com/usearch/manual/cmd_otutab_stats.html)生成`zotu`表的汇总信息

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch  -otutab_stats /project/dix-demo/data_analysis/zotu/report/zotu_table.txt \
-output /project/dix-demo/data_analysis/zotu/report/zotu_report.txt
```


### 11. 物种分类

**创建目录 :**

```bash
mkdir -p  /project/dix-demo/data_analysis/classify/classify
mkdir -p  /project/dix-demo/data_analysis/classify/zotus
```

** 物种分类 :**
使用`usearch` [sintax](https://drive5.com/usearch/manual/cmd_sintax.html) 根据`rdp_16s_v16_sp`进行预测分类水平

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -sintax      /project/dix-demo/data_analysis/zotu/report/zotus.fasta      \
    -db /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/rdp_16s_v16_sp.udb \
    -strand plus                                                              \
    -tabbedout /project/dix-demo/data_analysis/classify/classify/classify.txt \
    -log /project/dix-demo/data_analysis/classify/classify/classify.log       \
    -sintax_cutoff 0.8                                                        \
    -threads 72
```

`usearch -sintax`参数解析：  

    -db                 指定数据库（输入）
    -strand             指定plus或者both
    -tabbedout          物种分类预测结果列表（输出）
    -log                日志文件（输出）
    -sintax_cutoff      指定物种分类的预测准确性阈值
    -threads            线程数


`usearch -sintax`命令参考：[http://www.drive5.com/usearch/manual/cmd_sintax.html](http://www.drive5.com/usearch/manual/cmd_sintax.html)


**过滤ZOTU表 :**

使用`atlas-utils`  [filter](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/filter/filter.md) 将`classify/classify/classify.txt`的`OTU ID`和分类信息分别导出到`zotu_identifiers.txt`和`classify.txt`中

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils filter -r -t NONE  \
     /project/dix-demo/data_analysis/classify/classify/classify.txt                     \
  2> /project/dix-demo/data_analysis/classify/zotus/zotu_identifiers.txt                \
  1> /project/dix-demo/data_analysis/classify/zotus/classify.txt ;
```


### 12. ZOTU 注释


**抽取序列:**

使用`fastx-utils` [subseq](https://github.com/atlashub/biostack-suits-docs/blob/master/fastx-utils/subseq/subseq.md) 根据 `zotu_identifiers.txt` 提供的序列`ID`对`zotus.fasta`中的序列取交集

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils subseq  \
     /project/dix-demo/data_analysis/zotu/report/zotus.fasta                 \
     /project/dix-demo/data_analysis/classify/zotus/zotu_identifiers.txt     \
     >/project/dix-demo/data_analysis/classify/zotus/zotus.fasta
```


**抽取zotu表:**

使用`tsv-utils` [subset](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/subset/subset.md) 根据 `zotu_identifiers.txt`提供的序列`ID` 对 `zotu_table.txt`中的列表取交集

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subset  \
    /project/dix-demo/data_analysis/zotu/report/zotu_table.txt             \
    /project/dix-demo/data_analysis/classify/zotus/zotu_identifiers.txt    \
    >/project/dix-demo/data_analysis/classify/zotus/zotu_table.txt
```

**转换相对丰度:** 

使用`atlas-utils` [counts2freqs](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/counts2freqs/counts2freqs.md) 将`czotu_table.txt`的计数表转为相对丰度表

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  counts2freqs \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table.txt                  \
    >/project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs.txt
```


**zotu表注释:**

提取`classify.txt` 的序列`ID`和物种注释传递给 `atlas-utils` [annotation](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/annotation/annotation.md) 对 `zotu_table.txt` 进行注释

```bash
cut -f1,4 /project/dix-demo/data_analysis/classify/zotus/classify.txt | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils annotation - /project/dix-demo/data_analysis/classify/zotus/zotu_table.txt  \
    >/project/dix-demo/data_analysis/classify/zotus/zotu_table_ann.txt
```

**zotu表注释:**

提取`classify.txt`的序列`ID`和物种注释传递给`atlas-utils annotation`对`zotu_table_freqs.txt`进行注释

```bash
cut -f1,4 /project/dix-demo/data_analysis/classify/zotus/classify.txt | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils annotation - /project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs.txt \
    >/project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs_ann.txt ;
```

**格式转换:**

使用`tsv-utils` [tsv2xlsx](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/tsv2xlsx/tsv2xlsx.md) 将`zotu_table.txt` ，`zotu_table_ann.txt`，`zotu_table_freqs.txt`，`zotu_table_freqs_ann.txt`合并为`classify/zotus/zotu_table.xlsx` 

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils tsv2xlsx   /project/dix-demo/data_analysis/classify/zotus/zotu_table.xlsx \
otu_table:/project/dix-demo/data_analysis/classify/zotus/zotu_table.txt \
otu_table_ann:/project/dix-demo/data_analysis/classify/zotus/zotu_table_ann.txt \
otu_table_freqs:/project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs.txt \
otu_table_freqs_ann:/project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs_ann.txt
```


### 13. biom文件格式转换

**创建目录 :**

```bash
mkdir -p  /project/dix-demo/data_analysis/classify/report
```

使用[biom](http://biom-format.org/) `convert`转化 `zotu_table_ann.txt` 为 `zotu_table.biom` 格式

```bash
biom convert    -i /project/dix-demo/data_analysis/classify/zotus/zotu_table_ann.txt   \
                -o /project/dix-demo/data_analysis/classify/zotus/zotu_table.biom      \
                --table-type "OTU table"                                               \
                --to-json                                                              \
                --process-obs-metadata taxonomy
```

使用[biom](http://biom-format.org/) `summarize-table`对`zotu_table.biom`进行统计摘要

```bash
biom summarize-table  -i /project/dix-demo/data_analysis/classify/zotus/zotu_table.biom \
                      -o /project/dix-demo/data_analysis/classify/zotus/zotu_summary.txt
```

`biom`参数解析：

    convert                 文本表格与biom互转
    summarize-table         统计摘要
    -i                      输入文件
    -o                      输出文件
    --table-type            
    --to-json               转换经典表格为JSON格式
    --process-obs-metadata  带物种注释表格互转

`biom` 命令参考：https://blog.csdn.net/woodcorpse/article/details/84678543


**统计信息 :**

使用`atlas-utils` [summary](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/summary/summary.md) 计算`zotu_table.txt`的 `(Z)OTU` 数目以及`Tags`数目并使用`tsv-utils` [transpose](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/transpose/transpose.md) 进行转置

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils summary       \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table.txt |                \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils transpose - \
    >/project/dix-demo/data_analysis/classify/report/zotu.stats.txt
```


### 14. 标准化ZOTU表

**创建目录 :**

```bash
mkdir -p  /project/dix-demo/data_analysis/classify/partition
mkdir -p  /project/dix-demo/data_analysis/classify/report
```

** 获取样本中最小reads数 :**

使用`atlas-utils` [min_size](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/min_size/min_size.md) 获取 `zotu_table.txt` 的最小值

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils min_size \
/project/dix-demo/data_analysis/zotu/report/zotu_table.txt
```

    4220


**使用 USEARCH 标准化ZOTU表 :**

使用`usearch` [otutab_rare](http://www.drive5.com/usearch/manual/cmd_otutab_rare.html) 对`zotu_table.txt`进行标准化

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch  -otutab_rare  \
    /project/dix-demo/data_analysis/zotu/report/zotu_table.txt                  \
    -sample_size  4220                                                          \
    -randseed 11                                                                \
    -output /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.txt
```

`usearch -otutab_rare`参数解析

    -sample_size        样本读取次数
    -randseed           随机数种子
    -output             输出标准化文件（输出）

`usearch -otutab_rare`命令参考：[http://www.drive5.com/usearch/manual/cmd_otutab_rare.html](http://www.drive5.com/usearch/manual/cmd_otutab_rare.html)


如果因为`USEARCH 32bit` 受到限制,可以使用 `atlas-utils rare`  或者 分拆样本模式:

**使用 Atlas-utils 标准化ZOTU表 :**


使用`atlas-utils`  [rare](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/rare/rare.md) 对`zotu_table.txt`进行标准化

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils rare  -s 11 \
/project/dix-demo/data_analysis/zotu/report/zotu_table.txt 4220 \
>/project/dix-demo/data_analysis/classify/partition/zotu_table_norm.txt;
```

**分拆样本 :**


分拆: 使用`tsv-utils` [subcolumn](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/subcolumn/subcolumn.md) 提取 `#OTU ID` 和 `A-1` 列

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subcolumn -k \
/project/dix-demo/data_analysis/zotu/report/zotu_table.txt "#OTU ID",A-1    \
>/project/dix-demo/data_analysis/classify/partition/zotu_table.A-1.txt;
```

可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | grep -v "#" | perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subcolumn -k /project/dix-demo/data_analysis/zotu/report/zotu_table.txt "#OTU ID",$F[0] > /project/dix-demo/data_analysis/classify/partition/zotu_table.$F[0].txt;\n}' |bash
```


抽样: 使用`atlas-utils` [rare](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/rare/rare.md) 对`zotu_table.A-1.txt`进行标准化

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils rare  -s 11 \
/project/dix-demo/data_analysis/classify/partition/zotu_table.A-1.txt  4220     \
>/project/dix-demo/data_analysis/classify/partition/zotu_table_norm.A-1.txt;
```

可使用 `bash` 进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
   grep -v "#" | \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch  -otutab_rare  /project/dix-demo/data_analysis/classify/partition/zotu_table.$F[0].txt  -sample_size  4220  -randseed 11  -output /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.$F[0].txt;\n}' | \
   bash
```


抽取: 使用`fastx-utils` [view](https://github.com/atlashub/biostack-suits-docs/blob/master/fastx-utils/view/view.md)  抽取 `zotus.fasta` 的序列`ID`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/fastx-utils view  \
   /project/dix-demo/data_analysis/classify/zotus/zotus.fasta              \
   >/project/dix-demo/data_analysis/classify/zotus/zotus.txt
```


合并: 使用`tsv-utils` [reorder](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/reorder/reorder.md) 根据 `zotus.txt` 进行排序，然后使用 `sed` 将 `catalog` 替换为 `OTU ID`


```bash
cut -f1  /project/dix-demo/mapping_file.txt | \
   grep -v "#" |                              \
   perl -lne 'print qq{/project/dix-demo/data_analysis/classify/partition/zotu_table_norm.$_.txt}' | \
   paste -s |                                 \
   perl -lne 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils join -p 0 $_ | /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils reorder  - /project/dix-demo/data_analysis/classify/zotus/zotus.txt  |sed "s/catalog/OTU ID/">/project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt;\n}' |           \
bash
```

等价于:

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils join -p 0        \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.A-1.txt   \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.A-2.txt   \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.A-3.txt   \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.B-1.txt   \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.B-2.txt   \
       /project/dix-demo/data_analysis/classify/partition/zotu_table_norm.B-3.txt | \
   /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils reorder  -    \
      /project/dix-demo/data_analysis/classify/zotus/zotus.txt  |                   \
   sed 's/catalog/OTU ID/'  |                                                       \
   >/project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt
```


### 15. 统计样本信息


使用`tsv-utils` [join](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/join/join.md) 合并`trimming.stats.txt`，`mergepairs.stats.txt`，`pcrsearch.stats.txt`，`zotu.stats.txt` 并使用`tsv-utils cut`删除第`6,9`列，使用`tsv-utils` [add_headline](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/add_headline/add_headline.md)  添加标题行`"#label\tTrimmomatic\t\t\t\tmergepairs\t\tprimer_match\t\t\t\t\tZotu\t\t"`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils join      \
    /project/dix-demo/data_analysis/trimming/report/trimming.stats.txt       \
    /project/dix-demo/data_analysis/mergepairs/report/mergepairs.stats.txt   \
    /project/dix-demo/data_analysis/primer_strip/report/pcrsearch.stats.txt  \
    /project/dix-demo/data_analysis/classify/report/zotu.stats.txt  |        \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils cut -d -f6,9   -  |  \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils add_headline  "#label\tTrimmomatic\t\t\t\tmergepairs\t\tprimer_match\t\t\t\t\tZotu\t\t"  -   \
    >/project/dix-demo/data_analysis/classify/report/sample.stats.long.txt
```


简化: 使用`cut`提取`sample.stats.long.txt`第`1,2,6,7,9,14`的数据并使用`grep -v`选取不含`catalog`的行

```bash
cut -f1,2,6,7,9,14  \
    /project/dix-demo/data_analysis/classify/report/sample.stats.long.txt  | \
grep -v "catalog" \
    >/project/dix-demo/data_analysis/classify/report/sample.stats.short.txt
```


### 16. 不同物种分类水平的相对丰度计算

**创建目录 :**

```bash
mkdir -p  /project/dix-demo/data_analysis/taxonomy/classify
```

**计算不同水平的物种组成丰度 :**

使用`atlas-utils` [level](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/level/level.md) 计算`zotu_table_ann.txt`的门水平丰度并使用`tsv-utils view`转化为整数，使用`atlas-utils` [counts2freqs](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/counts2freqs/counts2freqs.md)  将`phylum.counts.txt`的计数转为相对丰度


```bash
declare -A taxonomy=( ["phylum"]="p" ["order"]="o" ["class"]="c" ["family"]="f" ["genus"]="g");
for i in phylum order class family genus
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils level -l   ${taxonomy[${i}]} /project/dix-demo/data_analysis/classify/zotus/zotu_table_ann.txt | \
      /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils view -r - \
      >/project/dix-demo/data_analysis/taxonomy/classify/${i}.counts.txt;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  counts2freqs \
      /project/dix-demo/data_analysis/taxonomy/classify/${i}.counts.txt            \
      >/project/dix-demo/data_analysis/taxonomy/classify/${i}.freqs.txt;
done
```

**合并文件 :**


将门纲目科属所有的`counts.txt` 和`freqs.txt` 文件合并为`taxonomy/classify/taxonomy.xlsx` 

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils tsv2xlsx           \
    /project/dix-demo/data_analysis/taxonomy/classify/taxonomy.xlsx                   \
    phylum.counts:/project/dix-demo/data_analysis/taxonomy/classify/phylum.counts.txt \
    phylum.freqs:/project/dix-demo/data_analysis/taxonomy/classify/phylum.freqs.txt   \
    order.counts:/project/dix-demo/data_analysis/taxonomy/classify/order.counts.txt   \
    order.freqs:/project/dix-demo/data_analysis/taxonomy/classify/order.freqs.txt     \
    class.counts:/project/dix-demo/data_analysis/taxonomy/classify/class.counts.txt   \
    class.freqs:/project/dix-demo/data_analysis/taxonomy/classify/class.freqs.txt     \
    family.counts:/project/dix-demo/data_analysis/taxonomy/classify/family.counts.txt \
    family.freqs:/project/dix-demo/data_analysis/taxonomy/classify/family.freqs.txt   \
    genus.counts:/project/dix-demo/data_analysis/taxonomy/classify/genus.counts.txt   \
    genus.freqs:/project/dix-demo/data_analysis/taxonomy/classify/genus.freqs.txt  ;
```


### 17. krona 可视化

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/taxonomy/krona/
```

**创建krona文件 :**

使用`atlas-utils` [krona](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/krona/krona.md) 生成样品`A-1`的`krona`兼容模式


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils krona       \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs_ann.txt A-1  \
    >/project/dix-demo/data_analysis/taxonomy/krona/A-1.txt
```

可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
  grep -v "#" | \
  perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils krona /project/dix-demo/data_analysis/classify/zotus/zotu_table_freqs_ann.txt $F[0] >/project/dix-demo/data_analysis/taxonomy/krona/$F[0].txt ;\n}' | \
bash
```

使用`tImportText`将`taxonomy/krona/`目录下的所有`.txt`文件生成多图层交互式饼图

```bash
ktImportText -o /project/dix-demo/data_analysis/taxonomy/krona/krona.html  \
    /project/dix-demo/data_analysis/taxonomy/krona/A-1.txt \
    /project/dix-demo/data_analysis/taxonomy/krona/A-2.txt \
    /project/dix-demo/data_analysis/taxonomy/krona/A-3.txt \
    /project/dix-demo/data_analysis/taxonomy/krona/B-1.txt \
    /project/dix-demo/data_analysis/taxonomy/krona/B-2.txt \
    /project/dix-demo/data_analysis/taxonomy/krona/B-3.txt 
```


### 18. 不同物种分类水平可视化

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/taxonomy/bars
mkdir -p /project/dix-demo/data_analysis/taxonomy/heatmap
```

**绘制柱状图 :**

使用`atlas-utils` [rank](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/rank/rank.md) 对 `freqs.txt`进行排序并显示指定行,使用`barplot.R`对`*.freqs.txt`绘制百分比柱状图,使用`pdf2png`导出`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils rank -r 10 -m -a \
    /project/dix-demo/data_analysis/taxonomy/classify/phylum.freqs.txt                \
    >/project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/barplot.R                       \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt                 \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.stack.pdf phylum ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                         \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.stack.pdf;
```

可使用`bash`进行循环迭代：

```bash
declare -A taxonomy=( ["phylum"]="10" ["order"]="15" ["class"]="20" ["family"]="25" ["genus"]="30");
for i in phylum order class family genus;
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils rank  \
    -r ${taxonomy[${i}]} -m -a \
    /project/dix-demo/data_analysis/taxonomy/classify/${i}.freqs.txt \
    >/project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.freqs.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/barplot.R \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.freqs.txt  \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.stack.pdf ${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.stack.pdf;
done
```

**格式化元数据 :**

提取`mapping_file.txt`的`1，2`列信息，使用`grep`选取不含`#`的行

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt | \
grep -v  "#">/project/dix-demo/data_analysis/taxonomy/bars/metadata.txt
```

**绘制热图 :**

使用`heatmap.R`绘制`*.freqs.txt`的热图并导出`pdf`格式，使用`pdf2png`导出`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/heatmap.R            \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt      \
    /project/dix-demo/data_analysis/taxonomy/heatmap/phylum.10.heatmap.pdf phylum;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png              \
    /project/dix-demo/data_analysis/taxonomy/heatmap/phylum.10.heatmap.pdf;
```

可使用`bash`进行循环迭代：


```bash
declare -A taxonomy=( ["phylum"]="10" ["order"]="15" ["class"]="20" ["family"]="25" ["genus"]="30");
for i in phylum order class family genus;
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/heatmap.R                         \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.freqs.txt      \
    /project/dix-demo/data_analysis/taxonomy/heatmap/${i}.${taxonomy[${i}]}.heatmap.pdf ${i};
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                           \
    /project/dix-demo/data_analysis/taxonomy/heatmap/${i}.${taxonomy[${i}]}.heatmap.pdf;
done
```

**绘制统计图 :**

使用`stats.R` 绘制 `*.freqs.txt`的百分比柱状图并导出`pdf`格式,使用`pdf2png`导出`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/stats.R -e T          \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt       \
    /project/dix-demo/data_analysis/taxonomy/bars/metadata.txt              \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.average.pdf phylum ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png               \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.average.pdf;
```

可使用`bash`进行循环迭代：

```bash
declare -A taxonomy=( ["phylum"]="10" ["order"]="15" ["class"]="20" ["family"]="25" ["genus"]="30");
for i in phylum order class family genus;
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/stats.R -e T                  \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.freqs.txt  \
    /project/dix-demo/data_analysis/taxonomy/bars/metadata.txt                      \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.average.pdf ${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                       \
    /project/dix-demo/data_analysis/taxonomy/bars/${i}.${taxonomy[${i}]}.average.pdf;
done
```


### 19. 构建系统进化树

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/phylogeny/clust/;
mkdir -p /project/dix-demo/data_analysis/phylogeny/report/;
```

**距离矩阵构建树 :**

使用`usearch`  [calc_distmx](https://drive5.com/usearch/manual/cmd_calc_distmx.html)  生成`zotus.fasta`的稀疏距离矩阵

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch -calc_distmx \
    /project/dix-demo/data_analysis/classify/zotus/zotus.fasta                \
    -tabbedout /project/dix-demo/data_analysis/phylogeny/clust/distmx.txt     \
    -maxdist 0.2  \
    -termdist 0.3 
```

`usearch -calc_distmx`参数解析

    -tabbedout        输出矩阵列表（输出）
    -maxdist          最大距离
    -termdist         终止计算的序列一致度阈值

`usearch -calc_distmx`命令参考：[http://www.drive5.com/usearch/manual/cmd_calc_distmx.html](http://www.drive5.com/usearch/manual/cmd_calc_distmx.html)


使用`usearch` [cluster_aggd](https://drive5.com/usearch/manual/cmd_cluster_aggd.html) 对距离矩阵 `distmx.txt` 输出`Newick`格式的进化树和簇文件

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -cluster_aggd /project/dix-demo/data_analysis/phylogeny/clust/distmx.txt \
    -treeout /project/dix-demo/data_analysis/phylogeny/report/zotus.tree \
    -clusterout /project/dix-demo/data_analysis/phylogeny/clust/clusters.txt \
    -id 0.80 \
    -linkage max;
```

`usearch -cluster_aggd`参数解析

    -treeout        输出Newick格式的进化树（输出）
    -clusterout     簇文件（输出）
    -id             序列一致度阈值
    -linkage        凝聚链接类型max（默认）、min 或者 avg.

`usearch -cluster_aggd`命令参考：[http://www.drive5.com/usearch/manual/cmd_cluster_aggd.html](http://www.drive5.com/usearch/manual/cmd_cluster_aggd.html)



### 20. alpha多样性

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/alpha/partition/
mkdir -p /project/dix-demo/data_analysis/alpha/diversity/
mkdir -p /project/dix-demo/data_analysis/alpha/report/
```

**计算alpha多样性指数 :**

使用`usearch` [alpha_div](http://www.drive5.com/usearch/manual/cmd_alpha_div.html) 对 `zotu_table.A-1.txt` 计算 `alpha`多样性

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -alpha_div /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt   \
    -output /project/dix-demo/data_analysis/alpha/partition/alpha.txt
```

`usearch` [alpha_div](https://drive5.com/usearch/manual/cmd_alpha_div.html) 参数解析

    -output        输出alpha文件（输出）

`usearch -alpha_div`命令参考：[http://www.drive5.com/usearch/manual/cmd_alpha_div.html](http://www.drive5.com/usearch/manual/cmd_alpha_div.html)

如果因为`USEARCH 32bit` 受到限制,可以使用下面模式:


**分拆样本 :**

使用`tsv-utils` [subcolumn](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/subcolumn/subcolumn.md) 提取`zotu_table.txt`的`"#OTU ID",A-1`列
使用`usearch` [alpha_div](http://www.drive5.com/usearch/manual/cmd_alpha_div.html) 对 `zotu_table.A-1.txt`计算`alpha`多样性

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subcolumn -k             /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt "#OTU ID",A-1        >/project/dix-demo/data_analysis/alpha/partition/zotu_table.A-1.txt;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch \
    -alpha_div /project/dix-demo/data_analysis/alpha/partition/zotu_table.A-1.txt \
    -output /project/dix-demo/data_analysis/alpha/partition/alpha.A-1.txt
```

可使用`bash`进行循环迭代：

```bash
cut -f1 /project/dix-demo/mapping_file.txt  | \
   grep -v "#" | \
   perl -ane 'print qq{/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subcolumn -k /project/dix-demo/data_analysis/classify/zotus/zotu_table.txt "#OTU ID",$F[0] >/project/dix-demo/data_analysis/alpha/partition/zotu_table.$F[0].txt;\n}' |\
   bash
```



**合并表格 alpha 多样性 :**

合并`alpha/partition/`中所有的`alpha`多样性文件，使用`tsv-utils` [view](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/view/view.md) 在第一行添加`#`，选取不以`Sample`开头的行然后根据`mapping_file.txt`的样品名顺序进行排序

```bash
cat /project/dix-demo/data_analysis/alpha/partition/alpha*.txt | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils view -c - | \
    grep -P -v "^Sample" | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils reorder - \
    <(cut -f1 /project/dix-demo/mapping_file.txt)  \
    >/project/dix-demo/data_analysis/alpha/diversity/alpha.long.txt
```

使用`tsv-utils`  [subcolumn](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/subcolumn/subcolumn.md)  选取`alpha.long.txt` 中`#Sample,richness,chao1,shannon_2,simpson,dominance,equitability`的列，重组alpha多样性文件

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils subcolumn \
   -k /project/dix-demo/data_analysis/alpha/diversity/alpha.long.txt         \
   "#Sample,richness,chao1,shannon_2,simpson,dominance,equitability"         \
   >/project/dix-demo/data_analysis/alpha/diversity/alpha.txt
```


### 21. 稀释曲线

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/alpha/rarefaction
```

**格式化zotu表 :**

删除`zotu_table.txt`中的#，并导出为`alpha/rarefaction/zotu_table.txt`

```bash
sed 's/#//' /project/dix-demo/data_analysis/classify/zotus/zotu_table.txt \
    >/project/dix-demo/data_analysis/alpha/rarefaction/zotu_table.txt
```

**格式化元数据 :**

提取`mapping_file.txt`的`1，2`列信息，使用`grep`选取不含`#`的行

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt | \
    grep -v '#' \
    >/project/dix-demo/data_analysis/alpha/report/metadata.txt
```

**构建稀释数据值 :**

使用`rtk`做稀释曲线。

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/rtk memory      \
    -i  /project/dix-demo/data_analysis/alpha/rarefaction/zotu_table.txt \
    -o  /project/dix-demo/data_analysis/alpha/rarefaction/rare.          \
    -ns                                                                  \
    -t 72                                                                \
    -d 500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000,10500,11000,11500,12000,12500,13000,13500,14000,14500,15000,15500,16000,16500,17000,17500,18000,18500,19000,19500,20000,20500,21000,21500,22000,22500,23000,23500,24000,24500,25000,25500,26000,26500,27000,27500,28000,28500,29000,29500,30000,30500,31000,31500,32000,32500,33000,33500,34000,34500,35000,35500,36000,36500,37000,37500,38000,38500,39000,39500,40000,40500,41000,41500,42000,42500,43000,43500,44000,44500,45000,45500,46000,46500,47000,47500,48000,48500,49000,49500,50000 \
    -r 50 ;
```

根据[rtk](https://rdrr.io/cran/rtk/man/rtk.html) 获得的结果绘制`richness`指数的稀释曲线，并导出`png`格式


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/rarefaction_curve.R -t T \
    /project/dix-demo/data_analysis/alpha/rarefaction/rare.alpha_richness.tsv  \
    /project/dix-demo/data_analysis/alpha/report/metadata.txt                  \
    /project/dix-demo/data_analysis/alpha/rarefaction/richness.rarefactions_curve.pdf richness ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/alpha/rarefaction/richness.rarefactions_curve.pdf ;
```


可使用`bash`进行循环迭代：

```bash
for i in richness shannon chao1 simpson ;
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/rarefaction_curve.R -t T \
    /project/dix-demo/data_analysis/alpha/rarefaction/rare.alpha_${i}.tsv      \
    /project/dix-demo/data_analysis/alpha/report/metadata.txt                  \
    /project/dix-demo/data_analysis/alpha/rarefaction/${i}.rarefactions_curve.pdf ${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/alpha/rarefaction/${i}.rarefactions_curve.pdf ;
done
```


### 22 beta多样性

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/beta/distmx/
mkdir -p /project/dix-demo/data_analysis/beta/report/
```

**计算beta多样性指数 :**

使用`usearch` [beta_div](https://drive5.com/usearch/manual/cmd_beta_div.html) 根据`zotu_table_norm.txt`构建`beta`多样性矩阵

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch      \
    -beta_div /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt  \
    -filename_prefix /project/dix-demo/data_analysis/beta/distmx/     \
    -tree /project/dix-demo/data_analysis/phylogeny/report/zotus.tree \
    -metrics jaccard,bray_curtis,euclidean,unifrac,unifrac_binary     \
    -log /project/dix-demo/data_analysis/beta/distmx/beta_div.log
```

`usearch -beta_div`参数解析

    -filename_prefix      输出文件的前缀, 可以指定为目录;
    -tree                 系统进化树;
    -metrics              beta多样性指数, 以`,`分割;
    -log                  日志文件;


`usearch -beta_div`命令参考：[https://drive5.com/usearch/manual/cmd_beta_div.html](https://drive5.com/usearch/manual/cmd_beta_div.html)


**格式化元数据 :**

获取`mapping_file.txt`的`1,2`列选取不含`#`的列

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt | \
grep -v '#'>/project/dix-demo/data_analysis/beta/report/metadata.txt
```

**绘制 upgma 与柱状图 :**

- 创建`beta/upgma/jaccard`目录
- 使用`distmx-fmt`转换`jaccard.sorted.txt`为`jaccard.txt`
- 使用`upgma.R`根据`phylum.10.freqs.txt`、`jaccard.txt`、`metadata.txt`绘制`jaccard`树和门水平百分比水平柱状图
- 转化为`png`格式

```bash
mkdir -p /project/dix-demo/data_analysis/beta/upgma/jaccard;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/distmx-fmt       \
    /project/dix-demo/data_analysis/beta/distmx/jaccard.sorted.txt     \
    >/project/dix-demo/data_analysis/beta/upgma/jaccard/jaccard.txt; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/upgma.R          \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt  \
    /project/dix-demo/data_analysis/beta/upgma/jaccard/jaccard.txt     \
    /project/dix-demo/data_analysis/beta/report/metadata.txt           \
    /project/dix-demo/data_analysis/beta/upgma/jaccard/jaccard.upgma.bar.pdf jaccard ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-demo/data_analysis/beta/upgma/jaccard/jaccard.upgma.bar.pdf ;
```

可使用`bash`进行循环迭代：


```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary
do
mkdir -p /project/dix-demo/data_analysis/beta/upgma/${i};
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/distmx-fmt       \
    /project/dix-demo/data_analysis/beta/distmx/${i}.sorted.txt        \
    >/project/dix-demo/data_analysis/beta/upgma/${i}/${i}.txt; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/upgma.R          \
    /project/dix-demo/data_analysis/taxonomy/bars/phylum.10.freqs.txt  \
    /project/dix-demo/data_analysis/beta/upgma/${i}/${i}.txt           \
    /project/dix-demo/data_analysis/beta/report/metadata.txt           \
    /project/dix-demo/data_analysis/beta/upgma/${i}/${i}.upgma.bar.pdf ${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-demo/data_analysis/beta/upgma/${i}/${i}.upgma.bar.pdf ;
done
```

**绘制 PCA 图 :**


- 创建`pca`目录
- 使用`PCA.R`根据`zotu_table_norm.txt`进行主成分分析
- 转换`pdf`为`png`格式

```bash
mkdir -p /project/dix-demo/data_analysis/beta/pca ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T  \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt \
    /project/dix-demo/data_analysis/beta/report/metadata.txt           \
    /project/dix-demo/data_analysis/beta/pca/zotu.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-demo/data_analysis/beta/pca/zotu.pca.pdf;
```

**绘制 PCoA 图 :**


- 创建目录
- 根据距离文件进行主坐标分析
- 转换`pdf`格式为`png`

```bash
mkdir -p /project/dix-seq/data_analysis/beta/pcoa/jaccard ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T \
    /project/dix-seq/data_analysis/beta/distmx/jaccard.txt             \
    /project/dix-seq/data_analysis/beta/report/metadata.txt            \
    /project/dix-seq/data_analysis/beta/pcoa/jaccard/jaccard.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-seq/data_analysis/beta/pcoa/jaccard/jaccard.pcoa.pdf ;
```

可使用`for`进行循环迭代：

```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary
do
mkdir -p /project/dix-seq/data_analysis/beta/pcoa/${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T \
    /project/dix-seq/data_analysis/beta/distmx/${i}.txt                \
    /project/dix-seq/data_analysis/beta/report/metadata.txt            \
    /project/dix-seq/data_analysis/beta/pcoa/${i}/${i}.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-seq/data_analysis/beta/pcoa/${i}/${i}.pcoa.pdf ;
done
```


**绘制 NMDS 图 :**

- 创建`beta/nmds/jaccard`目录
- 根据距离文件进行非度量多维尺度分析
- 转换`pdf格式`为`png格式`

```bash
 mkdir -p /project/dix-demo/data_analysis/beta/nmds/jaccard ; 
 /biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T \
    /project/dix-demo/data_analysis/beta/distmx/jaccard.txt             \
    /project/dix-seq/data_analysis/beta/report/metadata.txt             \
    /project/dix-demo/data_analysis/beta/nmds/jaccard/jaccard.nmds.pdf; 
 /biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-demo/data_analysis/beta/nmds/jaccard/jaccard.nmds.pdf ; 
```

可使用`for`进行循环迭代：

```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary; 
do 
mkdir -p /project/dix-demo/data_analysis/beta/nmds/${i} ; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T \
	/project/dix-demo/data_analysis/beta/distmx/${i}.txt                \
    /project/dix-seq/data_analysis/beta/report/metadata.txt             \
    /project/dix-demo/data_analysis/beta/nmds/${i}/${i}.nmds.pdf; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
	/project/dix-demo/data_analysis/beta/nmds/${i}/${i}.nmds.pdf ; 
done
```


### 23. 识别核心微生物群

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/core/report
mkdir -p /project/dix-demo/data_analysis/core/distmx
```

**计算距离矩阵 :**

使用`usearch`[calc_distmx](https://drive5.com/usearch/manual/cmd_calc_distmx.html) 根据`zotus.fasta`构建稀疏距离矩阵

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch              \
    -calc_distmx /project/dix-demo/data_analysis/classify/zotus/zotus.fasta   \
    -tabbedout /project/dix-demo/data_analysis/core/distmx/distmx.txt         \
    -maxdist 0.2                                                              \
    -termdist 0.3  ;
```

**删除未分类的序列 :**

```bash
perl -ane 'next if($#F != 3); print'                              \
     /project/dix-demo/data_analysis/classify/zotus/classify.txt  \
     >/project/dix-demo/data_analysis/core/report/classify.txt
```

**构建核心ZOTU :**

使用`usearch` [otutab_core](https://drive5.com/usearch/manual/cmd_otutab_core.html) 识别 `zotu_table_norm.txt` 中的核心微生物群


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/usearch                        \
    -otutab_core /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt     \
    -distmxin /project/dix-demo/data_analysis/core/distmx/distmx.txt                    \
    -sintaxin /project/dix-demo/data_analysis/core/report/classify.txt                  \
    -tabbedout /project/dix-demo/data_analysis/core/report/core.txt
```

`usearch -otutab_core`参数解析

    -distmxin       输入距离矩阵（输入）
    -sintaxin       输出核心OTU的物种分类（输出）
    -tabbedout      输出文件（输出）


`usearch -otutab_core`命令参考：[http://www.drive5.com/usearch/manual/cmd_alpha_div.html](http://www.drive5.com/usearch/manual/cmd_alpha_div.html)


### 24. rank_abundance曲线和物种累计曲线

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/rank_abundance/report
mkdir -p /project/dix-demo/data_analysis/specaccum_curve/report
```

**绘制rank_abundance 曲线:**

使用`rank_abundance.R`根据`zotu_table_norm.txt`绘制rank_abundance曲线，导出`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/rank_abundance.R   \
      /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt \
      /project/dix-demo/data_analysis/rank_abundance/report/rank_abundance.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png \
    /project/dix-demo/data_analysis/rank_abundance/report/rank_abundance.pdf ;
```

**绘制制物种累计曲线:**


使用`specaccum_curve.R`根据`zotu_table_norm.txt`绘制物种累计曲线，导出`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/specaccum_curve.R  \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt   \
    /project/dix-demo/data_analysis/specaccum_curve/report/specaccum_curve.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png \
    /project/dix-demo/data_analysis/specaccum_curve/report/specaccum_curve.pdf ;
```


### 25. anosim 组间差异

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/anosim/report
```

**计算anosim :**

- 使用`anosim-utils fmt`根据`mapping_file.txt`在jaccard bray_curtis euclidean unifrac unifrac_binary文件添加分组后缀；
- 使用`anosim.R`根据`jaccard.txt`输出不同组在jaccard bray_curtis euclidean unifrac unifrac_binary的组间差异；
- 转换`pdf`为`png`格式

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim-utils fmt  \
    /project/dix-demo/mapping_file.txt                                  \
    /project/dix-demo/data_analysis/beta/distmx/jaccard.txt             \
    >/project/dix-demo/data_analysis/anosim/report/jaccard.txt ; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim.R          \
    /project/dix-demo/data_analysis/anosim/report/jaccard.txt           \
    /project/dix-demo/data_analysis/anosim/report/jaccard jaccard ANOSIM;  
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png           \
    /project/dix-demo/data_analysis/anosim/report/jaccard.pdf; 
```

可使用`for`进行循环迭代：

```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary; 
do 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim-utils fmt \
    /project/dix-demo/mapping_file.txt                                 \
    /project/dix-demo/data_analysis/beta/distmx/${i}.txt               \
    >/project/dix-demo/data_analysis/anosim/report/${i}.txt ; 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim.R         \
    /project/dix-demo/data_analysis/anosim/report/${i}.txt             \
    /project/dix-demo/data_analysis/anosim/report/${i} ${i} ANOSIM;  
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png          \
    /project/dix-demo/data_analysis/anosim/report/${i}.pdf; 
done
```

**提取样本组之间样本编号:**


- 创建`anosim/subset/`目录;
- 创建`anosim/submatrix/`目录;
- 使用`anosim-utils paired` 根据`mapping_file.txt`取子矩阵

```bash
mkdir -p /project/dix-demo/data_analysis/anosim/subset/jaccard ;
mkdir -p /project/dix-demo/data_analysis/anosim/submatrix/jaccard ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim-utils paired \
    /project/dix-demo/mapping_file.txt                                    \
    /project/dix-demo/data_analysis/anosim/subset/jaccard ;
```

可使用`for`进行循环迭代：


```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary; 
do 
mkdir -p /project/dix-demo/data_analysis/anosim/subset/${i} ;
mkdir -p /project/dix-demo/data_analysis/anosim/submatrix/${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim-utils paired \
    /project/dix-demo/mapping_file.txt \
    /project/dix-demo/data_analysis/anosim/subset/${i} ;
done
```

**提取样本组之间距离矩阵:**


- 使用`tsv-utils` [submatrix](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/submatrix/submatrix.md) 获取分别获取`jaccard bray_curtis euclidean unifrac unifrac_binary`的在`A_B`组的子矩阵
- 使用`anosim.R`分别计算`A,B` `组间jaccard bray_curtis euclidean unifrac unifrac_binary`的显著性差异
- 转换`pdf`为`png`

`注意事项`：使用脚本实现如`A_B,A_C,B_C`等不同组间的循环，因`bash`太臃肿，故不在此展示。

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils submatrix \
    /project/dix-demo/data_analysis/anosim/report/jaccard.txt                \
    /project/dix-demo/data_analysis/anosim/subset/jaccard/A_B.txt            \
    >/project/dix-demo/data_analysis/anosim/submatrix/jaccard/A_B.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim.R               \
    /project/dix-demo/data_analysis/anosim/submatrix/jaccard/A_B.txt         \
    /project/dix-demo/data_analysis/anosim/submatrix/jaccard/A_B jaccard A_B;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                \
    /project/dix-demo/data_analysis/anosim/submatrix/jaccard/A_B.pdf ;
```

可使用`for`进行循环迭代：


```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary; 
do 
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils submatrix \
    /project/dix-demo/data_analysis/anosim/report/${i}.txt                   \
    /project/dix-demo/data_analysis/anosim/subset/${i}/A_B.txt               \
    >/project/dix-demo/data_analysis/anosim/submatrix/${i}/A_B.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/anosim.R               \
    /project/dix-demo/data_analysis/anosim/submatrix/${i}/A_B.txt            \
    /project/dix-demo/data_analysis/anosim/submatrix/${i}/A_B ${i} A_B;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png \
    /project/dix-demo/data_analysis/anosim/submatrix/${i}/A_B.pdf ;
done
```

**合并统计信息 :**

分别合并`jaccard bray_curtis euclidean unifrac unifrac_binary`的所有`.signif`文件

```bash
cat /project/dix-demo/data_analysis/anosim/submatrix/jaccard/*.signif    \
    >/project/dix-demo/data_analysis/anosim/report/jaccard.paired.signif;
```

可使用`for`进行循环迭代：


```bash
for i in jaccard bray_curtis euclidean unifrac unifrac_binary; 
do 
cat /project/dix-demo/data_analysis/anosim/submatrix/${i}/*.signif    \
    >/project/dix-demo/data_analysis/anosim/report/${i}.paired.signif;
done
```


### 26. DESeq2 差异分析

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/DESeq2/counts
mkdir -p /project/dix-demo/data_analysis/DESeq2/report
mkdir -p /project/dix-demo/data_analysis/DESeq2/stats
```

**格式化元信息 :**

提取`mapping_file.txt`中的1,2列并分别为`phylum class order family genus.counts.txt`文件添加分组表头

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt |                                 \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils groupline - \
    /project/dix-demo/data_analysis/taxonomy/classify/phylum.counts.txt        \
    >/project/dix-demo/data_analysis/DESeq2/counts/phylum.counts.txt ;
```

可使用`for`进行循环迭代：

```bash
for i in phylum class order family genus; 
do
cut -f1,2 /project/dix-demo/mapping_file.txt |                                 \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils groupline - \
      /project/dix-demo/data_analysis/taxonomy/classify/${i}.counts.txt        \
    >/project/dix-demo/data_analysis/DESeq2/counts/${i}.counts.txt ;
done
```

**分组差异分析 :**

- 分别在`phylum class order family genus`和不同分组之间实现循环
- 创建`DESeq2/DESeq2/A_B`目录;
- 使用`DESeq2.R`根据`phylum.counts.txt`计算门水平的A，B组间的显著性差异
- 使用`DESeq2-utils annotation` 根据`DESeq2/DESeq2/A_B/A_B_phylum.txt`在行末尾进行注释并导出`DESeq2/DESeq2/A_B/A_B_phylum.annotation.txt`
- 使用`DESeq2_ratio.R`根据`A_B_phylum.annotation.txt`绘制A，B间的散点图并导出`DESeq2/DESeq2/A_B/A_B_phylum.pdf`
- 使用`DESeq2-utils regulation` 进行排序并导出`DESeq2/DESeq2/A_B/A_B_phylum.list`
- 将`DESeq2/DESeq2/A_B/A_B_phylum.pdf`导出`png`格式

`注意事项`：使用脚本实现如A_B,A_C,B_C等不同组间的循环，因bash太臃肿，故不在此展示。

```bash
mkdir -p  /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2.R                      \
    /project/dix-demo/data_analysis/DESeq2/counts/phylum.counts.txt                 \
    --parallel=F A B /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2-utils annotation       \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.txt 0.05           \
    >/project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.annotation.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2_ratio.R                \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.annotation.txt A B \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2-utils regulation       \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.annotation.txt     \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.list ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                       \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_phylum.pdf ;
```

可使用`for`进行循环迭代：


```bash
for i in phylum class order family genus; 
do
mkdir -p  /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2.R                    \
    /project/dix-demo/data_analysis/DESeq2/counts/${i}.counts.txt                 \
    --parallel=F A B /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2-utils annotation     \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.txt 0.05           \
    >/project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.annotation.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2_ratio.R              \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.annotation.txt A B \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/DESeq2-utils regulation     \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.annotation.txt     \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.list ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                     \
    /project/dix-demo/data_analysis/DESeq2/DESeq2/A_B/A_B_${i}.pdf ;
done
```


### 27. tax4fun 功能预测

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/tax4fun/pipeline
mkdir -p /project/dix-demo/data_analysis/tax4fun/report
```

**序列比对 :**

使用`blastn`进行核酸比对

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/blastn \
    -task megablast \
    -db  /biostack/database/SILVA_123/SILVA_123_SSURef                          \
    -query /project/dix-demo/data_analysis/classify/zotus/zotus.fasta           \
    -num_threads 72     \
    -max_target_seqs 1  \
    -evalue 0.001       \
    -outfmt 6           \
    -perc_identity 90   \
    -out /project/dix-demo/data_analysis/tax4fun/pipeline/blastn.txt
```

`blastn`参数解析

```
-db                 指定blast搜索用的数据库
-query              输入序列（输入）
-num_threads        线程数
-max_target_seqs    设置最多的目标序列匹配数
-evalue             设置e值
-outfmt             可以自定义要输出哪些内容
-out                输出结果文件（输出）
```

**ZOTU注释 :**

提取`blastn.txt` 的`1,2`列传递给`tsv-utils` [annotation](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/annotation/annotation.md) 对第二列进行注释，提取结果的`1,3`列，传递给`atlas-utils  annotation`进行物种注释，使用`tsv-utils` [add_headline](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/add_headline/add_headline.md)为输出结果添加标题行`# Constructed from biom file`

```bash
cut -f1,2 /project/dix-demo/data_analysis/tax4fun/pipeline/blastn.txt |                   \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils  annotation -c 2   \
/biostack/database/SILVA_123/SILVA_123_SSURef.txt  - |                                    \
cut -f1,3 |                                                                               \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  annotation -        \
    /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt | \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils  add_headline          \
    "#         Constructed from biom file" -                                              \
    > /project/dix-demo/data_analysis/tax4fun/pipeline/zotu_table.txt
```

**Tax4fun预测 :**


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/Tax4fun.R      \
    /biostack/database/SILVA_123/Tax4fun                             \
    /project/dix-demo/data_analysis/tax4fun/pipeline/zotu_table.txt  \
    /project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt
```


**格式化输出 :**

```bash
perl -ne 's/^(\S+);.+?\t/$1\t/; print' /project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt | \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils view -c  -                    \
>/project/dix-demo/data_analysis/tax4fun/pipeline/tax4fun.txt
```

**格式化元数据 :**

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt |  \
grep -v '#' \
    >/project/dix-demo/data_analysis/tax4fun/report/metadata.txt
```

**分类等级提升 :**


- 提取`mapping_file.txt`的`1,2`列，使用`grep`提取不含`#`的行;
- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转化为`module` 表，使用`tsv-utils definition`根据`db/kegg/module-definition.txt`对结果进行注释;
- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转化为`pathway`表，使用`tsv-utils definition`根据`db/kegg/pathway-definition.txt`对结果进行注释;
- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转化为`catalog`表，使用`tsv-utils definition`根据`db/kegg/catalog-definition.txt`对结果进行注释;

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  kann   \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/ko-module.txt    \
    /project/dix-demo/data_analysis/tax4fun/pipeline/tax4fun.txt |           \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/module-definition.txt - \
    >/project/dix-demo/data_analysis/tax4fun/pipeline/module.txt
```

可使用`for`进行循环迭代：

```bash
for i in module pathway catalog; 
do
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  kann   \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/ko-${i}.txt      \
    /project/dix-demo/data_analysis/tax4fun/pipeline/tax4fun.txt |           \
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/${i}-definition.txt -    \
    >/project/dix-demo/data_analysis/tax4fun/pipeline/${i}.txt
done
```

**降唯可视化 :**


- 创建目录
- 分别对`ko`表、`module`表、`pathway`表和`catalog`表进行主成分分析、主坐标分析、非度量多维尺度。

```bash
mkdir -p /project/dix-demo/data_analysis/tax4fun/prediction/ko ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T            \
    /project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt                      \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T  -m bray  \
    /project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt                      \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray   \
    /project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt                      \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/ko/ko.nmds.pdf ;
```

可使用`bash`进行循环迭代：


```bash
for i in ko module pathway catalog; 
do
mkdir -p /project/dix-demo/data_analysis/tax4fun/prediction/${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T            \
    /project/dix-demo/data_analysis/tax4fun/pipeline/${i}.txt                    \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T  -m bray  \
    /project/dix-demo/data_analysis/tax4fun/pipeline/${i}.txt                    \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray   \
    /project/dix-demo/data_analysis/tax4fun/pipeline/${i}.txt                    \
    /project/dix-demo/data_analysis/tax4fun/report/metadata.txt                  \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                    \
    /project/dix-demo/data_analysis/tax4fun/prediction/${i}/${i}.nmds.pdf ;
done
```

**合并信息 :** 


合并`ko.txt`、`module.txt`、`pathway.txt`、`catalog.txt`为`tax4fun.kegg.xlsx`excel

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils tsv2xlsx  \
    /project/dix-demo/data_analysis/tax4fun/report/tax4fun.kegg.xlsx         \
    ko:/project/dix-demo/data_analysis/tax4fun/pipeline/ko.txt               \
    module:/project/dix-demo/data_analysis/tax4fun/pipeline/module.txt       \
    pathway:/project/dix-demo/data_analysis/tax4fun/pipeline/pathway.txt     \
    catalog:/project/dix-demo/data_analysis/tax4fun/pipeline/catalog.txt
```


### 28. picrust2 功能预测


`注意事项:`  删除目录（因目录存在会导致picrust2_pipeline.py报错，故先删除）

```bash
rm -rf /project/dix-demo/data_analysis/picrust2/pipeline
mkdir -p /project/dix-demo/data_analysis/picrust2/report
mkdir -p /project/dix-demo/data_analysis/picrust2/prediction
```

**picrust2预测功能组成 :**

```bash
picrust2_pipeline.py \
      -s  /project/dix-demo/data_analysis/classify/zotus/zotus.fasta         \
      -i  /project/dix-demo/data_analysis/classify/zotus/zotu_table_norm.txt \
      -o  /project/dix-demo/data_analysis/picrust2/pipeline                  \
      -p  72        \
      --stratified  \
      --wide_table  \
      --per_sequence_contrib
```

`picrust2_pipeline.py`参数解析

    -s                      OTU或ASV的序列文件
    -i                      序列丰度表
    -o                      输出目录（输出）
    -p                      线程数
    --stratified            在各层级产生分层的表，即功能对应物种来源
    --per_sequence_contrib  计算每条序列的贡献，即将计算每个个体的通路，只有当—coverage打开时才计算分层的覆盖度


`picrust2_pipeline.py` 命令参考：[https://github.com/picrust/picrust2/wiki/Full-pipeline-script](https://github.com/picrust/picrust2/wiki/Full-pipeline-script)


**格式化元数据 :**

提取`mapping_file.txt`的`1,2`列及不含`#`的列；


```bash
cut -f1,2 /project/dix-demo/mapping_file.txt | \
   grep -v '#'                                 \
   >/project/dix-demo/data_analysis/picrust2/report/metadata.txt
```

**添加注释信息 :**

使用`tsv-utils` [definition](https://github.com/atlashub/biostack-suits-docs/blob/master/tsv-utils/definition/definition.md) 根据`picrust2`的数据库进行注释

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/picrust2/ko_info.tsv \
    /project/dix-demo/data_analysis/picrust2/pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > \
    /project/dix-demo/data_analysis/picrust2/prediction/ko.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/picrust2/ec_level4_info.tsv  \
    /project/dix-demo/data_analysis/picrust2/pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz   \
    >/project/dix-demo/data_analysis/picrust2/prediction/enzyme.txt ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/picrust2/metacyc_pathways_info.txt\
    /project/dix-demo/data_analysis/picrust2/pipeline/pathways_out/path_abun_unstrat.tsv.gz  \
    >/project/dix-demo/data_analysis/picrust2/prediction/pathway.txt ;

```

**降维数据可视化 :**


- 创建目录
- 分别对`ko`表、`enzyme`表、`pathway`表进行主成分分析、主坐标分析、非度量多维尺度

```bash
mkdir -p /project/dix-demo/data_analysis/picrust2/prediction/ko ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T          \
    /project/dix-demo/data_analysis/picrust2/prediction/ko.txt                 \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T -m bray \
    /project/dix-demo/data_analysis/picrust2/prediction/ko.txt                 \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray \
    /project/dix-demo/data_analysis/picrust2/prediction/ko.txt                 \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/ko/ko.nmds.pdf ;
```

可使用`for`进行循环迭代：


```bash
for i in ko enzyme pathway; 
do
mkdir -p /project/dix-demo/data_analysis/picrust2/prediction/${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T          \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}.txt               \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T -m bray \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}.txt               \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}.txt               \
    /project/dix-demo/data_analysis/picrust2/report/metadata.txt               \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                  \
    /project/dix-demo/data_analysis/picrust2/prediction/${i}/${i}.nmds.pdf ;
done
```

**合并文件 :**

合并`ko.txt`、`enzyme.txt`、`pathway.txt`为`picrust2.metagenome.xlsx`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils tsv2xlsx   \
    /project/dix-demo/data_analysis/picrust2/report/picrust2.metagenome.xlsx  \
    ko:/project/dix-demo/data_analysis/picrust2/prediction/ko.txt             \
    enzyme:/project/dix-demo/data_analysis/picrust2/prediction/enzyme.txt     \
    pathway:/project/dix-demo/data_analysis/picrust2/prediction/pathway.txt 
```

### 29. kegg 功能关联

**创建目录 :**

```bash
mkdir -p /project/dix-demo/data_analysis/kegg/report
mkdir -p /project/dix-demo/data_analysis/kegg/annotation
```

**格式化元数据 :**


提取`mapping_file.txt`的`1,2`列及不含`#`的列；

```bash
cut -f1,2 /project/dix-demo/mapping_file.txt | \
  grep -v '#'                                  \
  >/project/dix-demo/data_analysis/kegg/report/metadata.txt
```

**格式化ko注释文件 :**

解压`pred_metagenome_unstrat.tsv.gz` 然后使用`tsv-utils view`去除注释行并在第一行添加`#`

```bash
gunzip -c /project/dix-demo/data_analysis/picrust2/pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz | \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils view -c - \
    >/project/dix-demo/data_analysis/kegg/report/ko.txt
```

**分类等级提升 :**


- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转为`module` 表然后使用`tsv-utils definition`根据`db/kegg/module-definition.txt` 进行注释；
- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转为`pathway`表然后使用`tsv-utils definition`根据`db/kegg/pathway-definition.txt`进行注释；
- 使用`atlas-utils` [kann](https://github.com/atlashub/biostack-suits-docs/blob/master/atlas-utils/kann/kann.md) 将`ko`表转为`catalog`表然后使用`tsv-utils definition`根据`db/kegg/catalog-definition.txt`进行注释；


```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  kann    \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/ko-module.txt     \
    /project/dix-demo/data_analysis/kegg/report/ko.txt |                      \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " " \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/module-definition.txt -      \
    >/project/dix-demo/data_analysis/kegg/annotation/module.txt
```


可使用`for`进行循环迭代：


```bash
for i in module pathway catalog; 
do
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/atlas-utils  kann          \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/ko-${i}.txt                 \
    /project/dix-demo/data_analysis/kegg/report/ko.txt |                                \
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils definition -d " "\
    /biostack/tools/protocols/dix-seq-0.0.2/bins/../db/kegg/${i}-definition.txt -       \
    >/project/dix-demo/data_analysis/kegg/annotation/${i}.txt
done
```

**降维数据可视化 :**

- 创建目录
- 分别对`module`表、`pathway`表、`catalog`表进行主成分分析、主坐标分析、非度量多维尺度

```bash
mkdir -p /project/dix-demo/data_analysis/kegg/annotation/module ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T           \
    /project/dix-demo/data_analysis/kegg/annotation/module.txt                  \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/module/module.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/module/module.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T -m bray  \
    /project/dix-demo/data_analysis/kegg/annotation/module.txt                  \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/module/module.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/module/module.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray  \
    /project/dix-demo/data_analysis/kegg/annotation/module.txt                  \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/module/module.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/module/module
    .nmds.pdf ;
```

可使用`for`进行循环迭代：


```bash
for i in module pathway catalog; 
do
mkdir -p /project/dix-demo/data_analysis/kegg/annotation/${i} ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCA.R -g T -t T           \
    /project/dix-demo/data_analysis/kegg/annotation/${i}.txt                    \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.pca.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.pca.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/PCoA.R -g T -t T -m bray  \
    /project/dix-demo/data_analysis/kegg/annotation/${i}.txt                    \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.pcoa.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.pcoa.pdf ;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/NMDS.R -g T -t T -m bray  \
    /project/dix-demo/data_analysis/kegg/annotation/${i}.txt                    \
    /project/dix-demo/data_analysis/kegg/report/metadata.txt                    \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.nmds.pdf;
/biostack/tools/protocols/dix-seq-0.0.2/bins/../utils/pdf2png                   \
    /project/dix-demo/data_analysis/kegg/annotation/${i}/${i}.nmds.pdf ;
done
```

**合并文件 :**

合并`module.txt`、`pathway.txt`、`catalog.txt`为`picrust2.kegg.xlsx`

```bash
/biostack/tools/protocols/dix-seq-0.0.2/bins/../binaries/tsv-utils tsv2xlsx \
    /project/dix-seq/data_analysis/kegg/report/picrust2.kegg.xlsx           \
    module:/project/dix-seq/data_analysis/kegg/annotation/module.txt        \
    pathway:/project/dix-seq/data_analysis/kegg/annotation/pathway.txt      \
    catalog:/project/dix-seq/data_analysis/kegg/annotation/catalog.txt
```


本文材料为 **BASE (Biostack Applied bioinformatic SEies ) 课程 Linux Command Line Tools for Life Scientists** 材料， 版权归 **上海逻捷信息科技有限公司** 所有。