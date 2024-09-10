## 百度云链接

### R语言常用包合集

Windows 10版4.0.x：链接：链接：https://pan.baidu.com/s/1sSQEasVGH2lUkmI0mLq62w 提取码：vfa1 

备用链接：http://210.75.224.110/db/R/4.0.zip

下载后解压至当前目录，可获得4.0目录，剪切并替换文档 - R - win-library - 4.0 目录即可

Mac OS版4.0.x：链接：链接：https://pan.baidu.com/s/1fh7I-HcqUoGYiu_3XVLIRg 提取码：0315

下载后解压至当前目录，可获得library目录，复制到 /Library/Frameworks/R.framework/Versions/4.0/Resources/library/


## 扩增子数据库

Windows系统下搭建扩增子分析平台所需软件合辑

链接：https://pan.baidu.com/s/1hdu_umrScnV3DpwZtnwIDQ 提取码：yjk9 

## 宏基因组数据库

### KneadData质量控制和去宿主

- 人类基因组Bowtie 2索引 Homo_sapiens_Bowtie2_v0.1.tar.gz (3.44GB) 链接：https://pan.baidu.com/s/1cX02cyoPFJSbh5Hxehwg-g 
提取码：l967

### HUMAnN2物种和功能注释

- metaphlan2数据库 链接：https://pan.baidu.com/s/1hIt-BQOg6J9AKyK14yHKVg 提取码：0315 
- 微生物泛基因组数据库 full_chocophlan_plus_viral.v0.1.1.tar.gz (5.37GB) 链接：https://pan.baidu.com/s/1ANH6s3ULjVemomofv7RwZw 提取码：i8qv
- 蛋白功能注释数据库 uniref90_annotated_1_1.tar.gz (5.87GB) 链接：https://pan.baidu.com/s/1dNKnVjBqPwYBYSC75AWJtQ 提取码：ospn

### HUMAnN3物种和功能注释

- 微生物泛基因组数据库 full_chocophlan.v201901.tar.gz (14.59GB) 链接：https://pan.baidu.com/s/1tQwWG79nGkKJD68Q4moGCQ 提取码：5opp
- 微生物泛基因组数据库 full_chocophlan.v296_201901.tar.gz 2020/7/8更新 (15.3GB) 链接：链接：https://pan.baidu.com/s/17J4fBPdUZhByRhbR5CZ1GQ 提取码：7ij2
- 蛋白功能注释数据库 uniref90_annotated_v201901.tar.gz (19.31GB) 链接：https://pan.baidu.com/s/1b0FwFXEGuLUjhkEJbMbunA 
提取码：0x24
- 蛋白功能注释数据库 uniref50_annotated_v201901.tar.gz (6.49GB) 链接：https://pan.baidu.com/s/1rHxFM0sNmzcEfQylRzfqqA  提取码：wo4r
- full_mapping_v201901.tar.gz 链接：https://pan.baidu.com/s/1YXOQw27qMMNKBLGI1wQ-1A 提取码：13gv

#### MetaPhlAn 3.0 物种注释

- metaphlan_databases (382.14MB) 链接：https://pan.baidu.com/s/1wCkALsoAzTsdJwUQbAcQIA 提取码：wdv9

上传文件夹内容至当前环境库中python3.7包中metaphlan的数据库文件夹中，如~/miniconda2/envs/humann3/lib/python3.7/site-packages/metaphlan/metaphlan_databases

### Kraken 2 数据库

- 标准数据库 kraken2.tar.gz.0/1/2 (34.59GB) 链接：https://pan.baidu.com/s/16pXppGAizMKdU3z5Ji3dzA 
提取码：0315

分卷压缩，下载后需要cat合并再解压 

	cat kraken2.tar.gz.* > merge.tar.gz
	tar xvzf  merge.tar.gz


### EggNOG 5.0数据库

2020年7月8日下载

- EggNOG数据库 eggnog.db.gz (7.87GB) 链接：https://pan.baidu.com/s/1K8LNeHGW-voAGkozHWVrEQ 提取码：rcds
- 蛋白序列的diamond索引 eggnog_proteins.dmnd.gz (4.85GB) 链接：https://pan.baidu.com/s/1WB0BOW08wtoPEYybHLAUUw 提取码：d2uk

### Metawrap 1.3 依赖数据库

#### NCBI_nt (2018/11/16)

NCBI_nt核酸数据库的BLAST索引，此版本可以支持metawrap默认的Blast

- NCBI_nt (41GB, 60个文件) 链接：https://pan.baidu.com/s/1tFpjzEORnV2qoCMY9pAcsA 提取码：0315 

- NCBI_tax (44M，1个文件) 链接：https://pan.baidu.com/s/14ROMoKUX9Tsd_gVNNe_Stg 提取码：0315 


    # 备用nt下载链接
	mkdir -p NCBI_nt && cd NCBI_nt
    wget -c http://210.75.224.110/db/metawrap/NCBI_nt_181116/filelist.txt
    for a in `cat filelist.txt`; do \
	  wget -c http://210.75.224.110/db/metawrap/NCBI_nt_181116/$a; done
    for a in nt.*.tar.gz; do tar xzf $a; done
    # 配套tax下载链接
	mkdir -p NCBI_tax && cd NCBI_tax
    wget -c http://210.75.224.110/db/metawrap/NCBI_tax_181116/taxdump.tar.gz
	tar xzf taxdump.tar.gz


## 备份公网链接

https://github.com/YongxinLiu/MicrobiomeStatPlot

同步备份站点 http://210.75.224.110/github/MicrobiomeStatPlot