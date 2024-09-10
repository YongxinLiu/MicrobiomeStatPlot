# 学习前准备，请将data目录复制到C:根目录，或服务器家目录(~/)

# 熟悉工作环境

# 进入C盘
cd /c

# 创建一个目录
mkdir -p amplicon

# 查看文件列表
ls

# 显示当前工作目录 print working directory
pwd 

# cd 切换工作目录
cd /c/amplicon/12Linux
pwd



# 编写一个shell程序 concatenate files and print on the standard output
cat > test.sh

# 输入如下内容，不包括开头的#和空格，按Ctrl+D结果编辑并保存
# #! /bin/bash
# echo 'Hello YSX!'

# 让文件变为程序 change file mode bits
chmod +x test.sh

# 运行程序
./test.sh



# 常用命令

# 显示当前文件夹文件 list directory contents
ls # ls是list的缩写
ls -l # 列表显示

# 新建文件夹 make directories 
mkdir test # 创建test目录

# 拷贝文件，原文件至目标位置
cp test.sh test/test.txt

# 进入文件夹
cd test

ls # 查看文件夹内容

# 切换至上级目录
cd .. 

# 拷贝文件，原文件至目标位置：cp是copy的缩写
cp test.sh file_temp.txt # 复制文件
cp test.sh test/ # 复制文件到指定目录

# 移动或改名文件:mv是move的缩写
mv test.sh temp.sh # 移动，不更新目录为改名

# 拷贝文件，原文件至目标位置
rm test/test.sh # 文件
rm -r test # 删除文件夹


# 快捷键

# Tab键补全
# cat f # 单选时自动补全
# ls s # 多选时提示侯选

# 中止命令 Ctrl+C
ping -t www.baidu.com



# fastq文件操作

# 按列表显示文件详细
ls -l example.fq.gz

# 解压缩
gunzip example.fq.gz

ls -l example.fq

# 显示文件前10行
head example.fq

# 按页查看文件,-S不换行，空格翻页，q退出
less -S example.fq

# 转换fastq为fasta，并按顺序重复命名序列
awk 'NR%4==2 {print ">E"NR/4+0.5"\n"$0}' example.fq > example.fa

# 显示fasta文件末尾10行
tail example.fa

# 查找某条序列
grep 'AAAACACAGGAACCTGGGTGAAAAC' example.fa | head 

# 统计序列条数
grep -c '>' example.fa

# 统计序列长度
grep -v '>' example.fa | awk '{print length($0)}'| head

# 统计序列长度分布
grep -v '>' example.fa | awk '{print length($0)}'|sort|uniq -c

# 清理临时文件
rm example.fa file_temp.txt temp.sh
# 压缩数据节约空间
gzip example.fq
