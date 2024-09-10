# 1. github项目删除部分文件

```bash
cd path/to/your/github/project  # 切换到你的GitHub项目目录
git rm file1.txt file2.txt     # 删除file1.txt和file2.txt文件
git commit -m "Remove some files"  # 提交更改
git push origin main          # 推送到GitHub（注意：根据你的分支名称可能需要修改，如master或其他分支名）

cd MicrobiomeStatPlot(master)
git rm -r 22ManhattanPlot
git commit -m "Remove specific folders"
git push origin master
#git pull
```

