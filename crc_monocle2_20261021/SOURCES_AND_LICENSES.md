# 来源、署名与公开范围

## 单细胞来源

输入来自已处理的公开GEO研究，保留来源编号及研究内样本/细胞标识；不包含直接身份资料。本交接不转移原作者权利，也不将第三方数据归为本仓库原创。

- [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)
- [GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735)
- [GSE161277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161277)
- [GSE166555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555)
- [GSE178341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341)
- [GSE200997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)
- [GSE231559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231559)
- [GSE245552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552)

附带CNV、身份和坐标是当前研究的派生结果，未作为最终论文结论验收。此公开上传用于用户指定的研究交接，不代表临床用途。

## 算法和参考实现

Monocle由Cole Trapnell及项目贡献者开发，安装包monocle 2.38.0声明Artistic-2.0。本包R兼容函数包含对其函数的适配；修改为进程内dplyr/igraph兼容、显式DFS、全细胞MST存储优化及根回退，不代表原作者发布版本。原始方法与版权归各自作者。

- [Monocle源码和版权](https://github.com/cole-trapnell-lab/monocle-release)
- [Monocle2论文，PMID28825705](https://pubmed.ncbi.nlm.nih.gov/28825705/)
- [Slingshot论文，PMID29914354](https://pubmed.ncbi.nlm.nih.gov/29914354/)
- 本包`reference/monocle_Artistic-2.0.txt`保留许可文本。

## Hallmark

`Hallmark_members.csv`、`Hallmark_database.csv.gz`来自msigdbr26.1.0分发的MSigDB2026.1，原数据库版本保留在表内。代码仅选择TASK_SPEC列出的六套Hallmark，不重新定义基因集合。署名Broad Institute、Massachusetts Institute of Technology及Regents of the University of California；MSigDB适用CC BY4.0及官网列明的集合特定条款。这里不打包KEGG、BioCarta或整篇文献PDF。

- [MSigDB许可](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp)
- [CC BY4.0](https://creativecommons.org/licenses/by/4.0/)
- [Hallmark集合说明](https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H)

本次传递的是既有导出文件；执行时记录实际命中/未命中与CRC内重新标准化，不改变资源归属。
