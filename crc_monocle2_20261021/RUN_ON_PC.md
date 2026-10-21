# PC环境与运行顺序

本文件给能访问本地终端和文件的GPT/Codex执行。普通网页聊天无法自行使用PC资源。哥哥的PC系统及软件路径尚未核实，不套用Mac或其他PC路径。

## 1. 获取输入

从 `https://github.com/COSMOS340/exchange` 克隆仓库，进入 `crc_monocle2_20261021`。安装GitHub CLI的PC可执行：

```bash
gh release download crc-monocle2-input-v1-20261021 --repo COSMOS340/exchange --pattern input_data.zip --pattern INPUT_ARCHIVE.json --dir .
```

也可以从Release网页下载这两个附件。读取`INPUT_ARCHIVE.json`中的SHA256/字节数，使用系统校验工具核对zip后再解压到本目录。不要解压到仓库根，目录应为`crc_monocle2_20261021/data/Overall`。随后运行入口的`verify-input`逐文件检查。`INPUT_MANIFEST.json`在代码目录，附件内部包含同一版本。

## 2. 环境检查

先列出本机实际解释器和版本。使用隔离环境，Python 3.11及以上，依赖见requirements.txt；记录最终`pip freeze`。本地准备工具已经在Python3.13及实际安装依赖上做软件测试，不表示PC环境已验证。

R使用4.5.x；Monocle2包名是`monocle`，版本2.38.0；DDRTree0.1.6、igraph2.2.2。采用对应Bioconductor版本/官方包归档，不随意装“最新”后取消兼容检查。还需要Matrix、Rcpp、data.table、jsonlite、rhdf5、dplyr、rlang、R.utils；Slingshot阶段额外需要slingshot及其依赖。历史R sessionInfo在reference下作为参照，不是PC已安装证明。

Windows的Rcpp需要与R版本兼容的Rtools；Linux需要C++编译工具。用实际Rscript测试Rcpp小函数。不要自动改全局R库、安装路径或系统配置。R与Python可以是Windows原生，也可以是已经存在的Linux环境，但同一次运行的路径格式必须统一，不猜盘符。

查阅官方包信息：[Monocle](https://bioconductor.org/packages/monocle/)、[Slingshot](https://bioconductor.org/packages/slingshot/)、[Rtools](https://cran.r-project.org/bin/windows/Rtools/)。先验证版本和函数，再运行附带图算法的合成等价性检查。

使用核实的Rscript执行`tests/check_graph_compatibility.R reference/Overall_Monocle2_generated.R`，保存输出JSON和stderr。该测试只用合成坐标/图验证MST总权重、零距离策略、DFS和补丁安装，不运行生物学拟合；实际模型还须核对根距离和树结构。

## 3. 运行

以下`python`表示已经核实的Python解释器。`--rscript`必须填写实际发现的Rscript完整路径；示例中不提供虚构机器路径。

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python Monocle2_PC.py verify-input
python Monocle2_PC.py prepare
```

检查`work/preparation_summary.json`及`WORK_MANIFEST.json`后，运行：

```text
python Monocle2_PC.py fit --rscript "实际核实的Rscript完整路径" --method Monocle2
```

Monocle2对象保存并检查后：

```text
python Monocle2_PC.py fit --rscript "同一Rscript完整路径" --method Slingshot
```

每条fit都在前台等待子进程结束，日志写在`results/04_日志和检查`。接手GPT可以用该机器可靠的后台会话管理运行，但必须记录PID、命令、目录和退出码，不能退出会话后误称仍运行。不要同时启动两份。

`prepare`已存在work时拒绝覆盖；失败先查日志和输入，不删除断点盲重跑。有效DDRTree断点可以在代码/输入/环境绑定一致时恢复；DDRTree内部未保存时，不能声称可以从任意迭代断点接续。

## 4. 后处理

依照TASK_SPEC实现到`results/code/`，由同一工作目录读取已保存表格/对象。参考文件仅用于移植，不能原样执行其全局变量，也不要导入原Mac大型主脚本。

本交接不提供“已跑通的全部后处理一键按钮”。接手GPT需完成图、动态、独立验证、逐图阅读指南和回传。已有表格即可绘图时不重新拟合。正式全量画图可以栅格化，不能随机减少点数。

64GB是用户提供的物理内存，不是成功保证。异常必须保留日志与现场信息；不能为了继续运行改为抽样或分块树。图排序补丁避免显式N×N存储，计算仍覆盖全体细胞。
