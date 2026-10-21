# CRC上皮Monocle2 PC输入包 v1

任务代码与说明：[crc_monocle2_20261021](https://github.com/COSMOS340/exchange/tree/main/crc_monocle2_20261021)

这是输入交接，不是计算结果。主拟合范围150,917个CRC来源上皮，正常/癌旁来源保留在223,736细胞的归档输入内。先读取START_HERE_FOR_GPT、AGENTS和TASK_SPEC。哥哥PC运行全部主细胞，不抽样、不分块独立建树。

附件：

- `input_data.zip`：计数H5、同序元数据/基因表、PCA/Harmony、Hallmark资源与逐文件清单。
- `INPUT_ARCHIVE.json`：整包SHA256及字节数，先核验再解压。

解压到仓库`crc_monocle2_20261021`内。遵循RETURN_PROTOCOL将输出放入新结果Release，不覆盖这个输入版本。许可和出处见SOURCES_AND_LICENSES。

本地8项Python软件测试、6项合成图权重/DFS检查已通过；未进行正式PC环境确认、CRC全量拟合或结果验收。后处理需要接手GPT按任务要求实施，参考文件不应原样执行。
