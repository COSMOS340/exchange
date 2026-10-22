完成全部 **223,736** 个上皮细胞的共同表达轨迹分析（正常来源72,819、CRC150,917）；原CRC-only结果保留。本次是横断面状态排序，不能解释为真实时间、恶性克隆演化或Normal→Adjacent→CRC癌变路线。

- Monocle2：全细胞一次DDRTree，主/替代根、7 State及各4条路径，3个完整RDS。
- Slingshot：同细胞30维Harmony，实际根粗簇C2/C4、11/10条路径，4个完整RDS，所有PT/权重/NA及原生UMAP投影保留。
- 动态基因、2个表达模块及GO全表，六套Hallmark，来源/CNV、旧CRC-only、全部根/路径及8队列敏感性均包含源表。
- PC数值核验通过；全部97张正式PNG及97个PDF的97页均实际打开检查，记录绑定最终文件SHA。
- 原生GLM/离散度/Quick-TRANSfer警告及Slingshot停止条件的不确定性如实保留；剔除GSE178341的MST因矩阵奇异不可评估，未将其填为成功。

结果分支：`pc-all-epithelial-results-20261022`。本Release精确代码提交：`c25fa02d6b89a33609e548b7db04069d1fc523b8`。实际Monocle2拟合提交为`bdcc8df4fdf927b9d2db11ad349c1e60da76b8d8`，Slingshot为`adeb6886ff58a7779dc36ac8e961eb6937112691`。

下载全部`results.zip.partNNN`及清单到**新目录**。分片仅用于传输，每片至多1.5GiB；按文件名顺序拼接即可还原一个完整ZIP：

```bash
cat results.zip.part* > results.zip
shasum -a 256 results.zip
unzip results.zip -d all_epithelial_received
```

总ZIP的字节数/SHA与每片校验值见`RETURN_MANIFEST.json`；全部结果文件见`RESULT_MANIFEST.json`，冻结任务/配置侧文件见`SIDECAR_MANIFEST.json`。本机已逐个读取ZIP成员核验CRC、字节和SHA，并核对分片拼接SHA，详见`LOCAL_TRANSPORT_VERIFICATION.json`。无需重复下载/上传原始1.18GB输入。

先读附件`RUN_SUMMARY.md`；解压后阅读`results/01_逐图阅读指南.md`及`results/02_方法与限制.md`。完整对象、表格、图、代码、版本、运行及失败重试日志均在包内。

Mac端仍需下载校验与接收复核，通过前不自动替换第14步正式入口。[旧CRC-only结果](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-pc-results-v3-20261022)保持不变。
