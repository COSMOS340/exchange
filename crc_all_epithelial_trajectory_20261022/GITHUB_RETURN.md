# 全上皮轨迹：GitHub交付

[完整结果Release](https://github.com/COSMOS340/exchange/releases/tag/crc-all-epithelial-trajectory-results-v1-20261022)已公开发布。分析代码/文档提交为`c25fa02d6b89a33609e548b7db04069d1fc523b8`，结果分支为`pc-all-epithelial-results-20261022`；Release标签精确指向该代码提交。本文件与回传收据是之后的元数据提交，未改写冻结分析包或拟合代码。

## 已完成

全部223,736上皮细胞一次共同分析，7个完整RDS模型、97组正式PNG/PDF（97页）、1,144个结果文件和20个绑定侧文件。PC独立数值核验、逐图逐页实际视觉检查及文件SHA均通过。完整ZIP逐成员CRC/字节/SHA、5分片拼接SHA也已通过。

GitHub共10个附件（5分片、4份清单/本地传输收据、RUN_SUMMARY.md），每个远端名称、字节数及GitHub服务端SHA256均与本地一致，发布后再次核对。精确标签提交也已核对。证据见[发布收据](PC_DELIVERY_RECEIPTS/GITHUB_PUBLICATION_RECEIPT.json)、[传输收据](PC_DELIVERY_RECEIPTS/LOCAL_TRANSPORT_VERIFICATION.json)和[返回清单](PC_DELIVERY_RECEIPTS/RETURN_MANIFEST.json)。

完整ZIP：6,835,510,263字节，SHA256：

```text
5bb3dce1e66e347cd6448f19e8c6bc0456d9e35773c5b82c388b7a73a84f0ffb
```

在新的接收目录下载Release中全部`results.zip.part001`至`results.zip.part005`，按顺序拼接：

```bash
cat results.zip.part* > results.zip
shasum -a 256 results.zip
unzip results.zip -d all_epithelial_received
```

分片及所有包内文件的SHA见清单。先读Release附件RUN_SUMMARY.md；解压后打开`results/01_逐图阅读指南.md`和`results/02_方法与限制.md`。归档内摘要保留封存时的发布待办，后续上传/传输状态以本交付文件和发布收据为准，不改动原包。

## 失败、限制与未做

初次上传在草稿创建成功后因按标签读取草稿返回404退出；已核对原草稿作者、ID、提交、正文与空附件状态后恢复同一草稿。之后一次本地代理连接重置，已保留真实失败/重试日志。未创建替代Release、覆盖旧附件或把失败记录改成成功。

科学限制见完整方法：原生GLM/参数离散度/Quick-TRANSfer警告保留；Slingshot无迭代历史，不能证明以阈值停止；剔除GSE178341的MST因矩阵奇异不可评估。PC数值一致性不等于所有原生优化器收敛，也不证明临床进程或因果。

Mac下载、来源/根/统计/图形接收复核及第14步正式入口更新仍待Mac端执行。旧CRC-only模型、结果分支和[Release](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-pc-results-v3-20261022)保留。
