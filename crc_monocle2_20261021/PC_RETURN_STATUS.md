# PC全量结果公开回传记录（2026-10-22）

本轮全部150,917个CRC来源上皮的Monocle2和Slingshot主/替代根计算已完成，结果已公开上传。原输入Release和旧结果保留，未修改main，也没有启动GitHub Actions。

- [完整结果Release](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-pc-results-v3-20261022)：`crc-monocle2-pc-results-v3-20261022`。
- [结果代码提交](https://github.com/COSMOS340/exchange/commit/4d8a13008e57dc5a30a6fb5b5fc1a48410151cce)：`4d8a13008e57dc5a30a6fb5b5fc1a48410151cce`。结果标签固定指向此计算代码提交。
- [回传分支](https://github.com/COSMOS340/exchange/tree/pc-results-20261021-wsl)：`pc-results-20261021-wsl`。发布后的说明和核验回执是后续文档提交，不改变结果标签及包内代码绑定。
- [输入Release](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-input-v1-20261021)：输入提交`4bb17d4ffc23912ad1cd82ad87f4760d5e0f4955`。
- [运行摘要](results/RUN_SUMMARY.md)、[逐图阅读指南](results/01_逐图阅读指南.md)、[远端附件核验回执](PUBLICATION_RECEIPT.json)。

## 完成与保留限制

7个真实RDS、全部表格、71张PNG和71份单页PDF、逐图指南、代码与环境绑定、原始日志和失败记录均在结果包中。本地独立数值检查和所有图页实际审阅完成；完整ZIP的735个文件已逐文件校验，并核对所有ZIP成员CRC、分片哈希及拼接整包哈希。

保留原生拟合警告，不能用退出0声称全局最优或曲线收敛。8项队列留出MST诊断中7项可计算；留出GSE178341后协方差奇异，该项不可评估，不能称全部敏感性检查通过。全量原生Slingshot基准未跑完，精确线段搜索加速的验证范围与限制详见运行摘要。全部结论为描述性；没有新增动态基因P值，拟时序不证明真实时间或因果。

## 附件与接收

Release共有6个附件：`results.zip.part001`、`results.zip.part002`、`results.zip.part003`、`RETURN_MANIFEST.json`、`RESULT_MANIFEST.json`、`RUN_SUMMARY.md`。已读取GitHub远端列表，逐项核对文件名、上传完成状态、字节数和服务端SHA256，全部一致，核验时刻为`2026-10-22T00:43:09.890577+08:00`。

完整ZIP共4,175,582,141字节，SHA256为：

```text
f8c5a6e36f8a875bbb414a947a532042fa20c0b43408fb820b85b0c7977b5a87
```

逐图指南应与解压后的图、表目录一起阅读；大图和RDS保存在Release附件中，代码仓库不重复存放。

分片仅用于传输，不拆分细胞模型。Mac请在新接收目录下载全部附件，按RETURN_MANIFEST顺序二进制拼接，核对整包SHA256后解压，并按RESULT_MANIFEST逐文件校验。接收步骤及科学复核要求见[回传协议](RETURN_PROTOCOL.md)。不要把旧OOM诊断包当成本轮结果。

包内`RUN_SUMMARY.md`、`code_version.json`、进度记录等保留打包时的原始快照，因此可能记载“尚未公开上传”；当前公开发布及远端附件核验状态以本记录、Release和回执为准。没有为更新发布状态而改写已冻结结果、重新打包或改变模型。

Mac接收与科学复核尚未执行，第14步正式入口未更新。GitHub不会自动通知Mac端GPT，请告诉接手端“PC已回传”，提供上面的结果Release与代码提交，再依照回传协议下载验收。
