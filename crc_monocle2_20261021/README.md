# CRC上皮Monocle2：PC计算交接

## PC最新交付状态（2026-10-22）

本轮全量计算和本地验收已完成，结果已公开发布：[结果Release](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-pc-results-v3-20261022)、[完成状态与接收说明](PC_RETURN_STATUS.md)、[远端附件核验回执](PUBLICATION_RECEIPT.json)。保留拟合警告，1项队列留出MST诊断不可评估；Mac接收复核仍待执行。

以下是原始交接说明，保留用于追溯；其中“尚未运行”等历史状态不代表当前PC结果状态。结果包内状态文件同样保留打包时快照。

这是用于哥哥64GB PC的计算交接包。输入与代码准备完成不等于拟时序分析完成。Mac端未运行本包的正式准备或拟合。

## 接手顺序

1. 先读 [AGENTS.md](AGENTS.md)、[任务要求](TASK_SPEC.md) 和 [环境与运行](RUN_ON_PC.md)。
2. 从仓库 Release `crc-monocle2-input-v1-20261021` 下载 `input_data.zip` 与 `INPUT_ARCHIVE.json`。附件解压到本目录，得到 `data/Overall`、`data/gene_sets`。
3. 核对 `INPUT_MANIFEST.json`，在PC上执行全量准备、Monocle2；保存主树后再串行运行Slingshot和后处理。
4. 完成数值与PNG/PDF视觉检查，并按 [回传协议](RETURN_PROTOCOL.md) 返回。代码需改动时必须保留修改记录。

## 本次范围

- 来源包：223,736个此前审核合格的上皮细胞，包含正常/癌旁72,819个与CRC150,917个。
- **主拟合：全部150,917个CRC来源上皮共同建树，不抽细胞、不按CNV组拆树、不按计算块拆树。**
- CRC范围：8队列、152份标本、163个建库单元、142个报告患者标签。标签数不能直接当作已核实的独立患者数。
- CNV组：Low 41,686、Intermediate 46,755、High 62,476。沿用上游标签，不重新切三等分。
- 方法：Monocle2主分析；Slingshot同细胞对照；主/替代计算参考起点。
- 标题使用“CRC来源上皮”及“相对RNA CNV”，不将组织来源或CNV高组当作恶性身份金标准。

## 文件状态

`Monocle2_PC.py` 提供校验、全量准备、两方法拟合和回传打包。生成R入口保留已有兼容/图排序修补，去掉Mac的1024列预处理实现。没有自动执行拟合的默认命令。

**后处理需要接手GPT按TASK_SPEC适配实现并运行**：提供了 `reference/postprocessing_reference.py` 作为已用公式和绘图方法参考，它不是可直接运行的PC工作流。`*_completed.json`只表示对应方法保存结束，不代表科学验收或交付。

[可选分析](OPTIONAL_ANALYSES.md) 单列，不因“可能需要”而擅自启动显著性检验、另换细胞范围或其他算法。

## 入口

- [给PC端GPT的完整指令](START_HERE_FOR_GPT.md)
- [本地打包与软件测试记录](HANDOFF_VALIDATION.json)
- [来源与第三方许可](SOURCES_AND_LICENSES.md)
- [输入范围核对](scope_summary.json)

仓库公开。输入只保留公开来源的研究标识和既有分析字段；不要上传账号凭据、患者直接身份信息、整篇版权PDF或机器环境中的秘密。
