# PC环境与运行顺序

## 获取输入

在PC克隆 `https://github.com/COSMOS340/exchange`，进入本目录。下载既有输入Release：

```bash
gh release download crc-monocle2-input-v1-20261021 \
  --repo COSMOS340/exchange \
  --pattern input_data.zip \
  --pattern INPUT_ARCHIVE.json \
  --dir .
```

先按 `INPUT_ARCHIVE.json` 核对字节数和SHA256，再解压到本任务目录，得到 `data/Overall` 与 `data/gene_sets`。不要复用旧任务的 `work/`、`results/` 或CRC-only准备结果。

## 环境

- 先读取机器实际R/Python路径、版本和可用内存，不猜盘符。
- 优先复用旧任务已经验证的隔离环境：Monocle2 2.38.0、DDRTree 0.1.6、igraph 2.2.2及相容R版本。
- 如果复用环境，仍需在本任务日志中保存R `sessionInfo()`、Python依赖和图算法合成测试结果。
- 不在GitHub Actions运行，不同时启动两份重计算，不因内存压力改为随机抽样或分块建树。

## 实施顺序

1. 建立本任务独立代码和结果目录。
2. 从旧任务移植输入校验和Monocle2兼容实现，只修改明确依赖CRC-only范围的代码。
3. 写并运行单元测试：来源三分类、223,736细胞保留、矩阵顺序、根排序、树结构、分箱、等权汇总和动态模块源表。
4. 运行全量输入检查与准备，保存准备摘要和清单。
5. 前台或可审计后台会话串行运行Monocle2；记录PID、命令、目录、开始/结束时间、退出码和峰值内存。
6. 主对象保存且数值检查通过后运行Slingshot。
7. 生成两套组织图、CNV/身份覆盖、动态基因热图、功能活动、敏感性和逐图阅读指南。
8. 完成数值与视觉验收，按RETURN_PROTOCOL打包和上传。

旧任务的 `Monocle2_PC.py` 不能原样运行，因为其准备阶段固定选择 `tissue_group == "Primary_CRC"`。如从旧代码派生，必须在本目录保存修改后的独立入口、差异说明和测试，且输出路径全部指向本任务。

