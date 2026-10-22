# CRC全上皮共同拟时序：PC计算任务

本目录是第14步的新任务，不覆盖 `crc_monocle2_20261021` 的CRC内部CNV分层分析。新任务把正常来源和CRC来源的全部223,736个已审核上皮细胞放入同一条Monocle2表达轨迹，用于描述跨组织上皮状态结构。

## 接手顺序

1. 先读 [AGENTS.md](AGENTS.md)、[任务要求](TASK_SPEC.md)、[组织来源核查](SOURCE_TISSUE_AUDIT.md) 和 [图形要求](REFERENCE_FIGURE_REQUIREMENTS.md)。
2. 复用输入Release [`crc-monocle2-input-v1-20261021`](https://github.com/COSMOS340/exchange/releases/tag/crc-monocle2-input-v1-20261021)，校验附件后解压。
3. 新建独立的准备、模型和结果目录。旧任务的准备结果只含CRC细胞，不能直接用于本任务。
4. 先完成全量Monocle2共同建树，再运行同细胞Slingshot对照和动态基因模块分析。
5. 完成数值检查、逐图视觉检查和科学边界核查后，按 [回传协议](RETURN_PROTOCOL.md) 返回。

## 冻结范围

- 15,635个基因，223,736个细胞。
- 全部正常来源上皮：72,819个；原发CRC来源上皮：150,917个。
- 8个GEO队列，229份标本标签，242个建库/捕获单元。
- 所有细胞共同拟合一棵树，不随机抽样，不按组织、CNV、队列或计算块拆树。
- 组织标签、CNV标签和细胞身份仅用于拟合后的展示与描述，不参与排序基因选择、建树或起点选择。

## 两套组织展示

1. **两分类**：`Normal_colorectal`（72,819）对 `Primary_CRC`（150,917）。
2. **三分类**：沿用项目既有来源审计的 `Normal`（64,654）/ `Adjacent`（8,165）/ `CRC`（150,917）。

三分类沿用既有来源审计的 `normal_origin`：`Normal_unspecified_distance`、`Explicit_adjacent`、`Not_normal`，图中显示为Normal、Adjacent、CRC；CRC同时必须满足 `tissue_group == "Primary_CRC"`。Normal表示取材距离未明确的正常结直肠组织，不等同于健康供者。各组来源、统计单位和队列构成见 [SOURCE_TISSUE_AUDIT.md](SOURCE_TISSUE_AUDIT.md)。

## 主要新增图

- 同一轨迹上的两分类和三分类组织来源图。
- 每组分面图使用相同坐标、比例和完整树背景。
- 组织来源、CNV、宽上皮身份、队列、细胞周期和UMI的轨迹覆盖图。
- 沿拟时序变化的基因热图：数据驱动的动态模块、代表基因、模块GO条目和组织/CNV组成注释。

轨迹是横断面表达状态的排序，不能直接写成正常组织在患者体内随时间发展为CRC。

