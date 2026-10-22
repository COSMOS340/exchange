# 任务需求与验收

## 科学问题

在全部正常来源与原发CRC来源上皮细胞的共同表达轨迹中，两类组织来源如何分布？进一步沿用项目既有来源审计后，Normal、Adjacent和CRC是否占据不同的状态区段？相对RNA CNV、宽上皮身份及表达程序的变化在多大程度上反映组织来源或细胞组成？

本分析描述横断面上皮状态结构。它不证明正常组织随真实时间发展为CRC，也不证明CNV沿拟时序逐步积累。

## 冻结输入

复用Release `crc-monocle2-input-v1-20261021`。`data/Overall/counts.h5` 为15,635基因 × 223,736细胞的CSC原生整数计数，H5键为 `data/indices/indptr/shape`。基因行对应 `genes.csv`，细胞列对应 `metadata.csv.gz`；`Harmony.npy`、`PCA.npy` 为223,736 × 30并与元数据同序。不得猜矩阵方向。

所有223,736个细胞均满足冻结的 `trajectory_eligible=True` 与 `cnv_summary_eligible=True`。本任务不再筛选 `tissue_group == "Primary_CRC"`，不得运行旧任务未修改的CRC-only准备步骤。

关键字段严格使用原名：`cell_id`、`dataset`、`GSM`、`sample_label`、`sample_id`、`capture_id`、`subject_key`、`tissue_group`、`broad_identity`、`coarse_state`、`CNV_relative_band`、`cnv_abs_log2`、`nCount_input`、`Phase`。`orig.ident`在冻结输入中为常量，不能当作建库编号。

## A. 输入与来源核查

- 验证Release附件SHA256、全部文件清单、矩阵形状、非负整数、唯一cell_id/基因名、列UMI和元数据顺序。
- 核对总细胞223,736、正常来源72,819、CRC来源150,917、捕获单元242。
- 依据 `SOURCE_TISSUE_AUDIT.md` 和 `prior_specimen_source_audit.csv` 复用旧 `normal_origin` 值 `Normal_unspecified_distance`、`Explicit_adjacent`、`Not_normal`，同时核对原始 `tissue_group`，要求0个未匹配细胞；保存逐细胞映射表和按队列/捕获单元/标本的计数。
- 两分类直接使用原始 `tissue_group`。三分类图显示为Normal、Adjacent、CRC；Normal不得改写为健康供者或远端正常。
- 保存输入、来源规则和代码哈希。任何计数不一致都停止，不自动放宽。

## B. 全量共同Monocle2轨迹

- 全部223,736个细胞一次共同拟合一棵Monocle2 DDRTree；不抽样、不分块独立建树、不按组织或CNV拆树。
- 沿用输入中的冻结排序基因清单，并重新报告其在全量分析池的检出覆盖。不得使用组织差异基因、CNV差异基因或动态基因反向选择排序基因。
- 沿用已验证的稀疏矩阵、离散度、DDRTree和图排序兼容实现。保存完整参数、目标函数、收敛信息、警告和峰值内存。
- 根据共同图的端点、宽上皮身份覆盖和多个标本支持选择计算主根；另选结构上合理的替代根。根选择不使用组织标签、CNV或功能分数。
- 主根可以由normal占比较高的端点支持，但不能以“必须产生normal到CRC顺序”为条件。保存请求根、实际root_state、root cell_id、回退记录及根覆盖表。
- 主根和替代根共享同一DDRTree坐标与树，只重新排序，不重复拟合制造不同树。
- 保存原始CDS、主/替代根对象、逐细胞坐标/Pseudotime/State、中心坐标、边表、路径表和目标函数。

## C. 组织来源与CNV图

所有点图使用相同全量坐标、圆点、栅格化和固定比例，不随机抽细胞。每张图同时输出PNG和PDF及源表。

1. 共同树按两分类 `Normal_colorectal` / `Primary_CRC` 着色。
2. 两分类分面，共用完整树背景和坐标。
3. 共同树按三分类来源着色。
4. 三分类分面，共用完整树背景和坐标。
5. 同一树按Pseudotime、State、CNV低中高、连续CNV、宽上皮身份、队列、细胞周期和UMI着色。
6. 各组织来源沿拟时序的细胞密度和组成；同时给出合并细胞描述与捕获单元/标本等权描述，分母写清。
7. 逐队列组织分布及留一队列敏感性。三分类中的Adjacent仅来自GSE245552，必须单独展示队列构成，不能将Adjacent差异与队列差异分开声称。
8. 在共享宽上皮身份中重复关键组织比较，避免吸收细胞、干/祖细胞等构成差异被误写为连续状态转变。
9. CRC来源子集叠加旧CRC-only结果，用cell_id精确对齐比较Pseudotime、State/路径和CNV趋势；新分析不覆盖旧对象。
10. 主根与替代根对照；实际根相同或方向翻转时按实报告。

## D. 拟时序动态基因与功能模块

严格执行 `REFERENCE_FIGURE_REQUIREMENTS.md`：

- 生成基因 × 拟时序箱的Z-score热图，按数据驱动的动态模式分模块。
- 标注代表基因、模块GO Biological Process条目、组织来源/CNV/宽身份/队列组成。
- 动态基因检验以捕获单元/标本等权的分箱表达为核心，不将细胞作为独立生物学重复。
- 保存每个基因的覆盖、效应幅度、统计结果、模块、平滑表达及每个GO条目的完整命中信息。
- 动态模块必须进行主根/替代根、Monocle2/Slingshot和逐队列稳定性检查。

## E. 既有功能活动

沿用旧任务冻结的六套Hallmark：OXIDATIVE_PHOSPHORYLATION、INTERFERON_GAMMA_RESPONSE、INTERFERON_ALPHA_RESPONSE、TNFA_SIGNALING_VIA_NFKB、E2F_TARGETS、MYC_TARGETS_V2。

- 使用 `nCount_input` 计算log1p(CP10K)，在全量223,736个细胞分析池重新计算基因均值/总体标准差后取基因Z值均值。
- 至少10个可变命中基因且覆盖原集合至少50%，不足则记为不可评估。
- 每条主根到叶路径独立缩放0至1并固定分箱。先捕获单元内汇总，再标本/有依据的报告患者标签等权。
- 展示六套活动、连续CNV、CNV三组比例、细胞周期和UMI对照。空箱不补0，覆盖不足不连线。
- 这是表达活动，不是代谢通量、蛋白活性或新GSEA。

## F. Slingshot同细胞对照

- 使用同223,736个细胞、同Harmony表示及同一来源审计字段。
- 起点依据与Monocle2一致且不使用组织/CNV/功能值；保存主/替代起点全部曲线、权重、未分配和NA。
- 曲线投影到原状态UMAP/共同坐标时验证拟时序值不被改变。
- 方法间按cell_id一对一比较共同覆盖、方向、路径对应和动态模块，不挑最高相关路径后称为独立验证。

## G. 统计单位与解释

- 细胞级图用于显示结构，不对细胞直接做独立样本显著性检验。
- 主要趋势和组间描述以 `capture_id` 或标本等权；报告患者标签等权结果只纳入有明确链接支持的记录，并标为报告标签汇总。
- 展示每个箱的细胞数、捕获单元数、标本数、报告患者标签数和队列数。覆盖不足留空。
- 默认不加星号。若执行推断检验，必须预先写明统计单位、模型、协变量、多重校正和效应量，并与描述图分开。
- 若全上皮轨迹出现CNV随拟时序变化而CRC-only轨迹没有，应解释为组织来源/上皮身份构成相关现象，除非CRC内部等权敏感性也支持同方向。

## H. 验收

- 223,736个cell_id与输入一一对应，无重复、随机丢失或标签改写。
- 三分类0个未匹配；两分类和三分类计数与 `scope_summary.json` 一致。
- 坐标和拟时序有限；树连通无环；独立复算根到细胞的加权距离并与拟时序比较。
- 全量对象、表格、动态热图源数据、活动矩阵、树边、路径、环境和日志均可复算。
- 逐图打开PNG、逐页渲染PDF，记录字体、遮挡、裁切、空白、图例和色标检查。
- `results/RUN_SUMMARY.md`逐项写完成/失败/未做，列出实际根、State、路径、细胞/统计单位数量、关键警告和可保留结论。
- `results/01_逐图阅读指南.md`用相对路径嵌图，逐图写问题、轴/颜色、统计单位、实际观察、队列构成和解释边界。

