# 组织来源核查与既有三分类规则

## 沿用的项目口径

本任务不重新解释“正常黏膜”或按数据集介绍推断取材距离，直接沿用第9B步已经冻结的标本来源审计：

- `Normal_unspecified_distance`：图中显示为 `Normal`。
- `Explicit_adjacent`：图中显示为 `Adjacent`。
- `Not_normal` 且 `tissue_group == "Primary_CRC"`：图中显示为 `CRC`。

旧来源审计文件为第9B步 `02_specimen_source_audit.csv`。本任务已冻结其必要列为 [prior_specimen_source_audit.csv](prior_specimen_source_audit.csv)，共235条标本记录，SHA256为 `71b626ffe5e42c051ae27d293ea87a0b394814182e353943c70b17a4c93542e8`。本次已按 `dataset + specimen_id + tissue_group` 与全上皮冻结元数据做一对多精确关联，223,736个细胞全部匹配，未匹配为0；审计表中不在本次全上皮输入的记录不进入结果。

## 数量

两分类仍保留原始 `tissue_group`：

- `Normal_colorectal`：72,819个细胞，79个捕获单元。
- `Primary_CRC`：150,917个细胞，163个捕获单元。

三分类为：

- `Normal` / `Normal_unspecified_distance`：64,654个细胞，75个捕获单元，72个报告受试者标签，来自7个队列。
- `Adjacent` / `Explicit_adjacent`：8,165个细胞，4个捕获单元，4个报告受试者标签，仅来自GSE245552。
- `CRC` / `Primary_CRC`：150,917个细胞，163个捕获单元，142个报告受试者标签，来自8个队列。

## 精确赋值规则

PC端优先使用随任务冻结的旧来源审计记录与元数据按 `dataset + specimen_id + tissue_group` 精确关联。验收要求 `many_to_one` 且0个未匹配细胞。禁止从cell_id、样本名后缀或字母N/T自行推断。

在核对既有审计后，所得规则应等价于：

```text
if tissue_group == "Primary_CRC" and normal_origin == "Not_normal":
    source_origin = "Not_normal"
    display_3class = "CRC"
else if normal_origin == "Explicit_adjacent":
    source_origin = "Explicit_adjacent"
    display_3class = "Adjacent"
else if normal_origin == "Normal_unspecified_distance":
    source_origin = "Normal_unspecified_distance"
    display_3class = "Normal"
else:
    STOP and report the unmatched exact values
```

两分类字段直接使用原始 `tissue_group`，不得从三分类反向重建。

## 队列分布

| 数据集 | Normal细胞 | Adjacent细胞 | CRC细胞 | 旧来源证据 |
|---|---:|---:|---:|---|
| GSE132465 | 1,041 | 0 | 14,843 | Normal mucosa；距离未明确 |
| GSE144735 | 1,116 | 0 | 4,572 | Normal mucosa；距离未明确 |
| GSE161277 | 3,827 | 0 | 7,504 | normal；距离未明确 |
| GSE166555 | 16,497 | 0 | 12,310 | 样本标题为normal；距离未明确 |
| GSE178341 | 35,512 | 0 | 82,001 | 标本字段为N；距离未明确 |
| GSE200997 | 2,690 | 0 | 6,911 | 样本标题为Normal；距离未明确 |
| GSE231559 | 3,971 | 0 | 9,548 | Colon normal；距离未明确 |
| GSE245552 | 0 | 8,165 | 13,228 | colon adjacent tissue；明确癌旁 |

详细数量见 [tissue_source_audit.csv](tissue_source_audit.csv)。官方记录：

- [GSE132465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)
- [GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735)
- [GSE161277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161277)
- [GSE166555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555)
- [GSE178341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341)
- [GSE200997](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997)
- [GSE231559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231559)
- [GSE245552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552)

## 解释边界

- 两分类图回答全部正常来源与CRC来源上皮在共同表达轨迹上的分布。
- 三分类图沿用项目旧口径，回答Normal、明确Adjacent和CRC的分布。
- Normal是 `Normal_unspecified_distance`，不能写成远端正常或健康供者。
- Adjacent只来自GSE245552的4个捕获单元，来源与队列完全重合。任何Adjacent与其他组的差异都必须同时报告这一限制，并提供逐队列和宽上皮身份敏感性图。
- 三类在拟时序轴上的顺序不是同一患者的连续病程。

