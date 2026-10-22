# 给PC端GPT/Codex的完整指令

请接手 `COSMOS340/exchange` 仓库中的 `crc_all_epithelial_trajectory_20261022`。

本任务重新进行第14步主分析：将输入Release中的223,736个正常来源与CRC来源上皮细胞放入同一Monocle2轨迹。旧目录 `crc_monocle2_20261021` 的150,917个CRC内部拟时序保留为敏感性分析，不能覆盖，也不能把旧CRC筛选后的work目录当成本任务输入。

严格执行以下要求：

1. 下载并校验Release `crc-monocle2-input-v1-20261021` 的 `input_data.zip` 与 `INPUT_ARCHIVE.json`。
2. 读取 `data/Overall` 的完整矩阵和元数据，核对15,635 × 223,736、cell_id顺序和所有冻结标签。
3. 按 `SOURCE_TISSUE_AUDIT.md` 建立两分类与三分类来源字段。只允许使用文件中列明的精确规则，不从样本名后缀自行推断。
4. 在全部223,736个细胞上一次拟合Monocle2 DDRTree。排序基因和起点选择不使用组织来源、CNV或功能值。
5. 同一坐标展示：全部正常来源对CRC；沿用旧来源审计的Normal、Adjacent、CRC三类；CNV低中高；宽上皮身份；队列；细胞周期和UMI。
6. 制作 `REFERENCE_FIGURE_REQUIREMENTS.md` 规定的拟时序动态基因热图，包含动态模块、代表基因、模块GO条目和完整源表。
7. 用同一批细胞运行Slingshot作方法敏感性，不能用另一个细胞范围制造一致性。
8. 以建库/标本等权进行趋势摘要与敏感性检查；不能将223,736个细胞当作223,736个独立生物学重复。
9. 完成数值与视觉验收后，上传新结果Release `crc-all-epithelial-trajectory-results-v1-20261022` 或在标签已存在时使用递增的新标签，禁止覆盖旧Release。

无需向用户重复确认上述主任务。若必须改变细胞范围、输入标签、核心算法或统计单位，应停止并说明具体证据，不得静默修改。

