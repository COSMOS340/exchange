# 全上皮共同轨迹：PC结果摘要

PC全量计算、独立数值核验和实际视觉检查已完成；完整返回包尚需单独传输核验/发布记录，Mac接收验收仍待执行。旧CRC-only结果与对象保留。

## A. 输入与来源核查

完成：223,736细胞、15,635基因、497,539,060非零计数一次共同分析；无随机细胞抽样、按组织拆树或新增QC。输入ZIP/全部文件SHA、矩阵方向、基因/细胞顺序及原字段逐字核对通过。两分类72,819正常来源/150,917CRC，三分类64,654 Normal/8,165 Adjacent/150,917CRC，0未匹配。8队列/242 captures/229 specimens/149报告标签（非核实独立患者）。

## B. 全量共同Monocle2

完成：1,998个冻结排序基因；146中心/145边；223,736节点投影树连通无环；两根同几何，均7 State、4根到叶路径。全部细胞的加权根距离独立复算最大误差0，size factor独立误差0。保留原始checkpoint、主根和替代根3个完整RDS。

- primary：请求端点Y_125/State7；实际 `GSM5075664::p009n2:TCCGGGATCCGGACTG`，C2/State7。
- alternative：请求端点Y_54/State1；实际 `GSE178341::C157_N_0_0_0_c1_v3_id-ACCCAAAAGCTGGCCT`，C4/State1。

原生GLM/参数离散度非收敛、Quick-TRANSfer上限等警告完整保留；15,621基因保存离散度均有限，不等于所有原生拟合已收敛。DDRTree最后相对目标变化0.000131，含一次上升，非全局最优证明。Monocle约56m22s、峰值RSS45.9GiB；详见原生日志与警告审阅。

## C. 来源/CNV及旧CRC-only比较

完成：两/三分类共同图与分面、PT/State/CNV/身份/队列/Phase/UMI、pooled及capture→specimen等权、6个共享宽身份、全部8队列描述性留一；全部150,917旧CRC细胞精确对齐，原CNV/身份字段一致，全部根/路径配对保存。

可保留：组织来源、身份和队列在共同结构上分布不均且重叠。主根/替代根全池PT秩相关约−0.696，存在分支关系而非全树单纯反向。旧/新CRC排序四配对rho约0.587/−0.458/0.326/−0.528；CRC内部CNV趋势随换根/几何变化，不能据全池图证明统一Normal→Adjacent→CRC恶性进展。Adjacent与单一队列完全混杂。

## D. 动态基因、模块与GO

- primary：11,902基因可测，1,787达到冻结q/幅度阈值；3,116标本-bin观测，208个可提供标本内变化的标本、122报告/未知链接簇，29/30箱支持。完整229标本资格表保留。
- alternative：11,902基因可测，2,205达到冻结q/幅度阈值；2,655标本-bin观测，209个可提供标本内变化的标本、123报告/未知链接簇，27/30箱支持。完整229标本资格表保留。

完成：主根K=2模块（607/1,180基因），候选K=2..8全基因silhouette、整链接单位bootstrap10次；根选择交集1,472基因、Jaccard=0.584；全部8队列留一、全部方法/根/路径的固定模块与代表基因敏感性保存。完整热图和两张全部成员放大图、所有基因/箱/GO全表可追溯。

GO背景11,902测量且检测合格基因，另有10,713 BP注释覆盖背景敏感性；完整47,948模块/术语组合行，其中25,148有限P值纳入BH，其余没有合格背景成员而记不可测。逐项独立复算通过。显示注释为：

- Module1：oxidative phosphorylation；RNA splicing；purine ribonucleoside triphosphate metabolic process。
- Module2：actin cytoskeleton organization；epithelium development；regulation of actin filament-based process。

动态统计为条件于推断轨迹的标本内关联，未传播轨迹不确定性，也非细胞独立样本检验。GO不参与模块选择且不证明通路/蛋白活性或因果。

## E. 冻结六套Hallmark

完成：同223,736全池重算基因总体均值/SD，按nCount_input计算CP10K；命中和变异资格按冻结规则审核。每条真实路径12箱，主体capture→specimen→有支持报告标签等权（>=3标签/bin），缺少患者链接的代理标本单独保留；CNV连续/三组、六活动、cycling及UMI全部输出，空箱不补0、不跨接。

## F. 同细胞Slingshot对照

- primary：同223,736细胞/30维Harmony，实际起点C2；11路径，223,736至少一条PT有限、0完全未分配。全部原权重/NA保留，UMAP投影前后PT和权重逐项一致。
- alternative：同223,736细胞/30维Harmony，实际起点C4；10路径，223,736至少一条PT有限、0完全未分配。全部原权重/NA保留，UMAP投影前后PT和权重逐项一致。

Slingshot两根连同原生UMAP投影实际耗时9,206.388秒（2小时33分26秒），峰值RSS3,111,264KiB（约2.97GiB），退出码0。原生日志未见警告，但对象没有保存迭代次数或hasConverged；最大30次、阈值0.001不等于已证实以阈值停止。原生曲线延伸及二维交叉保留。

- MST PCA_vs_Harmony：边Jaccard=0.882。
- MST leave_out_GSE132465：边Jaccard=0.882。
- MST leave_out_GSE144735：边Jaccard=0.882。
- MST leave_out_GSE161277：边Jaccard=0.882。
- MST leave_out_GSE166555：边Jaccard=0.882。
- MST leave_out_GSE178341：不可评估，原因：system is computationally singular: reciprocal condition number = 2.0322e-18。
- MST leave_out_GSE200997：边Jaccard=0.882。
- MST leave_out_GSE231559：边Jaccard=1.000。
- MST leave_out_GSE245552：边Jaccard=0.882。

逐队列动态模型敏感性（固定共同图与检测背景）：

- GSE132465：可测，选择1,775基因，选择Jaccard=0.952，主成员分区ARI=0.982。
- GSE144735：可测，选择1,769基因，选择Jaccard=0.978，主成员分区ARI=0.998。
- GSE161277：可测，选择1,784基因，选择Jaccard=0.984，主成员分区ARI=0.993。
- GSE166555：可测，选择1,866基因，选择Jaccard=0.886，主成员分区ARI=0.962。
- GSE178341：可测，选择1,812基因，选择Jaccard=0.791，主成员分区ARI=0.925。
- GSE200997：可测，选择1,818基因，选择Jaccard=0.967，主成员分区ARI=0.998。
- GSE231559：可测，选择1,803基因，选择Jaccard=0.975，主成员分区ARI=0.998。
- GSE245552：可测，选择1,768基因，选择Jaccard=0.968，主成员分区ARI=0.995。

所有路径对照均报告，不挑最高相关路径；MST层面PCA/Harmony与逐队列剔除另存，失败/不可评估保持NA及原因。没有独立外部生物学验证，也没有全量未加速Slingshot重复基线。

## G. 单位、未做与解释

细胞只作结构显示；主趋势使用捕获/标本/有支持报告标签分层等权。未把149标签称为149独立患者，未加细胞级显著性星号，未新增GSEA、DNA CNV、代谢通量、真实时间、恶性克隆或因果推断。固定图上的队列留一不等于重训Harmony/DDRTree/全部曲线。

## H. 验收与交付

完成：最终数值证据链、7个完整RDS哈希、全部97组PNG/PDF及逐页实际观察通过PC验收。全部源表、生成代码、参数/软件版本、拟合提交绑定、失败重试与图形修正均保留。

失败/修复记录：GO下载早期HTTP403后按完整字节范围重试并验证SHA；独立来源验证器的重复分组键已修复后全量通过；独立R复算器把CSV的True/False读为字符而中断，已显式解析布尔值后对全部真实基因重跑QR/F/BH核验通过，原始统计模型未改动，失败与重试日志保留；图注/实际根/色标/空缺区间/孤立均值点等实际视觉问题已保留旧版后修复并重新打开核对。原生优化警告属于保留限制，不伪称修复收敛。其他实际不可评估项目以对应源表状态和日志为准。

仍待：完整ZIP/分片本地CRC与SHA传输核验、新Release上传后远端名称/字节/SHA及精确代码提交核验；Mac下载与第14步接收复核。后续回传收据单独记录这些状态，不将PC验收等同于Mac正式接收。

入口：[逐图阅读指南](01_逐图阅读指南.md)、[方法与限制](02_方法与限制.md)、`04_日志和检查/numeric_acceptance.json`、`04_日志和检查/visual_audit.json`。
