# PC计算进度与阶段验收

当前主Monocle2已完成；Slingshot仍在全150,917细胞主起点拟合中，整个任务尚未验收或发布最终结果包。

## 已完成

- A 输入/环境：输入Release `crc-monocle2-input-v1-20261021`，ZIP与全部7个输入文件哈希、CRC范围与UMI/元数据顺序核验通过。冻结150,917个CRC来源上皮；8队列、152标本、163建库、142个报告患者标签（不等于已验证独立患者）。Low/Intermediate/High=41,686/46,755/62,476；原始连续CNV逐值不变；正常/癌旁72,819细胞保留在来源输入，不参加本次树拟合。
- B 全量Monocle2：第二轮用时1370.3秒、退出0。全部细胞共同拟合一棵2维DDRTree，原生size factor/离散度，排序基因1,996；`~dataset`、maxIter100、tol0.001；未抽样或拆模型。3个原始RDS、全细胞坐标/PT/State、中心树、投影图边、路径和目标函数已保存。
- 实际主根 `GSM5075666::p009t2:TGTTTGTAGGATTTCC`（C1、State1）；替代根 `GSM5075665::p009t1:TTGTGTTTCATGCAGT`（C11、State2）。两根各7个State、4条路径，共享完全相同DDRTree几何；144个中心、143条中心树边；投影图150,917节点、150,916边。
- C 主图：13组PNG/PDF均已实际查看，每个PDF1页，标签/布局问题修复后重新查看。原完整坐标与全部细胞保留。独立计算到实际根的图距离与导出PT最大差0；原生size factor最大相对差0，路径成员也独立重建通过。图排序/根回退合成软件测试在本机通过。
- D 功能评分和Monocle2动态：六套固定Hallmark均可评估；原始`nCount_input`分母、全CRC池基因总体标准化。全部150,917细胞独立重算评分最大差5.44e-15，CNV组成独立汇总通过。8条Monocle2路径各12箱，建库>=10细胞后依次建库/标本/有支持报告标签等权；不足3标签留空，未支持患者链接另图。逐层分箱及标签内首末观测变化独立复算通过。16组动态PNG/PDF均实际查看。
- 目前共29组PNG/PDF；逐图真实观察、单位和限制见`01_逐图阅读指南.md`。全部实际视觉记录带文件SHA256，见`main_visual_review.json`和`Monocle2_dynamics_visual_review.json`。

## 仍在进行/待做

- E Slingshot：正式第二轮始于2026-10-21 18:27:44（UTC+08）。同细胞、继承30维Harmony、同C1/C11参考、maxit30、approx_points150，保持原入口；当前主起点尚未保存对象。等待其正常完成，勿重复启动。
- Slingshot保存后：重新加载核验全部权重/NA、UMAP embedCurves与PT不变，制作两起点UMAP图和逐路径动态，完整报告所有方法/路径配对。PCA/Harmony及逐队列留出仅作MST敏感性，未重训Harmony或重拟合曲线。
- F 最终整体数值验收、剩余图实际视觉审核、完整回传包与新的结果Release尚未完成。

## 原生警告与限制

原生离散度GLM/参数拟合出现非收敛、NaN及步长截断警告；保存的15,617个离散度估计均有限，但不证明参数最优。`Removing 48 outliers`是离散度拟合的基因异常点，不是删除细胞。kmeans Quick-TRANSfer达到7,545,850步上限。DDRTree目标函数5项、两次上升，最后相对变化6.4120576e-5满足原生tol0.001停止条件，不能声称全局最优。完整判读见`warning_review.json`。

配置seed20261020；Monocle2原生reduceDimension内部设置seed2016；未改原生行为，不把配置seed当全部底层随机数已控制的证据。

三组CNV在共同树中广泛重叠；路径的CNV和活动方向不一致，更换起点会改变排序（全细胞PT Spearman0.635213）。State3只有50细胞，统计尤其脆弱。E2F/MYC与周期、UMI并行展示；功能评分是表达活动而非代谢通量、新GSEA或因果机制。报告标签并非确认的独立患者，箱间组成变化不能等同同一患者进展；共享树干路径不是独立重复。没有新增显著性检验、星号、BEAM或tradeSeq。

## 环境、失败与修复记录

WSL2 Ubuntu26.04，R4.5.3、Python3.13.16、Monocle2 2.38.0、DDRTree0.1.6、实际igraph2.2.2、Slingshot2.18.0。igraph从CRAN源码安装，Conda元数据保留初始2.2.3，必须以R实际加载版本为准。完整版本、sessionInfo、源代码及输入绑定见日志和`code/README.md`。

第一次Monocle2在WSL原31GiB内存的size factor阶段OOM（R -9），无断点。用户授权修改48GB RAM+24GB磁盘交换空间并自行重启，验证后第二轮成功，原失败日志保留。Slingshot首轮缺DelayedMatrixStats，在新环境v2仅补装两个依赖且不替换现有包，合成全API检查通过后才重启该对照；主模型未重算或改动。独立后处理验证器修复了pandas factorize列名保留问题；失败日志保留，复算通过。

输入commit `4bb17d4ffc23912ad1cd82ad87f4760d5e0f4955`；回传分支 `pc-results-20261021-wsl`。原始输入与已有失败诊断包保留。当前完整成功结果包尚未生成，不能用旧OOM诊断包代替本轮结果。
