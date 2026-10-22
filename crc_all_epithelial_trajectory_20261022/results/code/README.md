# 全上皮共同轨迹执行代码

本目录所有输出以本任务根目录为准，统计范围为全部223,736个已冻结上皮细胞。代码从上一轮已验证实现移植并进行独立范围/来源测试，具体差异见任务根目录`IMPLEMENTATION_CHANGES.md`。

- `pipeline_common.py`：固定路径、输入清单、全细胞校验、来源精确关联和汇总规则。
- `prepare_all.py`、`build_engine.py`：全量准备、原元数据文本保留及兼容引擎派生。
- `select_endpoint_roots.R`：共同树上的端点覆盖选择，不接收组织/CNV/功能列。
- `run_fit_all.py`：单重任务锁、参数/源码绑定、实际PID、开始结束、退出码和GNU time峰值内存记录。
- `exact_projection_acceleration.R`、`exact_projection_bvh.cpp`、`fit_slingshot_exact_entry.R`：上一轮的精确候选线段搜索，所有原始顶点和细胞保留。
- `validate_source_all.py`、`validate_monocle_saved.R`、`validate_slingshot_saved.R`：独立原字段检查及已保存真实对象的检查；不拟合或删除细胞。
- `score_functions_crc.py`、`dynamics_crc.py`：文件名保留来源方便追溯，本任务内已明确改为AllEpithelial和223,736细胞；六套活动在整个新分析池重新标准化，路径汇总使用capture_id。
- `render_pdf_audit.py`只负责渲染，不代表已经人工查看；实际视觉记录和科学验收需另行完成。

具体软件版本、实际解释器路径和运行状态以`04_日志和检查`中的实际记录为准。不得把准备完成、拟合退出0、存在RDS或打包成功当作科学验收完成。

## 动态基因与完整证据链

- `aggregate_dynamic_expression.py`保存全部基因的capture/bin与specimen/bin表达；`validate_dynamic_expression.py`从原始CSC以另一种直接细胞加权方式复算所有均值、检出和总体标准化。
- `dynamic_model.py`/`analyze_dynamic_genes.py`实现冻结的标本截距WLS与链接簇CR1；`validate_dynamic_model_full.py`使用完整标本dummy的独立QR设计，核对两根每个真实基因系数、协方差秩/F、R的P/BH、幅度及选择。
- `derive_dynamic_modules.py`/`cluster_profiles.py`保存全部K=2..8、整链接单位bootstrap与精确全基因silhouette。`validate_source_profiles_modules.py`独立复算来源/共享身份/队列留一描述、全部主模块silhouette及代表基因；`dynamic_sensitivity.py`保留每个根/全部队列的结果与失败，`validate_dynamic_sensitivity.py`另以完整dummy QR复算所有留一系数与幅度。
- `fetch_GO_sources.py`/`build_GO_snapshot.py`保存官方来源及可复算的人类注释/原始本体快照；`enrich_dynamic_modules.py`在测量合格背景上检验全部BP条目，`validate_GO_independent.py`从原始注释独立重建关系并用R复算全部概率/BH。显示去冗余规则先冻结于`GO_display_policy.json`，不参与选基因或聚类。
- `module_path_dynamics.py`保存全细胞固定模块/代表基因值与每个返回路径的完整单位汇总；`validate_module_activity_paths.py`从原始CSC逐基因Z直接重算每个细胞，并核对全部路径成员、权重阈值、逐层均值、IQR、代理标本及NA。
- `compare_previous_CRC.py`只读取/复制旧CRC表，精确匹配150,917细胞，保留全部根/State/路径/CNV对照；不重新拟合旧对象。

## 执行次序与实际运行记录

在任务根运行`AllEpithelium_PC.py verify-input`、`prepare`，再分别`fit --method Monocle2`和`fit --method Slingshot`。正式重计算串行，不重复启动活动模型。已存在的原始输出会受到不覆盖保护；重试前先保存失败版本及原因。

实际Monocle2拟合代码提交为`bdcc8df4fdf927b9d2db11ad349c1e60da76b8d8`，Slingshot为`adeb6886ff58a7779dc36ac8e961eb6937112691`；完整命令、解释器、参数与SHA在对应binding/process日志。两者后的新图表/验证代码在独立后续提交中完成，不倒填为拟合时已使用的代码。

`post_monocle_serial.py`在真实对象保存后串行检查/计分/主图/全基因汇总与统计。`post_slingshot_serial.py`等待真实Slingshot保存再串行核验对象、生成全部路径活动、核对实际基因模型、根/队列/方法模块敏感性与MST诊断。实际早已启动的该监督进程编译自`97f240a`（9步骤）；后加的`plot_dynamic_sensitivity.py`步骤在本次运行中另行执行，并在日志保存实际命令，不能只凭最终脚本猜执行历史。

其后单独运行`validate_module_activity_paths.py`、`validate_dynamic_sensitivity.py`、`plot_dynamic_sensitivity.py`。逐步真实返回码、失败/重试与资源记录为准；本README不是完成证明。

路径默认为当前PC的已核验环境。另一台机器重算前需明确配置`pipeline_common.ENV`及相应Python/R前缀，根据冻结环境记录重建依赖；WSL环境二进制不应直接作为Mac环境。Mac可直接接收完整RDS和源表，原始输入按原Release及SHA另取。

## 视觉、文档与回传

各`plot_*`输出PNG/PDF与源表清单。图形修正支持`--only`的脚本只重绘指定图，旧版及实际问题保留在日志attempts。使用隔离视觉环境运行`render_pdf_audit.py`后，必须实际打开每张PNG和每页PDF渲染结果，并将具体观察及最终文件SHA写入`actual_visual_observations.json`；不允许由渲染成功自动填通过。

完成实际警告/收敛审阅后，`final_acceptance_all.py`检查所有独立数值证据、真实对象/哈希及逐图逐页记录；`build_deliverable_docs.py`生成A–H摘要和相对路径逐图指南。另保留完整方法/限制、软件原生日志、许可证、图形修正及失败记录。

提交最终代码/小型汇总/文档后，`AllEpithelium_PC.py pack-return --label NEW_VERSION`只打包通过验收的本次results和绑定侧文件，分片仅用于传输、每片1.5GiB。用`verify_return_archive.py RETURN_FOLDER`逐个读取所有ZIP成员验证CRC/字节/SHA，并核对分片拼接SHA。`publish_new_release.py`只创建新的结果分支Release，锁定打包代码提交，逐附件核对远端名称/字节/SHA后发布；保存的草稿状态支持续传，不覆盖旧Release或main。Mac正式接收验收仍需独立进行。
