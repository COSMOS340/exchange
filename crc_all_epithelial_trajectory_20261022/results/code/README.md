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
