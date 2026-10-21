# 参考文件使用边界

- `Overall_Monocle2_generated.R`：来源版本的完整R实现，开头仍含Mac预处理补丁。主入口`build_engine`会按已核实标志去掉这两段及调用，再保存PC的`work/fit_generated.R`。不要直接运行该参考脚本。
- `postprocessing_reference.py`：从旧第14步抽出的评分、动态和绘图函数。全局路径、Overall/Normal+CRC范围、旧数量及图例需按TASK_SPEC适配。旧函数的I/O读取块不是独立模型；PC可用全矩阵稀疏计算实现同公式。新的报告标签支持限制须显式实现。
- `Overall_Monocle2_session.txt`：来源机器历史包版本，不是PC环境确认。
- `iterative_DFS_equivalence.json`、`preparation_validation.json`及逐文库核对CSV：历史审核证据，只确认所列来源版本，不代表PC拟合完成。
- `monocle_Artistic-2.0.txt`：保留Monocle许可文本。

所有参考文件SHA在REFERENCE_MANIFEST.json中。不要原样执行旧Mac脚本，不上传整篇文献PDF或其他阶段的完整项目。
