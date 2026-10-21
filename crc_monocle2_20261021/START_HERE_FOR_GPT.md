# 给哥哥PC上的GPT/Codex

请接手 https://github.com/COSMOS340/exchange 中 `crc_monocle2_20261021` 的CRC上皮Monocle2任务。

先完整读取该目录AGENTS.md、README.md、TASK_SPEC.md、RUN_ON_PC.md及RETURN_PROTOCOL.md。读取实际环境，不猜盘符或软件路径。使用本PC计算，不在GitHub Actions运行。输入在Release `crc-monocle2-input-v1-20261021`；先校验SHA，再解压和运行。

主分析要求全部150,917个CRC来源上皮共同拟合一棵Monocle2树，无随机抽样，无分块独立建树。保留既有CNV低/中/高标签，在同一轨迹比较三组；不称所有CRC上皮均为恶性，不强设低CNV为起点。完成计算参考根和替代根、原始对象/表、主图、功能动态以及Slingshot同细胞对照。优先交付主Monocle2，遇到后续失败保留主结果并报告。

入口已经提供准备和拟合，但全部后处理还需你按TASK_SPEC实施。reference文件仅提供移植依据，不要原样执行旧路径或照抄旧细胞数。保存每项代码修改、实际环境、日志及输入/输出哈希；拟合完成标记不等于验收通过。全部PNG和PDF页面均要实际看图，修复遮挡、留白和裁切，再写逐图发现。无可靠发现或覆盖不足时按实说明。

OPTIONAL_ANALYSES只作扩展清单。无需确认的主任务继续完成；需要改变细胞范围、新增统计检验或算法的扩展先向用户确认。不要自动删除细胞或改用其他方法解决错误。

每个主要阶段更新results/RUN_SUMMARY.md，明确真正保存的进度。最终按RETURN_PROTOCOL提交代码回传分支，并将模型、表格、图和日志打包上传新的结果Release。保留本地输出，回报准确的commit和Release链接，方便Mac端读取、验收和归入原项目。
