# WSL重启后继续

当前没有模型在运行。2026-10-21第一次Monocle2在原生size factor阶段触发Linux OOM；R退出信号为SIGKILL（-9），没有DDRTree断点或模型结果。

已获用户授权，将Windows用户的`.wslconfig`更新为`memory=48GB`、`swap=24GB`，保留其原有网络等设置，并备份原文件。没有执行WSL重启。配置修改记录和备份路径见`04_日志和检查/wsl_hardware_and_proposal.json`。

用户方便时，先保存其他WSL任务，在**Windows PowerShell**执行：

```powershell
wsl --shutdown
```

该命令会停止所有WSL发行版。重新打开Ubuntu/Codex后，回到本聊天发送“已重启，继续”。无需重新下载输入、重装环境或重做准备。

接手检查：

1. 用`free -h`确认48GB内存和24GB交换空间配置实际生效。核对其他同项目进程及日志，避免双启动。
2. 项目目录是`/home/anran/hehe/exchange/crc_monocle2_20261021`；R/Python隔离环境是同目录的`.venv`。先激活现有环境，验证真实版本，不升级或重装现有环境。
3. 输入7个文件及`WORK_MANIFEST.json`已经验证；`work/preparation_summary.json`记录150,917细胞、1,996排序基因、C1/C11参考群。再次核验绑定后直接使用准备结果；不要重新执行`prepare`覆盖。
4. 检查第一轮失败诊断及对象目录。当前没有DDRTree断点，不能声称断点续算。把首轮Monocle2日志、绑定和进程记录归档到明确的新尝试目录，保留原文件和失败包；新尝试必须有独立记录，不能静默覆盖首轮证据。
5. 激活方式：`source /home/anran/miniforge3/etc/profile.d/conda.sh`，然后`conda activate /home/anran/hehe/exchange/crc_monocle2_20261021/.venv`。之后用`python results/code/run_fit_local.py Monocle2`串行启动；该监督入口发现旧进程记录会拒绝覆盖，必须先完成第4项。
6. 主模型保存后完成独立数值核验和主图，再按TASK_SPEC实施Slingshot、功能动态、图形审阅及交付。`results/code/validate_monocle_saved.R`已准备，当前仅做语法检查，没有真实对象可供执行；不应把该文件存在当作验收通过。

48GB内存和24GB交换空间不是后续拟合必定成功的保证；遇到新错误先保留日志和检查实际资源，不改变细胞范围或擅改原生算法。

Microsoft说明：https://learn.microsoft.com/en-us/windows/wsl/wsl-config
