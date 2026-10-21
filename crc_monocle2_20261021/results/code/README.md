# 本机新增代码与环境记录

- `run_fit_local.py`：串行监督一个拟合入口，保存实际命令、目录、PID、起止时间与退出码，并使用文件锁避免双启动。第一次Monocle2尝试已用此入口执行，失败现场保留在日志目录。
- `validate_monocle_saved.R`：拟合保存后的独立图遍历、CNV/细胞对应和主/替代根几何核验。当前没有真实模型，仅完成语法检查，尚未执行数值验收。
- `wslconfig.proposed`：本次已授权应用的WSL配置内容，保留原有设置。实际修改、备份路径及哈希见日志目录。
- `environment_bootstrap_explicit.txt`：隔离环境的初始Conda包锁定清单。注意它包含临时依赖igraph 2.2.3，不能据此直接启动正式拟合。随后从CRAN官方归档安装了指定的igraph 2.2.2，并补装了zlib开发包和pkg-config；R实际版本核验通过后才运行。

igraph源码地址：`https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.2.2.tar.gz`。源码SHA256及实际环境见`04_日志和检查/environment_provenance.json`、`R_environment_validation.log`、`python_freeze.txt`和`conda_packages_explicit.txt`。Conda元数据仍记录初始igraph 2.2.3；R实际加载的版本为源码安装的2.2.2，不能只看Conda元数据推断版本。

不在此目录存放输入数据、大型模型或认证信息。尚未实现或执行功能动态及绘图后处理。
