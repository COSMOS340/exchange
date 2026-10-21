# 本机新增后处理与验证代码

全部在原生150,917个CRC来源上皮模型上运行；不改变输入细胞、CNV标签、排序基因、拟合参数或生成R引擎。计算失败与修复过程保留在日志目录。

- `run_fit_local.py`、`monitor_local.py`：单个拟合的文件锁、精确命令、目录、PID、起止和退出码；观察内存与原生进度。已有进程记录时拒绝静默覆盖。
- `validate_monocle_saved.R`：从真实保存对象独立遍历树，核验细胞/CNV、原始UMI、原生size factor、投影图距离、路径成员和两根相同几何。已通过。修复验证器的整数节点名格式，避免100000变为科学计数法；模型未改。
- `check_root_fallback.R`：合成图的实际根回退及加权距离核验，已通过。此类软件测试不是正式细胞下采样。
- `plot_main_crc.py`：13组共同轨迹主图与CNV组成/覆盖表。所有拟合细胞保留；只随机绘图顺序以减少遮挡。`--only <图名>`只重绘指定图。修复表格字体/间距、队列面板空白及PNG/PDF布局变化；最新PNG与每页PDF均已实际查看。
- `score_functions_crc.py`：全部CRC池的六套固定基因集评分，使用原生`nCount_input`分母及总体标准差；稀疏分批读取只用于I/O，不拆模型。已完成。
- `dynamics_crc.py`：每路径12箱，按建库>=10细胞、建库/标本/有支持报告标签逐层等权；不足3标签留空。缺少患者链接的标本另图。`--method Monocle2`与`--method Slingshot`可分别后处理保存对象，两者齐备后生成全路径配对及合并表。Monocle2的8条路径已完成，16组PNG/PDF均已实际查看；Slingshot等待其正式拟合保存。
- `validate_postprocess_crc.py`：独立的CSC密集面板算术及NumPy factorization/add.at汇总。`--stage scores-and-compositions`已在全部150,917细胞通过，最大评分差5.44e-15；`--stage dynamics --method Monocle2`已核验8条路径的分箱、逐层权重、空箱规则和标签内首末观测变化。`--stage dynamics`最终核验全方法路径配对，包括低覆盖配对，不挑最大相关。
- `check_slingshot_environment.R`：120个合成细胞的getLineages/getCurves/embedCurves/PT/权重/MST接口检查，已通过。
- `validate_slingshot_saved.R`：尚待正式Slingshot对象保存后执行；验证30维Harmony、全部权重/NA、元数据对应、导出曲线和UMAP投影前后拟时序与权重一致。
- `plot_slingshot_crc.py`：尚待执行，只展示正式embedCurves保存的UMAP曲线与全部细胞。
- `topology_sensitivity.R`：尚待Slingshot完成后串行执行；固定上游PCA/Harmony及继承簇，仅比较MST与逐队列留出结构，不重训Harmony或重拟合曲线。
- `render_pdf_audit.py`：用独立视觉审核环境PyMuPDF逐页渲染全部PDF；只生成待审清单，不会自动声称人工视觉通过。真实审阅另存`main_visual_review.json`和`Monocle2_dynamics_visual_review.json`，包含已实际查看文件的SHA256。

运行脚本以本目录反推交接根路径；不需要原来源机器路径。Python绘图可设置`MPLCONFIGDIR`为工作区可写缓存，数值后处理设置BLAS/OMP/MKL线程为1并低优先级运行。

# 环境与原生警告

`environment_bootstrap_explicit.txt`是初始Conda锁，包含临时igraph 2.2.3元数据；随后从CRAN官方归档源码安装实际igraph **2.2.2**，并添加zlib开发包与pkg-config。必须以R实际加载版本核验，不能只看Conda元数据。隔离环境R4.5.3、Monocle2 2.38.0、DDRTree 0.1.6、Slingshot 2.18.0。

Slingshot首轮因缺DelayedMatrixStats退出。环境v2只添加`DelayedMatrixStats 1.32.0`及`sparseMatrixStats 1.22.0`，原有包未替换。新增包精确下载清单为`Slingshot_environment_v2_additions.txt`，实际版本与完整合成API检查见日志。主Monocle2对象保留在原环境版本下的绑定记录中。

`wslconfig.proposed`为已获用户授权应用的48GB内存、24GB磁盘交换空间设置；用户已重启，实际配置已核验生效。

原生Monocle2内部`reduceDimension`设置seed2016；配置seed20261020不是所有底层随机数的证明。原生离散度非收敛与kmeans最大步数警告没有隐藏；目标函数5项、两次上升，末次相对变化6.412e-5满足原生tol0.001停止条件，但不证明全局最优。完整判断与限制见`warning_review.json`。

igraph源码：https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.2.2.tar.gz 。源码SHA、实际环境和所有修复记录见日志；本目录不存输入、大型模型、认证令牌或系统密钥。
