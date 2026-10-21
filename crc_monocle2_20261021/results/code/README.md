# 本机新增后处理与验证代码

全部在原生150,917个CRC来源上皮模型上运行；不改变输入细胞、CNV标签、排序基因、拟合参数或生成R引擎。计算失败与修复过程保留在日志目录。

- `run_fit_local.py`、`monitor_local.py`：单个拟合的文件锁、精确命令、目录、PID、起止和退出码；观察内存与原生进度。已有进程记录时拒绝静默覆盖。
- `finish_saved_postprocess.py`：已完成6个串行步骤，全部退出0；核验真实对象、MST敏感性、绘图/动态、全方法独立复算及逐页PDF渲染。MST脚本成功执行不等于所有诊断可计算：GSE178341留出协方差奇异，源表保留NA和错误。记录每步命令、代码SHA及退出码；失败即停止，不会重启拟合、自动通过人工视觉审核或发布结果。
- `validate_monocle_saved.R`：从真实保存对象独立遍历树，核验细胞/CNV、原始UMI、原生size factor、投影图距离、路径成员和两根相同几何。已通过。修复验证器的整数节点名格式，避免100000变为科学计数法；模型未改。
- `check_root_fallback.R`：合成图的实际根回退及加权距离核验，已通过。此类软件测试不是正式细胞下采样。
- `plot_main_crc.py`：13组共同轨迹主图与CNV组成/覆盖表。所有拟合细胞保留；只随机绘图顺序以减少遮挡。`--only <图名>`只重绘指定图。修复表格字体/间距、队列面板空白及PNG/PDF布局变化；最新PNG与每页PDF均已实际查看。
- `score_functions_crc.py`：全部CRC池的六套固定基因集评分，使用原生`nCount_input`分母及总体标准差；稀疏分批读取只用于I/O，不拆模型。已完成。
- `dynamics_crc.py`：每路径12箱，按建库>=10细胞、建库/标本/有支持报告标签逐层等权；不足3标签留空。缺少患者链接的标本另图。`--method Monocle2`与`--method Slingshot`可分别后处理保存对象，两者齐备后生成全路径配对及合并表。Monocle2的8条与Slingshot的20条路径均已完成；56组动态PNG/PDF已逐张/逐页实际查看。
- `validate_postprocess_crc.py`：独立的CSC密集面板算术及NumPy factorization/add.at汇总。`--stage scores-and-compositions`已在全部150,917细胞通过，最大评分差5.44e-15；`--stage dynamics --method Monocle2`已核验8条路径的分箱、逐层权重、空箱规则和标签内首末观测变化。`--stage dynamics`已核验全部28条路径、80个全方法路径配对，包括低覆盖配对，不挑最大相关；逐层均值最大误差3.02e-14。
- `validate_unit_hierarchy.py`：冻结建库→标本→报告标签及支持标志的一致性检查通过；163建库、152标本、142报告标签、其中126有对应支持，无冲突；不证明独立患者身份。
- `validate_cnv_exports.py`：逐行比较冻结CNV十进制文本与三组指示值；150,917细胞导出值完全一致，最大十进制差0。冻结元数据及基因去向的原文件副本另存表格目录并逐字节核验。
- `inspect_slingshot_initial_mst.R`：轻量MST初始化诊断，所有150,917细胞和30维Harmony、两个起点各10条候选路径；不执行曲线拟合。原生初始化投影工作量可达数百亿点段配对，不能把150曲线点参数当150个拟合细胞。
- `check_slingshot_environment.R`：120个合成细胞的getLineages/getCurves/embedCurves/PT/权重/MST接口检查，已通过。
- `validate_slingshot_saved.R`：已重新加载4个真实对象，30维Harmony、全部权重/NA、元数据对应、导出曲线及UMAP投影前后PT/权重均核验通过；PT导出误差最多5.26e-13，投影前后对象值完全相同。
- `plot_slingshot_crc.py`：已生成并实际审阅两根正式embedCurves UMAP曲线及全部细胞的PNG/PDF。
- `topology_sensitivity.R`：已执行固定PCA/Harmony及8个队列留出MST诊断，7项留出可计算，GSE178341留出协方差奇异、保留NA与错误。未重训Harmony/曲线，未改距离/抽样来凑齐成功结果。
- `render_pdf_audit.py`：用独立视觉审核环境PyMuPDF逐页渲染全部PDF；只生成待审清单，不会自动声称人工视觉通过。实际审阅分方法记录在4份`*_visual_review.json`；71张PNG与71页PDF均已查看，与全部142个正式图文件SHA256绑定。

运行脚本以本目录反推交接根路径；不需要原来源机器路径。Python绘图可设置`MPLCONFIGDIR`为工作区可写缓存，数值后处理设置BLAS/OMP/MKL线程为1并低优先级运行。

# 环境与原生警告

`environment_bootstrap_explicit.txt`是初始Conda锁，包含临时igraph 2.2.3元数据；随后从CRAN官方归档源码安装实际igraph **2.2.2**，并添加zlib开发包与pkg-config。必须以R实际加载版本核验，不能只看Conda元数据。隔离环境R4.5.3、Monocle2 2.38.0、DDRTree 0.1.6、Slingshot 2.18.0。

Slingshot首轮因缺DelayedMatrixStats退出。环境v2只添加`DelayedMatrixStats 1.32.0`及`sparseMatrixStats 1.22.0`，原有包未替换。新增包精确下载清单为`Slingshot_environment_v2_additions.txt`，实际版本与完整合成API检查见日志。主Monocle2对象保留在原环境版本下的绑定记录中。

`wslconfig.proposed`为已获用户授权应用的48GB内存、24GB磁盘交换空间设置；用户已重启，实际配置已核验生效。

原生Monocle2内部`reduceDimension`设置seed2016；配置seed20261020不是所有底层随机数的证明。原生离散度非收敛与kmeans最大步数警告没有隐藏；目标函数5项、两次上升，末次相对变化6.412e-5满足原生tol0.001停止条件，但不证明全局最优。完整判断与限制见`warning_review.json`。

igraph源码：https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.2.2.tar.gz 。源码SHA、实际环境和所有修复记录见日志；本目录不存输入、大型模型、认证令牌或系统密钥。

后处理Harmony/PCA读取显式声明`header=TRUE`，避免数字列名0..29被自动识别为一行数据。已保留轻量诊断首轮失败与修复记录，重读150,917×30成功；正式生成R引擎原本已显式声明表头，其SHA与拟合绑定保持不变。

# Slingshot精确线段搜索v3

原生Slingshot的初始化会把全部输入细胞投影到含数万顶点的折线上，再压缩为配置中的150个曲线点；本机原生尝试计算超过3小时仍未保存曲线。已核查princurve 2.1.6官方C++源码，点×线段×维度循环有明确边界，不能把耗时称为卡死。

`exact_projection_bvh.cpp`对原始线段建立区间包围盒，只有保守距离下界已经超过当前最近距离时才跳过候选。保留全部输入细胞、原曲线顶点和线段索引；投影标量运算、相同距离时选最早原线段、原排序及弦长累加与princurve参考一致。包围盒使用long double下界及宽松浮点误差保护；不作近似最近邻、不抽样，也不修改已安装包。代码含princurve GPL-3归属说明，许可证在`licenses/`。

`test_exact_projection_acceleration.R`通过27个随机、重复、回折、交叉、等距和微小转折软件情形，坐标/PT/距离最大差0；360个合成细胞的双根联合曲线PT/权重差0。单次合成投影基准0.882秒降至0.070秒，不能据此承诺全流程加速倍数。`test_exact_entry_end_to_end.R`用2,400个从零生成的软件细胞运行实际冻结引擎，两根PT、权重、投影PT及UMAP曲线均与原始入口完全相同。它们不是正式生物学下采样；样例对象仅在被忽略的`.runtime/`中。

`fit_slingshot_exact_entry.R`仅在本R进程的Slingshot导入绑定安装已验证索引，随后执行完全相同的冻结R引擎与配置。小于1000顶点或非零stretch调用仍用原始投影。`run_slingshot_exact_local.py`核对输入/工作清单、已验证源码SHA、包版本与Rcpp，记录全部有效源码、环境、命令和退出状态；持有同一拟合锁，拒绝覆盖已有尝试或对象。投影审计逐次记录完整行数、保留顶点数和实际候选比较数。

最初的“删除共线顶点”原型未通过回折情形，已拒绝并移到`rejected_polyline/`，未用于正式计算；不能将几何相同当作原生拟时序相同。正式v3仅使用已通过上述零差异检查的全部线段精确搜索。旧原生尝试的进程、日志和退出状态已单独归档；Monocle2原始断点和两根模型保留原绑定。


# 最终交付记录工具

- `record_viewed_pdfs.py`：仅在通过图像查看工具实际查看指定PDF渲染页后，按显式图名记录审阅声明；不执行或推断视觉审核，没有“全图自动通过”选项。
- `assemble_slingshot_reviews.py`：要求42张新图的PNG/PDF均有明确实际查看声明、观察与匹配SHA，才整理审阅记录，附真实路径细胞/支持箱/代理标本数量；原Monocle2审核保留。
- `build_reviewed_guide.py`：从真实观察生成71组逐图指南；要求所有正式PNG/PDF都有审阅记录、当前SHA一致且所有PDF页面被实际审阅，源表路径存在；不以渲染成功冒充实际查看。
- `capture_final_provenance.py`：记录最终实际R版本/Python依赖、交付commit和所有有效源码SHA，复核7个真实对象及原3个Monocle2对象字节不变。历史preflight只表示其记录当时状态。
- `finalize_acceptance.py`：综合全部独立检查、源码绑定、实际视觉记录与警告审阅，明确列出GSE178341留出不可评估，接受范围为本地核心计算/交付，不声称全部敏感性均通过或Mac已经接收验收。
- `verify_return_archive.py`：流式检查完整ZIP与分片SHA/字节数、全部ZIP成员CRC、RESULT_MANIFEST每个文件SHA及路径；不解包大型模型、不把传输成功叫分析成功。首次已在旧小型OOM诊断包验证工具，完整结果包须单独运行。

Slingshot v3最终用时5950.6秒、退出0，两根各10条路径、150,917细胞、有限PT未分配0。原生对象未保存完整逐轮收敛状态；退出0和maxit30不证明收敛或全局最优。详见最终RUN_SUMMARY、warning_review及本地验收JSON。
