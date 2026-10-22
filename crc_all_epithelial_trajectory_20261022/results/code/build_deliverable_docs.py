"""Build the final figure guide and A-H readout from verified saved results."""
import json
from pathlib import Path
import pandas as pd
from pipeline_common import ROOT,TAB,LOG,save_json,digest
from final_acceptance_all import visual_acceptance


def read(name):return json.loads((LOG/name).read_text())


def question(name):
    if name.startswith('30_') or name.startswith('31_'):return '哪些基因呈现共同的动态形状，模块注释与组成同时如何变化？'
    if name.startswith('32_'):return '候选模块数的分离程度和整链接单位bootstrap稳定性支持哪个K？'
    if name.startswith('33_'):return '保留主模块成员时，换根后的表达形状及选择重叠如何？'
    if name.startswith('34_'):return '逐一剔除全部队列后，固定共同图上的模块结果是否保留？'
    if name.startswith('40_'):return '相同CRC细胞在旧CRC-only与新共同树上的排序是否一致？'
    if name.startswith('41_'):return '旧、新拟合的State区段有哪些精确细胞对应？'
    if name.startswith('42_'):return '相同CRC细胞的CNV趋势对共同建树与换根是否敏感？'
    if name.startswith('50_'):return '全部路径上固定模块与代表基因的等权表达趋势是否相同？'
    if name.startswith('27_'):return '连续CNV色阶、UMI和细胞周期组成怎样影响共同树的直观判断？'
    if name.startswith('13_'):return '三种组织来源在同一个Slingshot全量拟合和投影中如何分布？'
    return '本图描述的共同结构、组成与覆盖是否支持所报告的结果？'


def encoding(name):
    if 'heatmap' in name or '_zoom' in name:
        return '横轴为主根30箱相对拟时序；纵轴为全部选中基因，白线分模块。红/蓝为逐基因调整样条Z（±2.5仅显示限幅），灰列为不足覆盖；上方五条为capture→specimen等权组成，右侧为代表基因及GO/BH q。'
    if name.startswith('50_'):
        return '各面板横轴为每条路径自身的0..1拟时序，纵轴包含全部返回路径；色值为有支持报告标签等权均值，模块单位为成员基因Z、代表基因为log1p(CP10K)，灰色为不足3标签覆盖。每特征跨方法/根共用色标。'
    if 'supported_labels' in name or 'proxy_specimens' in name:
        return '各面板横轴为路径内12箱相对拟时序；纵轴依次为连续CNV、三组比例、六套Hallmark平均基因Z、cycling比例与log1p(UMI)。主体为等报告标签均值/IQR，代理图每线为单个标本；NA断开，阴影不是置信区间。'
    if 'UMAP' in name or name.startswith('13_'):
        return '横纵轴为继承状态UMAP；完整细胞按图例着色，曲线来自30维Harmony拟合的原生投影，二维交叉不等于转化。来源分面保留完整背景与所有同一拟合曲线。'
    if name.startswith('40_'):return '横轴旧CRC-only、纵轴新共同树CRC子集的全池缩放拟时序；四面板展示所有根配对，hexbin聚合全部150,917细胞，色条为log刻度细胞数。'
    if name.startswith('41_'):return '纵轴旧State、横轴新State；每格为匹配细胞数，蓝色色阶为旧State内比例。原生State在各自拟合内换根不变，因此四根配对列联表经核对一致后只显示一次。'
    if name.startswith('42_'):return '横轴各次拟合全池的相对拟时序；纵轴为相同CRC细胞的原CNV，颜色区分旧/新拟合，实线capture→specimen等权均值、虚线pooled、阴影标本间IQR。'
    if name.startswith('32_'):return '横轴候选K；左图完整选中基因平均silhouette，右图10次完整链接单位bootstrap中位ARI，横虚线0.75，纵虚线标示依冻结规则选择的K。'
    if name.startswith('33_') or name.startswith('34_'):return '横轴为固定共同图上的相对拟时序，纵轴为固定主模块成员平均基因Z；颜色仅标示模块。换根图分别显示两根；队列留一图实线为留一模型，点虚线为完整池参考。'
    if 'composition' in name:return '横轴为既有CNV三组，纵轴State或路径；色阶0..1。左右区分pooled与capture→specimen等权，覆盖表列出细胞/标本数。共享树干路径覆盖不能相加。'
    if 'distributions' in name:return '横轴各根30箱相对拟时序，纵轴为该来源内每箱概率质量；颜色对应Normal/Adjacent/CRC，分别显示pooled与capture→specimen等权。真实0不代表缺失。'
    if 'leaveout' in name:return '横轴相对拟时序，纵轴capture→specimen等权来源比例；各面板依次排除一个队列，图例区分三种来源。固定树/表示，不重拟合；剔除唯一Adjacent队列后其比例为0。'
    if name.startswith('26_'):return '每行一宽身份（含全池），左右两根；横轴统一0..1，纵轴相对RNA CNV，每行共用尺度。圆点为至少3标本箱的capture→specimen等权均值，阴影标本间IQR；孤立点及NA间隙保留。'
    if name.startswith('27_'):return '四面板共用完整DDRTree坐标/树：连续CNV全范围、明确1%-99%颜色范围、log1p输入UMI与继承Phase。只有右上颜色限幅，原始数值与CNV分组未变。'
    if name=='10_root_orientation':return '前两面板共用DDRTree几何与PT色标，黑星是实际根细胞，标题区分请求State与实际根身份；右图为同细胞两根PT的全部点/描述性秩相关。'
    return '横纵轴为同一全量DDRTree坐标；完整黑树与细胞范围保留。类别颜色、连续色条、分面标题及n按图中图例解释；灰点为完整池背景，不是额外拟合。'


def main():
    numeric=read('numeric_acceptance.json');assert numeric['passed']
    records,observations=visual_acceptance()
    lines=['# 全上皮共同轨迹：逐图阅读指南','',
        '图按结果包内相对路径嵌入；正式PNG已实际打开、PDF全部页面已渲染并实际打开。观察与最终文件SHA绑定于 `04_日志和检查/actual_visual_observations.json`，完整核验见 `visual_audit.json`。',
        '', '所有图的组织来源均按冻结审计：Normal为取材距离未明确的正常组织，Adjacent为明确癌旁（仅GSE245552），不含已证实健康供者分类。共8队列、242捕获单元、229标本；149报告标签不是已核实的独立患者数。图中PT为计算参考距离/排序，不是真实时间。',
        '', '全量共同树上的来源、CNV、身份与队列分布须联合解读；共享输入的方法一致性不是独立验证，IQR不是置信区间，NA不补0。进一步方法及原生警告见 [方法与限制](02_方法与限制.md)。','']
    for name,row in sorted(records.items()):
        unit=row.get('statistical_unit')
        if not unit:
            unit=('Capture -> specimen -> supported reported label equal means; unsupported proxies separate' if name.startswith('50_') else
                'Specimen-intercept WLS; linkage-cluster covariance; all selected genes displayed' if name.startswith(('30_','31_')) else
                'Fixed-graph descriptive/module sensitivity; statistical unit stated in source tables')
        lines += ['## '+name,'',f'![{name}]({row["png"]})','',
            '**问题：** '+row.get('question',question(name)),'',
            '**轴与颜色：** '+encoding(name),'', '**统计单位：** '+unit+'。','',
            '**实际观察：** '+observations[name]['observation'],'',
            '**队列构成与解释边界：** 原始池的队列/来源审计不变；本图局部资格和路径/bin的细胞、capture、specimen、报告标签及队列数以源表为准。Adjacent与GSE245552无法分开识别；宽身份构成、共享树干和换根会影响排序，不能据图认定癌变方向或独立患者进展。','',
            '**源表与PDF：** '+', '.join(f'[{Path(p).name}]({p})' for p in row.get('sources',[]))+f'；[PDF]({row["pdf"]})。','']
    (ROOT/'results/01_逐图阅读指南.md').write_text('\n'.join(lines)+'\n')
    mono=read('Monocle2_independent_validation.json');sling=read('Slingshot_independent_validation.json')
    roots=json.loads((TAB/'Monocle2_roots.json').read_text())
    sensitivity=read('dynamic_sensitivity_summary.json');model={}
    for root in ['primary','alternative']:model[root]=json.loads((TAB/f'dynamic_genes/Monocle2_{root}/dynamic_model_summary.json').read_text())
    go=pd.read_csv(TAB/'dynamic_genes/Monocle2_primary/dynamic_module_GO_BP_display.csv')
    function=pd.read_csv(TAB/'function_gene_coverage.csv')
    summary=['# 全上皮共同轨迹：PC结果摘要','',
        'PC全量计算、独立数值核验和实际视觉检查已完成；完整返回包尚需单独传输核验/发布记录，Mac接收验收仍待执行。旧CRC-only结果与对象保留。','',
        '## A. 输入与来源核查','',
        '完成：223,736细胞、15,635基因、497,539,060非零计数一次共同分析；无随机细胞抽样、按组织拆树或新增QC。输入ZIP/全部文件SHA、矩阵方向、基因/细胞顺序及原字段逐字核对通过。两分类72,819正常来源/150,917CRC，三分类64,654 Normal/8,165 Adjacent/150,917CRC，0未匹配。8队列/242 captures/229 specimens/149报告标签（非核实独立患者）。','',
        '## B. 全量共同Monocle2','',
        '完成：1,998个冻结排序基因；146中心/145边；223,736节点投影树连通无环；两根同几何，均7 State、4根到叶路径。全部细胞的加权根距离独立复算最大误差0，size factor独立误差0。保留原始checkpoint、主根和替代根3个完整RDS。','']
    for root in ['primary','alternative']:
        r=roots[root];summary.append(f'- {root}：请求端点{r["requested_endpoint_center"]}/State{r["requested_State"]}；实际 `{r["actual_root_cell"]}`，{r["actual_root_group"]}/State{r["actual_root_State"]}。')
    summary += ['', '原生GLM/参数离散度非收敛、Quick-TRANSfer上限等警告完整保留；15,621基因保存离散度均有限，不等于所有原生拟合已收敛。DDRTree最后相对目标变化0.000131，含一次上升，非全局最优证明。Monocle约56m22s、峰值RSS45.9GiB；详见原生日志与警告审阅。','',
        '## C. 来源/CNV及旧CRC-only比较','',
        '完成：两/三分类共同图与分面、PT/State/CNV/身份/队列/Phase/UMI、pooled及capture→specimen等权、6个共享宽身份、全部8队列描述性留一；全部150,917旧CRC细胞精确对齐，原CNV/身份字段一致，全部根/路径配对保存。','',
        '可保留：组织来源、身份和队列在共同结构上分布不均且重叠。主根/替代根全池PT秩相关约−0.696，存在分支关系而非全树单纯反向。旧/新CRC排序四配对rho约0.587/−0.458/0.326/−0.528；CRC内部CNV趋势随换根/几何变化，不能据全池图证明统一Normal→Adjacent→CRC恶性进展。Adjacent与单一队列完全混杂。','',
        '## D. 动态基因、模块与GO','']
    for root,m in model.items():summary.append(f'- {root}：{m["tested_genes"]:,}基因可测，{m["selected_genes"]:,}达到冻结q/幅度阈值；{m["observations"]:,}标本-bin观测，{m["specimens"]}个可提供标本内变化的标本、{m["linkage_clusters"]}报告/未知链接簇，{m["supported_bins"]}/30箱支持。完整229标本资格表保留。')
    summary += ['',f'完成：主根K=2模块（607/1,180基因），候选K=2..8全基因silhouette、整链接单位bootstrap10次；根选择交集{sensitivity["root_selection_overlap"]:,}基因、Jaccard={sensitivity["root_selection_Jaccard"]:.3f}；全部8队列留一、全部方法/根/路径的固定模块与代表基因敏感性保存。完整热图和两张全部成员放大图、所有基因/箱/GO全表可追溯。', '',
        'GO背景11,902测量且检测合格基因，另有10,713 BP注释覆盖背景敏感性；47,948术语检验逐项独立复算通过。显示注释为：']
    for module,d in go.groupby('module'):summary.append('- '+module+'：'+'；'.join(d.term.tolist())+'。')
    summary += ['', '动态统计为条件于推断轨迹的标本内关联，未传播轨迹不确定性，也非细胞独立样本检验。GO不参与模块选择且不证明通路/蛋白活性或因果。','',
        '## E. 冻结六套Hallmark','',
        '完成：同223,736全池重算基因总体均值/SD，按nCount_input计算CP10K；命中和变异资格按冻结规则审核。每条真实路径12箱，主体capture→specimen→有支持报告标签等权（>=3标签/bin），缺少患者链接的代理标本单独保留；CNV连续/三组、六活动、cycling及UMI全部输出，空箱不补0、不跨接。','',
        '## F. 同细胞Slingshot对照','']
    for root in ['primary','alternative']:
        s=sling[root];summary.append(f'- {root}：同{s["cells"]:,}细胞/{s["dimensions"]}维Harmony，实际起点{s["requested_root"]}；{s["paths"]}路径，{s["assigned_cells"]:,}至少一条PT有限、{s["unassigned_cells"]:,}完全未分配。全部原权重/NA保留，UMAP投影前后PT和权重逐项一致。')
    topology=pd.read_csv(TAB/'topology_sensitivity.csv')
    for row in topology.itertuples():
        summary.append(f'- MST {row.comparison}：'+(f'边Jaccard={row.edge_Jaccard:.3f}。' if pd.notna(row.edge_Jaccard) else '不可评估，原因：'+str(row.detail)+'。'))
    summary += ['', '逐队列动态模型敏感性（固定共同图与检测背景）：']
    for row in pd.read_csv(TAB/'dynamic_genes/sensitivity/all_cohort_leaveout_summary.csv').itertuples():
        summary.append('- '+row.excluded_dataset+'：'+(f'可测，选择{int(row.selected_genes):,}基因，选择Jaccard={row.selection_Jaccard:.3f}，主成员分区ARI={row.primary_genes_partition_ARI:.3f}。' if row.status=='assessed' else '不可评估，原因：'+str(row.reason)+'。'))
    summary += ['', '所有路径对照均报告，不挑最高相关路径；MST层面PCA/Harmony与逐队列剔除另存，失败/不可评估保持NA及原因。没有独立外部生物学验证，也没有全量未加速Slingshot重复基线。','',
        '## G. 单位、未做与解释','',
        '细胞只作结构显示；主趋势使用捕获/标本/有支持报告标签分层等权。未把149标签称为149独立患者，未加细胞级显著性星号，未新增GSEA、DNA CNV、代谢通量、真实时间、恶性克隆或因果推断。固定图上的队列留一不等于重训Harmony/DDRTree/全部曲线。','',
        '## H. 验收与交付','',
        f'完成：最终数值证据链、7个完整RDS哈希、全部{numeric["formal_figure_pairs"]}组PNG/PDF及逐页实际观察通过PC验收。全部源表、生成代码、参数/软件版本、拟合提交绑定、失败重试与图形修正均保留。', '',
        '失败/修复记录：GO下载早期HTTP403后按完整字节范围重试并验证SHA；独立来源验证器的重复分组键已修复后全量通过；图注/实际根/色标/空缺区间/孤立均值点等实际视觉问题已保留旧版后修复并重新打开核对。原生优化警告属于保留限制，不伪称修复收敛。其他实际不可评估项目以对应源表状态和日志为准。','',
        '仍待：完整ZIP/分片本地CRC与SHA传输核验、新Release上传后远端名称/字节/SHA及精确代码提交核验；Mac下载与第14步接收复核。后续回传收据单独记录这些状态，不将PC验收等同于Mac正式接收。','',
        '入口：[逐图阅读指南](01_逐图阅读指南.md)、[方法与限制](02_方法与限制.md)、`04_日志和检查/numeric_acceptance.json`、`04_日志和检查/visual_audit.json`。']
    (ROOT/'results/RUN_SUMMARY.md').write_text('\n'.join(summary)+'\n')
    save_json(LOG/'deliverable_document_hashes.json',{name:digest(ROOT/'results'/name) for name in ['RUN_SUMMARY.md','01_逐图阅读指南.md','02_方法与限制.md']})
    print('Final A-H summary and per-figure source-backed guide written.',flush=True)


if __name__=='__main__':
    main()
