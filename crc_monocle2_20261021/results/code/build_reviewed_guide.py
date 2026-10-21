"""Assemble the guide only from actually reviewed, hash-bound figure files."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
RESULT=ROOT/'results'
LOG=RESULT/'04_日志和检查'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    manifests=['main_figure_manifest.json','slingshot_figure_manifest.json','dynamics_manifest.json']
    figures=[]
    for name in manifests:
        figures.extend(json.loads((LOG/name).read_text()))
    reviewed={}
    for file in sorted(LOG.glob('*_visual_review.json')):
        for record in json.loads(file.read_text())['records']:
            assert record['actual_visual_inspection'] and record['status']=='passed'
            assert digest(RESULT/record['file'])==record['sha256'],record['file']
            assert record.get('observations'),record['file']
            record['review_record']=file.name
            reviewed[record['file']]=record
    actual={str(p.relative_to(RESULT)) for p in (RESULT/'01_图形').rglob('*')
            if p.suffix in ['.png','.pdf']}
    described={item[kind] for item in figures for kind in ['png','pdf']}
    assert actual==described and actual==set(reviewed)
    pages=json.loads((LOG/'pdf_page_manifest.json').read_text())
    for page in pages:
        record=reviewed[page['pdf']]
        assert page['page'] in record['pages_reviewed']
        assert (RESULT/page['rendered_png']).is_file()
        page.update(visual_review='passed_actual_view',review_record=record['review_record'],
                    pdf_sha256=record['sha256'])
    (LOG/'pdf_page_manifest.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2)+'\n')
    parts=['# 逐图阅读指南\n\n',
      f'本指南覆盖全部{len(figures)}组正式PNG/PDF；全部PNG与{len(pages)}页PDF已实际查看，审阅记录与当前文件SHA256逐一绑定。\n\n',
      '拟合池为150,917个CRC来源上皮，原CNV值和三分组保持不变。CRC来源不等于恶性；拟时序是表达结构的计算排序，不证明真实时间、克隆演化、恶性转化或耐药机制。\n\n',
      'State为Monocle2树区段，C编号为继承簇，宽身份为上游表达注释。共有142个报告患者标签，其中126个有对应支持；它们并非已验证的独立患者。动态先按每建库每箱至少10细胞，再建库/标本/报告标签逐层等权，主体每箱至少3个有支持标签，缺少患者链接的标本另图；空箱断线。四分位区间不是置信区间。\n\n',
      '每条路径独立缩放至0–1、12箱；路径共享树干，计数不能相加为独立细胞。箱间标签/标本组成可能改变，均值端点差不等于同一患者纵向进展。表达活动是全CRC池基因Z均值，不是代谢通量或新GSEA；CP10K使用原始nCount_input分母。\n\n',
      'Slingshot以继承30维Harmony拟合，embedCurves投影至原UMAP；高维拟时序与权重保留。精确线段索引保留全部细胞、全部原顶点和线段索引；软件原实现对照见日志，原生全CRC基准未算完，不声称两套正式全量拟合已逐值比较。没有新增P值、星号或动态基因显著性。\n\n',
      'MST敏感性使用冻结PCA/Harmony与继承簇，不重训Harmony或曲线；8个队列逐一留出中7项可计算，留出GSE178341后协方差奇异，保留为不可评估，不能声称全部队列留出均稳健。详见[敏感性源表](02_表格/CRC/topology_sensitivity.csv)及[警告审阅](04_日志和检查/warning_review.json)。\n\n']
    for item in figures:
        record=reviewed[item['png']]
        sources=list(item['sources'])
        axis=item.get('axis_and_color','')
        if '/dynamics/' in item['png']:
            proxy=item['name'].endswith('_proxy_specimens')
            stem=item['name'].removesuffix('_proxy_specimens').removesuffix('_supported_labels')
            sources=[f'02_表格/CRC/dynamics/{stem}_unsupported_specimens.csv.gz' if proxy else
                     f'02_表格/CRC/dynamics/{stem}_summary.csv',
                     f'02_表格/CRC/dynamics/{stem}_library_bins.csv.gz',
                     f'02_表格/CRC/dynamics/{stem}_specimen_bins.csv.gz',
                     f'02_表格/CRC/dynamics/{stem}_reported_label_bins.csv.gz',
                     f'02_表格/CRC/dynamics/{stem}_boundaries.csv']
            if not axis:
                axis='横轴为本路径0–1缩放拟时序的12箱；纵轴依次为连续RNA CNV平均绝对log2、三组/周期比例、六套平均基因Z活动及log1p原始总UMI。'+('各彩线为个别标本均值，身份以源表为准；缺箱断线。' if proxy else '蓝线为报告标签等权均值，阴影为标签间四分位区间，覆盖不足留空。')
        elif item['name'].startswith('12_Slingshot') and not axis:
            axis='横、纵轴为继承UMAP坐标，无物理时间单位；细胞颜色为冻结CNV三组，线条为30维拟合后投影的路径。二维交叉不证明转换。'
        assert axis,item['name']
        for source in sources:assert (RESULT/source).is_file(),source
        parts.extend([f'## {item["name"]}\n\n',
          f'![{item["name"]}]({item["png"]})\n\n[PDF]({item["pdf"]})\n\n',
          f'问题：{item["question"]}\n\n轴和颜色：{axis}\n\n',
          f'统计单位：{item["statistical_unit"]}。\n\n实际观察及限制：{record["observations"]}\n\n',
          '源表：'+'、'.join(f'[{Path(p).name}]({p})' for p in sources)+'。\n\n'])
    (RESULT/'01_逐图阅读指南.md').write_text(''.join(parts).rstrip()+'\n')
    report=dict(status='passed_all_actual_visual_reviews_bound_to_files',png_count=len(figures),
       pdf_count=len(figures),pdf_page_count=len(pages),all_formal_files_covered=True,
       all_file_hashes_match_actual_review=True,records=list(reviewed.values()))
    (LOG/'visual_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(f'Guide and actual visual audit assembled: {len(figures)} PNGs / {len(pages)} PDF pages')


if __name__=='__main__':
    main()
