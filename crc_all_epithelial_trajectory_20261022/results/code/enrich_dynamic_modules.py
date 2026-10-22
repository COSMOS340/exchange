"""Module BP overrepresentation against the measured detection-eligible universe."""
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from dynamic_model import bh_adjust
from pipeline_common import ROOT,TAB,LOG,save_json,digest


def main():
    folder=TAB/'dynamic_genes/Monocle2_primary'
    destination=folder/'dynamic_module_GO_BP_full.csv'
    if destination.exists():raise RuntimeError('GO results already exist; preserve/version before rerunning')
    policy=dict(created_before_enrichment=True,display_q=.05,max_terms_per_module=3,
        exclude_root_terms=['GO:0008150'],minimum_background_hits=10,
        maximum_background_hits=500,minimum_module_hits=3,
        redundancy='Ancestor/descendant or hit GeneID Jaccard >=0.5, greedy BH then P then GO_ID',
        statistical_tests_and_gene_module_selection_unchanged=True)
    policy_file=ROOT/'results/code/GO_display_policy.json'
    if policy_file.exists():assert json.loads(policy_file.read_text())==policy
    else:save_json(policy_file,policy)
    snapshot=TAB/'GO_snapshot'
    assert json.loads((snapshot/'GO_snapshot_provenance.json').read_text())['passed']
    coverage=pd.read_csv(TAB/'dynamic_genes/dynamic_gene_detection_coverage.csv')
    background_features=set(coverage.loc[coverage.detection_eligible,'gene'])
    mapping=pd.read_csv(snapshot/'measured_gene_symbol_mapping.csv',dtype=str,keep_default_na=False)
    status_by_gene=dict(zip(mapping.gene,mapping.status))
    gene_to_id={r.gene:(r.GeneID if r.GeneID else 'unmapped_feature::'+r.gene) for r in mapping.itertuples()}
    resolved=[r.GeneID for r in mapping.itertuples() if r.GeneID]
    assert len(resolved)==len(set(resolved)), 'Many measured features map to the same GeneID; audit before enrichment'
    background={gene_to_id[x] for x in background_features}
    reverse={gene_to_id[x]:x for x in background_features}
    terms=pd.read_csv(snapshot/'GO_BP_terms_and_ancestors.csv.gz',keep_default_na=False).set_index('GO_ID')
    associations=pd.read_csv(snapshot/'human_measured_gene_BP_associations.csv.gz',dtype={'GeneID':str,'GO_ID':str})
    term_sets=defaultdict(set)
    for row in associations.itertuples():
        if row.GeneID in background:term_sets[row.GO_ID].add(row.GeneID)
    annotated=set().union(*term_sets.values())
    assignments=pd.read_csv(folder/'dynamic_gene_module_assignment.csv')
    assignments=assignments[assignments.selected]
    assert set(assignments.gene)<=background_features
    rows=[];module_counts=[]
    for module,d in assignments.groupby('module',sort=True):
        members={gene_to_id[x] for x in d.gene}
        n=len(members);na=len(members&annotated);N=len(background);Na=len(annotated)
        module_counts.append(dict(module=module,module_measured_genes=len(d),module_units=n,
            module_BP_annotated_units=na,background_measured_genes=len(background_features),
            background_units=N,background_BP_annotated_units=Na))
        for go,term in terms.iterrows():
            universe=term_sets.get(go,set());hit=members&universe;k=len(hit);K=len(universe)
            p=float(hypergeom.sf(k-1,N,K,n)) if K else np.nan
            pa=float(hypergeom.sf(k-1,Na,K,na)) if K and na else np.nan
            rows.append(dict(module=module,GO_ID=go,term=term['name'],
                module_size=n,background_size=N,module_hits=k,background_hits=K,
                fold_enrichment=(k/n)/(K/N) if K else np.nan,p_value=p,
                BP_annotated_module_size=na,BP_annotated_background_size=Na,
                p_value_annotated_background=pa,hit_GeneIDs=';'.join(sorted(hit)),
                hit_gene_symbols=';'.join(sorted(reverse[x] for x in hit)),
                background_hit_gene_symbols=';'.join(sorted(reverse[x] for x in universe)),
                status='tested' if K else 'no_measured_eligible_background_members'))
    full=pd.DataFrame(rows)
    full['BH_q_value']=bh_adjust(full.p_value.to_numpy())
    full['BH_q_annotated_background']=bh_adjust(full.p_value_annotated_background.to_numpy())
    full.to_csv(destination,index=False)
    pd.DataFrame(module_counts).to_csv(folder/'dynamic_GO_background_coverage.csv',index=False)
    pd.DataFrame([dict(gene=x,universe_id=gene_to_id[x],BP_annotated=gene_to_id[x] in annotated,
                       mapping_status=status_by_gene[x]) for x in sorted(background_features)]).to_csv(
        folder/'dynamic_GO_background_genes.csv',index=False)
    display=[]
    for module,d in full.groupby('module',sort=True):
        candidates=d[d.BH_q_value.le(.05)&d.background_hits.between(10,500)&d.module_hits.ge(3)&~d.GO_ID.eq('GO:0008150')]
        chosen=[];chosen_hits=[]
        for row in candidates.sort_values(['BH_q_value','p_value','GO_ID']).to_dict('records'):
            go=row['GO_ID'];hit=set(row['hit_GeneIDs'].split(';'))
            ancestor=set(str(terms.loc[go,'ancestors']).split(';'))
            redundant=any(old['GO_ID'] in ancestor or go in set(str(terms.loc[old['GO_ID'],'ancestors']).split(';')) or
                          len(hit&old_hit)/len(hit|old_hit)>=.5 for old,old_hit in zip(chosen,chosen_hits))
            if redundant:continue
            chosen.append(row);chosen_hits.append(hit)
            if len(chosen)==3:break
        if chosen:
            for row in chosen:
                row['display_status']='significant_nonredundant'
                row['annotation_background_sensitive']=not (row['BH_q_annotated_background']<=.05)
                display.append(row)
        else:
            any_significant=bool(d.BH_q_value.le(.05).any())
            display.append(dict(module=module,GO_ID='',term='no significant term meets display criteria' if any_significant else 'no enriched GO term',
                                display_status='no_eligible_display_term' if any_significant else 'no_significant_enrichment'))
    pd.DataFrame(display).to_csv(folder/'dynamic_module_GO_BP_display.csv',index=False)
    save_json(folder/'dynamic_GO_enrichment_summary.json',dict(completed=True,
        background_detection_eligible_features=len(background_features),background_universe_units=len(background),
        BP_annotated_background_units=len(annotated),modules=module_counts,
        BP_terms_in_full_table=len(terms),multiple_testing='BH across all finite module-term P values',
        no_whole_human_genome_background=True,annotation_covered_background_sensitivity_saved=True,
        source_snapshot_sha256=digest(snapshot/'GO_snapshot_provenance.json'),
        display_policy_sha256=digest(policy_file),independent_validation_pending=True,
        interpretation='Overrepresentation of annotations in expression-derived modules; no causal or pathway-activity proof'))
    (ROOT/'results/code/licenses/GO_attribution.txt').write_text(
        'Gene Ontology Consortium; ontology data release 2026-07-26.\n'
        'Copyright 1999-2026 Gene Ontology Consortium. Licensed under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/\n'
        'Provided as-is, without warranty. Official source: https://current.geneontology.org/ontology/go-basic.obo\n'
        'Terms are linked by their persistent GO IDs. Frozen ontology and human NCBI annotation subsets are included with source SHA256.\n'
        'Citation policy: https://www.geneontology.org/docs/go-citation-policy/\n'
        'Original publication DOI: https://doi.org/10.1038/75556\n'
        'GO knowledgebase update DOI: https://doi.org/10.1093/nar/gkaf1292\n')
    print('GO module overrepresentation and annotation-coverage sensitivity saved.',flush=True)


if __name__=='__main__':main()
