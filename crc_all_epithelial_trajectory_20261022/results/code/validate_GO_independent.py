"""Rebuild GO from raw OBO/human rows; use R for independent tail/BH math."""
import csv,gzip,json,re,shutil,subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
from pipeline_common import ROOT,TAB,LOG,digest,save_json,ENV


def main():
    snap=TAB/'GO_snapshot';folder=TAB/'dynamic_genes/Monocle2_primary'
    original=ROOT/'data/GO_sources/Homo_sapiens.gene_info.gz'
    frozen=snap/'Homo_sapiens.gene_info_original.gz'
    if not frozen.exists():shutil.copy2(original,frozen)
    assert digest(original)==digest(frozen)
    mapping=pd.read_csv(snap/'measured_gene_symbol_mapping.csv',dtype=str,keep_default_na=False)
    official,aliases=defaultdict(set),defaultdict(set)
    with gzip.open(frozen,'rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['#tax_id']!='9606':continue
            symbol=row['Symbol_from_nomenclature_authority']
            symbol=row['Symbol'] if symbol=='-' else symbol
            official[symbol].add(row['GeneID'])
            for alternate in set(row['Synonyms'].split('|')+[row['Symbol']])-{'-',symbol}:
                aliases[alternate].add(row['GeneID'])
    for row in mapping.itertuples():
        matches=official.get(row.gene,set()) or aliases.get(row.gene,set())
        assert set(filter(None,row.candidates.split(';')))==matches
        assert row.GeneID==(next(iter(matches)) if len(matches)==1 else '')
    with gzip.open(snap/'go-basic_original.obo.gz','rt') as f:text=f.read()
    terms={};alternate={}
    for block in text.split('[Term]\n')[1:]:
        block=block.split('\n[',1)[0]
        pairs=defaultdict(list)
        for line in block.splitlines():
            if ': ' in line:
                key,value=line.split(': ',1);pairs[key].append(value)
        if pairs.get('namespace')!=['biological_process'] or pairs.get('is_obsolete')==['true']:continue
        go=pairs['id'][0]
        parents={x.split()[0] for x in pairs.get('is_a',[])}
        parents.update(x.split()[1] for x in pairs.get('relationship',[]) if x.startswith('part_of '))
        terms[go]=parents
        for alias in pairs.get('alt_id',[]):alternate[alias]=go
    expanded={}
    for go in terms:
        # Iterative traversal, independent of the author's recursive memoization.
        todo=[go];seen=set()
        while todo:
            raw=todo.pop();node=alternate.get(raw,raw)
            if node in seen or node not in terms:continue
            seen.add(node);todo.extend(terms[node])
        expanded[go]=seen
    saved_terms=pd.read_csv(snap/'GO_BP_terms_and_ancestors.csv.gz',keep_default_na=False)
    assert set(saved_terms.GO_ID)==set(terms)
    for row in saved_terms.itertuples():assert set(row.ancestors.split(';'))==expanded[row.GO_ID]
    direct=defaultdict(set)
    with gzip.open(snap/'human_measured_gene2go_original_rows.tsv.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if row['Category']!='Process' or row['Evidence']=='ND' or 'NOT' in row['Qualifier'].split('|'):continue
            go=alternate.get(row['GO_ID'],row['GO_ID'])
            if go in terms:direct[row['GeneID']].add(go)
    rebuilt={(gene,ancestor) for gene,gos in direct.items() for go in gos for ancestor in expanded[go]}
    saved=pd.read_csv(snap/'human_measured_gene_BP_associations.csv.gz',dtype=str)
    assert rebuilt==set(zip(saved.GeneID,saved.GO_ID))
    detection=pd.read_csv(TAB/'dynamic_genes/dynamic_gene_detection_coverage.csv')
    bg_features=set(detection.loc[detection.detection_eligible,'gene'])
    id_by_gene={x.gene:x.GeneID or 'unmapped_feature::'+x.gene for x in mapping.itertuples()}
    universe={id_by_gene[g] for g in bg_features}
    term_members=defaultdict(set)
    for gene,go in rebuilt:
        if gene in universe:term_members[go].add(gene)
    annotated=set().union(*term_members.values())
    assignment=pd.read_csv(folder/'dynamic_gene_module_assignment.csv')
    modules={m:{id_by_gene[g] for g in d.gene} for m,d in assignment[assignment.selected].groupby('module')}
    full=pd.read_csv(folder/'dynamic_module_GO_BP_full.csv',keep_default_na=False)
    for row in full.itertuples():
        hits=modules[row.module]&term_members[row.GO_ID]
        assert row.module_hits==len(hits) and row.background_hits==len(term_members[row.GO_ID])
        assert row.background_size==len(universe) and row.module_size==len(modules[row.module])
        assert set(filter(None,row.hit_GeneIDs.split(';')))==hits
        assert row.BP_annotated_background_size==len(annotated)
        assert row.BP_annotated_module_size==len(modules[row.module]&annotated)
    numeric=LOG/'GO_independent_numeric_input.csv'
    full[['module_hits','background_size','background_hits','module_size','BP_annotated_background_size','BP_annotated_module_size']].to_csv(numeric,index=False)
    result=LOG/'GO_independent_R_probabilities.csv'
    subprocess.run([str(ENV/'bin/Rscript'),str(ROOT/'results/code/validate_GO_probabilities.R'),str(numeric),str(result)],check=True)
    reference=pd.read_csv(result)
    comparisons={}
    for key in ['p_value','BH_q_value','p_value_annotated_background','BH_q_annotated_background']:
        actual=pd.to_numeric(full[key],errors='coerce').to_numpy();expected=reference[key].to_numpy()
        assert np.array_equal(np.isfinite(actual),np.isfinite(expected))
        assert np.allclose(actual,expected,rtol=2e-9,atol=1e-300,equal_nan=True),key
        comparisons[key]=float(np.nanmax(np.abs(actual-expected)))
    display=pd.read_csv(folder/'dynamic_module_GO_BP_display.csv')
    for module,d in display.groupby('module'):
        chosen=[];chosen_hits=[]
        for row in d.itertuples():
            assert row.BH_q_value<=.05 and 10<=row.background_hits<=500 and row.module_hits>=3
            assert row.GO_ID!='GO:0008150'
            hit=modules[module]&term_members[row.GO_ID]
            for previous,old in zip(chosen,chosen_hits):
                assert previous not in expanded[row.GO_ID] and row.GO_ID not in expanded[previous]
                assert len(old&hit)/len(old|hit)<.5
            chosen.append(row.GO_ID);chosen_hits.append(hit)
        assert len(chosen)<=3
    save_json(LOG/'GO_independent_validation.json',dict(passed=True,
        mapping_genes_checked=len(mapping),ontology_BP_terms=len(terms),original_human_rows_rebuilt=True,
        exact_gene_term_pairs=len(rebuilt),all_module_term_tests_checked=len(full),
        background_features=len(universe),R_phyper_and_p_adjust_max_absolute_errors=comparisons,
        original_gene_info_sha256=digest(frozen),display_policy_checked=True))
    print('Independent raw GO, universe, all term counts and R tail/BH validation passed.',flush=True)


if __name__=='__main__':main()
