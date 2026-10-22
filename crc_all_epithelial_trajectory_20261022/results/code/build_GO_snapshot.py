"""Extract a reproducible human annotation/ontology snapshot from official data."""
import csv
import gzip
import json
from collections import defaultdict

import pandas as pd

from pipeline_common import ROOT, TAB, LOG, digest, save_json


def read_ontology(path):
    terms, aliases, current = {}, {}, None
    version = None
    def commit():
        if current and 'id' in current and current.get('namespace')=='biological_process' and not current.get('obsolete'):
            terms[current['id']] = current.copy()
            for alias in current['aliases']:
                aliases[alias] = current['id']
    for line in path.read_text().splitlines():
        if line.startswith('data-version: '):
            version=line.split(': ',1)[1]
        if line=='[Term]':
            commit();current=dict(parents=[],aliases=[])
        elif line.startswith('['):
            commit();current=None
        elif current is not None:
            for key in ['id','name','namespace']:
                if line.startswith(key+': '):current[key]=line.split(': ',1)[1]
            if line.startswith('is_a: '):current['parents'].append(line.split()[1])
            if line.startswith('relationship: part_of '):current['parents'].append(line.split()[2])
            if line.startswith('alt_id: '):current['aliases'].append(line.split()[1])
            if line=='is_obsolete: true':current['obsolete']=True
    commit()
    memo, active = {}, set()
    def ancestors(term):
        if term in memo:return memo[term]
        if term in active:raise ValueError('Unexpected cycle in GO basic BP graph')
        active.add(term); result={term}
        for parent in terms[term]['parents']:
            parent=aliases.get(parent,parent)
            if parent in terms:result.update(ancestors(parent))
        active.remove(term);memo[term]=result;return result
    for term in terms:ancestors(term)
    return terms,aliases,memo,version


def main():
    downloads=json.loads((LOG/'GO_source_downloads.json').read_text())
    if downloads['status']!='all_downloads_complete_pending_annotation_validation':
        raise RuntimeError('Wait for complete immutable annotation downloads')
    folder=ROOT/'data/GO_sources'
    for item in downloads['files']:
        assert digest(folder/item['name'])==item['sha256']
    genes=pd.read_csv(ROOT/'work/AllEpithelial/genes.csv').gene.tolist()
    exact,alias=defaultdict(set),defaultdict(set)
    with gzip.open(folder/'Homo_sapiens.gene_info.gz','rt') as f:
        reader=csv.DictReader(f,delimiter='\t')
        for row in reader:
            if row['#tax_id']!='9606':continue
            official=row['Symbol_from_nomenclature_authority']
            if official=='-':official=row['Symbol']
            identifier=row['GeneID']
            exact[official].add(identifier)
            alternatives=set(row['Synonyms'].split('|')) | {row['Symbol']}
            for name in alternatives-{'-',official}:
                alias[name].add(identifier)
    rows, ids = [], set()
    for gene in genes:
        candidates=exact.get(gene,set())
        status='exact_official' if candidates else 'unambiguous_alias'
        if not candidates:candidates=alias.get(gene,set())
        if len(candidates)!=1:status='unmapped' if not candidates else 'ambiguous_symbol'
        identifier=next(iter(candidates)) if len(candidates)==1 else None
        if identifier is not None:ids.add(identifier)
        rows.append(dict(gene=gene,GeneID=identifier,status=status,
                         candidates=';'.join(sorted(candidates))))
    mapping=pd.DataFrame(rows)
    output=TAB/'GO_snapshot';output.mkdir(exist_ok=False)
    mapping.to_csv(output/'measured_gene_symbol_mapping.csv',index=False)
    terms,aliases,ancestors,version=read_ontology(folder/'go-basic.obo')
    direct=defaultdict(set);counts=defaultdict(int)
    with gzip.open(folder/'gene2go.gz','rt') as f:
        for line in f:
            if not line.startswith('9606\t'):continue
            values=line.rstrip('\r\n').split('\t')
            counts['human_rows']+=1
            _,identifier,go,evidence,qualifier,term,pubmed,category=values
            if identifier not in ids:continue
            counts['measured_mapped_rows']+=1
            if 'NOT' in qualifier.split('|') or evidence=='ND' or category!='Process':
                counts['excluded_NOT_ND_or_non_BP']+=1;continue
            go=aliases.get(go,go)
            if go not in terms:
                counts['unknown_or_obsolete_BP']+=1;continue
            direct[identifier].add(go)
    associations=[]
    for identifier,gos in direct.items():
        expanded=set().union(*(ancestors[go] for go in gos))
        associations.extend(dict(GeneID=identifier,GO_ID=go,direct=go in gos) for go in sorted(expanded))
    assert len(associations)>0 and counts['human_rows']>0, 'Human BP annotation parsing produced no associations'
    pd.DataFrame(associations).to_csv(output/'human_measured_gene_BP_associations.csv.gz',index=False)
    pd.DataFrame([dict(GO_ID=t,name=d['name'],parents=';'.join(d['parents']),
                       ancestors=';'.join(sorted(ancestors[t]))) for t,d in terms.items()]).to_csv(
        output/'GO_BP_terms_and_ancestors.csv.gz',index=False)
    save_json(output/'GO_snapshot_provenance.json',dict(passed=True,ontology_version=version,
        input_gene_count=len(genes),mapping_status=mapping.status.value_counts().to_dict(),
        mapped_GeneIDs=len(ids),BP_annotated_GeneIDs=len(direct),BP_terms=len(terms),counts=dict(counts),
        source_download_manifest_sha256=digest(LOG/'GO_source_downloads.json'),
        positive_propagation_relations=['is_a','part_of'],NOT_ND_excluded=True,
        documentation=downloads['sources_documentation'],
        citation_and_license='https://www.geneontology.org/docs/go-citation-policy/'))
    print('Official human GO snapshot complete:',len(associations),'gene/term associations;',version,flush=True)


if __name__=='__main__':
    main()
