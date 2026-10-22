"""Derive a new task engine, preserving proven fitting and graph algorithms."""
from pipeline_common import ROOT, LOG, digest, save_json


def once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Unexpected inherited engine fragment: ' + old[:80])
    return text.replace(old, new)


def main():
    source = ROOT / 'reference/CRC_verified_fit_engine.R'
    original = source.read_text()
    text = original
    start = text.index(' choose=function(label,object){')
    end = text.index(' summaries=list()', start)
    text = text[:start] + """ source(file.path(normalizePath(file.path(input,'../..')),'results/code/select_endpoint_roots.R'))
 center_tree_for_roots=minSpanningTree(base)
 nearest_for_roots=base@auxOrderingData[['DDRTree']]$pr_graph_cell_proj_closest_vertex
 root_decisions=select_endpoint_roots(center_tree_for_roots,
   V(center_tree_for_roots)$name[nearest_for_roots[,1]],
   m[,c('cell_id','specimen_id','dataset','coarse_state','broad_identity')],
   as.character(pData(base)$State))
 fwrite(root_decisions$coverage,file.path(tab,'Monocle2_endpoint_root_coverage.csv'))
 fwrite(root_decisions$identity_coverage,file.path(tab,'Monocle2_endpoint_identity_coverage.csv'))
 write_json(root_decisions[setdiff(names(root_decisions),c('coverage','identity_coverage'))],
   file.path(tab,'Monocle2_endpoint_root_decisions.json'),pretty=TRUE,auto_unbox=TRUE)
""" + text[end:]
    text = once(text, "    requested=if(name=='primary')spec$root else spec$alternative\n    state=choose(requested,base)",
                "    decision=root_decisions[[name]]\n    requested=decision$requested_root_group\n    state=decision$state")
    text = once(text, 'summaries[[name]]=list(requested_root_group=requested,requested_State=state,',
                'summaries[[name]]=list(requested_root_group=requested,requested_State=state,\n       requested_endpoint_center=decision$center,endpoint_supporting_specimens=decision$supporting_specimens,\n       endpoint_support_threshold_met=decision$support_threshold_met,')
    text = once(text, ' roots=list(primary=spec$root,alternative=spec$alternative)',
                " mono_roots=fromJSON(file.path(tab,'Monocle2_roots.json'))\n roots=list(primary=mono_roots$primary$actual_root_group,alternative=mono_roots$alternative$actual_root_group)\n stopifnot(all(unlist(roots) %in% m$coarse_state))\n write_json(roots,file.path(tab,'Slingshot_actual_Monocle2_reference_groups.json'),pretty=TRUE,auto_unbox=TRUE)")
    # Only root policy and metadata differ. All algorithms before options(warn=1)
    # and the native DDRTree fit call remain byte-identical.
    boundary = original.index('options(warn=1)')
    if text[:boundary] != original[:boundary]:
        raise ValueError('Inherited compatibility or sparse graph code changed')
    native_fit = "cds=reduceDimension(cds,max_components=2,reduction_method='DDRTree',\n      residualModelFormulaStr='~dataset',maxIter=cfg$monocle_max_iter,tol=cfg$monocle_tolerance)"
    if native_fit not in text:
        raise ValueError('Native full-cell DDRTree call changed')
    destination = ROOT / 'work/fit_generated.R'
    destination.write_text(text)
    save_json(LOG / 'engine_derivation.json', dict(
        original_verified_engine_sha256=digest(source), derived_engine_sha256=digest(destination),
        original_graph_and_compatibility_functions_byte_identical=True,
        native_DDRTree_call_unchanged=True,
        changes=['Root requests selected after common graph by endpoint/specimen coverage',
                 'Slingshot reference groups read from actual saved Monocle2 roots'],
        cell_selection_in_engine=False, block='AllEpithelial',
        native_fit_parameters_changed=False))


if __name__ == '__main__':
    main()
