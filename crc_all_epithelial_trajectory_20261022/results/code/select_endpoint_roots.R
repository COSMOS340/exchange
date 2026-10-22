# Root selection sees only topology, IDs and coverage. Tissue/CNV/function
# columns are deliberately excluded at the call site.
select_endpoint_roots <- function(center_tree, nearest_center, metadata,
                                  base_states, minimum_cells=5L,
                                  minimum_specimens=5L) {
  stopifnot(igraph::is_tree(center_tree),length(nearest_center)==nrow(metadata),
            length(base_states)==nrow(metadata),!anyDuplicated(metadata$cell_id))
  ids=igraph::V(center_tree)$name
  leaves=sort(ids[igraph::degree(center_tree)==1L])
  stopifnot(length(leaves)>=2L)
  allowed=c('cell_id','specimen_id','dataset','coarse_state','broad_identity')
  stopifnot(identical(names(metadata),allowed))
  rows=list();identity_rows=list()
  total=table(base_states)
  for(leaf in leaves){
    inside=nearest_center==leaf
    if(!any(inside)){
      rows[[leaf]]=data.frame(center=leaf,cells=0L,specimens=0L,
        supporting_specimens=0L,datasets=0L,requested_State=NA_character_,
        requested_root_group=NA_character_,broad_identities=0L)
      next
    }
    d=metadata[inside,,drop=FALSE]
    counts=table(d$specimen_id)
    state_counts=table(as.character(base_states[inside]))
    st=data.frame(state=names(state_counts),n=as.integer(state_counts))
    st$fraction=st$n/as.numeric(total[st$state])
    st=st[order(-st$fraction,-st$n,as.integer(st$state)),]
    groups=sort(table(d$coarse_state),decreasing=TRUE)
    rows[[leaf]]=data.frame(center=leaf,cells=nrow(d),specimens=length(counts),
      supporting_specimens=sum(counts>=minimum_cells),datasets=length(unique(d$dataset)),
      requested_State=st$state[1],requested_root_group=names(groups)[1],
      broad_identities=length(unique(d$broad_identity)))
    by_identity=aggregate(cell_id~broad_identity+specimen_id,data=d,FUN=length)
    names(by_identity)[3]='cells';by_identity$center=leaf
    identity_rows[[leaf]]=by_identity
  }
  coverage=do.call(rbind,rows);rownames(coverage)=NULL
  eligible=coverage[coverage$cells>0 & !is.na(coverage$requested_State),,drop=FALSE]
  stopifnot(nrow(eligible)>=2L)
  supported=eligible[eligible$supporting_specimens>=minimum_specimens,,drop=FALSE]
  primary_pool=if(nrow(supported)>0) supported else eligible
  primary_pool=primary_pool[order(-primary_pool$supporting_specimens,
    -primary_pool$datasets,-primary_pool$specimens,-primary_pool$cells,primary_pool$center),]
  primary=primary_pool[1,,drop=FALSE]
  alternatives=eligible[eligible$center!=primary$center,,drop=FALSE]
  support_alt=alternatives[alternatives$supporting_specimens>=minimum_specimens,,drop=FALSE]
  if(nrow(support_alt)>0) alternatives=support_alt
  distinct=alternatives[alternatives$requested_State!=primary$requested_State,,drop=FALSE]
  if(nrow(distinct)>0) alternatives=distinct
  alternatives$distance_from_primary=as.numeric(igraph::distances(center_tree,
    v=primary$center,to=alternatives$center,weights=igraph::E(center_tree)$weight))
  alternatives=alternatives[order(-alternatives$distance_from_primary,
    -alternatives$supporting_specimens,-alternatives$datasets,alternatives$center),]
  alternative=alternatives[1,,drop=FALSE]
  choice=function(d) list(center=as.character(d$center),state=as.character(d$requested_State),
    requested_root_group=as.character(d$requested_root_group),
    supporting_specimens=as.integer(d$supporting_specimens),
    support_threshold_met=as.logical(d$supporting_specimens>=minimum_specimens))
  list(primary=choice(primary),alternative=choice(alternative),coverage=coverage,
       identity_coverage=do.call(rbind,identity_rows),
       same_requested_State=identical(primary$requested_State,alternative$requested_State),
       policy='Most specimen-supported graph endpoint; alternative maximizes tree distance among supported distinct-State endpoints when available. Identity only audited, never ranked. Native root_state/fallback determines actual root, always reported.')
}
