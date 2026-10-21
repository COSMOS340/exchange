# Synthetic software checks only. Does not fit or subset biological data.
args = commandArgs(TRUE)
stopifnot(length(args) == 1)
source_text = paste(readLines(args[1], warn=FALSE), collapse="\n")
begin = regexpr('patch_monocle2_dplyr <- function()', source_text, fixed=TRUE)[1]
finish = regexpr('options(warn=1)', source_text, fixed=TRUE)[1]
stopifnot(begin > 0, finish > begin)
eval(parse(text=substr(source_text, begin, finish-1)))
checks = list()
set.seed(721)
for (case in c("distinct", "duplicates", "grid")) {
  z = switch(case, distinct=matrix(rnorm(80), ncol=2),
             duplicates=rbind(matrix(rnorm(76), ncol=2), c(0,0), c(0,0)),
             grid=as.matrix(expand.grid(1:5,1:5)))
  rownames(z) = paste0("v", seq_len(nrow(z)))
  for (shifted in c(FALSE, TRUE)) {
    d = as.matrix(dist(z))
    if (shifted) d = d + min(d[d>0])
    diag(d) = 0
    dense = mst(graph_from_adjacency_matrix(d, mode="undirected", weighted=TRUE, diag=FALSE))
    small = complete_tree(z, shifted=shifted, skip_zero=!shifted)
    stopifnot(is_tree(small), vcount(small)==nrow(z), ecount(small)==nrow(z)-1,
              abs(sum(E(dense)$weight)-sum(E(small)$weight))<1e-8)
    tr = iterative_dfs(small, V(small)$name[1])
    native = dfs(small, root=1, parent=TRUE, unreachable=FALSE)
    stopifnot(identical(as.integer(tr$order), as.integer(native$order)),
              identical(as.integer(tr$parent), as.integer(native$parent)))
    checks[[paste(case,shifted)]] = TRUE
  }
}
patch_monocle2_dplyr(); patch_monocle2_igraph(); memory_patch(); patch_root_fallback()
stopifnot(grepl("ROOT_FALLBACK", paste(deparse(get("orderCells", asNamespace("monocle"))),collapse="")))
jsonlite::write_json(list(passed=TRUE, tests=checks,
  compatibility_installation=TRUE, root_fallback_text_installed=TRUE,
  real_model_root_distance_check_pending=TRUE, formal_biological_fit_performed=FALSE),
  stdout(), pretty=TRUE, auto_unbox=TRUE)
