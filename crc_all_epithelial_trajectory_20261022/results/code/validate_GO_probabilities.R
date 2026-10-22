a=commandArgs(TRUE);d=read.csv(a[1]);ok=d$background_hits>0
p=rep(NA_real_,nrow(d));pa=p
p[ok]=phyper(d$module_hits[ok]-1,d$background_hits[ok],d$background_size[ok]-d$background_hits[ok],d$module_size[ok],lower.tail=FALSE)
pa[ok]=phyper(d$module_hits[ok]-1,d$background_hits[ok],d$BP_annotated_background_size[ok]-d$background_hits[ok],d$BP_annotated_module_size[ok],lower.tail=FALSE)
q=qa=p;q[is.finite(p)]=p.adjust(p[is.finite(p)],method='BH');qa[is.finite(pa)]=p.adjust(pa[is.finite(pa)],method='BH')
write.csv(data.frame(p_value=p,BH_q_value=q,p_value_annotated_background=pa,BH_q_annotated_background=qa),a[2],row.names=FALSE)
