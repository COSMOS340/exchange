a=commandArgs(TRUE);d=read.csv(a[1]);p=rep(NA_real_,nrow(d));ok=d$tested & is.finite(d$F_statistic)
p[ok]=pf(d$F_statistic[ok],5,as.integer(a[3]),lower.tail=FALSE)
q=p;q[ok]=p.adjust(p[ok],method='BH')
write.csv(data.frame(p_value=p,BH_q_value=q),a[2],row.names=FALSE)
