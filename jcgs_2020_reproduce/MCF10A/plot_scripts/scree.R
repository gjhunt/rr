RRsvds = lapply(SVD_no_center,"[[","RR")
RRdu = lapply(RRsvds,"[[","du")
RRdu = lapply(RRdu,function(x)x^2/sum(x^2))
RRdu = data.frame(RRdu)
RRdu$K <- 1:192
mRRdu = data.table(melt(RRdu,id.vars='K'))
colnames(mRRdu) <- c("K","Feature","Pct")
mRRdu$Trans <- "RR"

NTsvds = lapply(SVD_no_center,"[[","NT")
NTdu = lapply(NTsvds,"[[","du")
NTdu = lapply(NTdu,function(x)x^2/sum(x^2))
NTdu = data.frame(NTdu)
NTdu$K <- 1:192
mNTdu = data.table(melt(NTdu,id.vars='K'))
colnames(mNTdu) <- c("K","Feature","Pct")
mNTdu$Trans <- "NT"

du = rbind(mRRdu,mNTdu)


plt = ggplot(data=du[Feature%in%ex_ss,],mapping=aes(x=K,y=Pct,color=Feature))+geom_point()+scale_y_log10()+labs(x="Num. Singular Vectors",y="Pct. of Variance (log-scale)")+ggtitle("Pct of Variance Captured by Singular Vectors")+facet_wrap(~Trans)
ggsave(paste0(plotdir,"scree.pdf"),plt,width=15,height=12)
