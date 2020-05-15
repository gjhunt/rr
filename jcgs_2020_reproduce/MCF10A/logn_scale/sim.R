library('reshape2')
library('ggplot2')
library('grid')
library('gridExtra')
library('cowplot')
library('viridis')
library('RColorBrewer')

theme_set(theme_bw())

set.seed(3251991)

make_Y = function(){
  u = rlnorm(100)
  v = rlnorm(100)
  e = rnorm(100^2,0,1)
  Y = u%*%t(v) + abs(array(e,c(100,100)))
  return(Y)
}

fctr = c(rep(0,100),rep(1,100))

Y = rbind(make_Y(),1.5+make_Y())
logY = log(Y)
#Y = scale(Y,scale=FALSE)
#logY = scale(logY,scale=FALSE)
#Y = scale(Y,scale=FALSE)
#logY = scale(logY,scale=FALSE)

sigY = svd(Y,nu=200,nv=100)$d^2
sigY = data.frame(index=1:100,"Cumul. Pct. of Variance"=cumsum(sigY/sum(sigY)),"Transformation"="Untransformed")
siglogY = svd(logY,nu=200,nv=100)$d^2
siglogY = data.frame(index=1:100,"Cumul. Pct. of Variance"=cumsum(siglogY/sum(siglogY)),"Transformation"="log")
sd = rbind(sigY,siglogY)
colnames(sd) = c("Index","Cumul. Pct. of Variance","Transformation")
pct_plt = ggplot(data=sd,mapping=aes(x=Index,y=`Cumul. Pct. of Variance`,group=Transformation,lty=Transformation))+geom_hline(yintercept=1,lty=2)+geom_line(lwd=2)+ggtitle("Recovery of Total Variance by PCs")+labs(x="Num. PCs (k)")+ylim(0,1)+xlim(1,5)+scale_linetype_manual(values = c(1,3))
ggsave("pct.pdf",pct_plt)

uY = svd(Y,nu=200,nv=100)$u
ulogY = svd(logY,nu=200,nv=100)$u
vY = svd(Y,nu=200,nv=100)$v
vlogY = svd(logY,nu=200,nv=100)$v

d=melt(list("Untransformed"=Y,"log"=logY))
d$f = as.factor(fctr)
levels(d$f)=c("Group 1", "Group 2")
colnames(d) = c("row","col","value","trans","group")
d$trans = as.factor(d$trans)
d$trans = ordered(d$trans,levels=c("Untransformed","log"))
d$value[d$value>50] <- NA
d_plts = ggplot(data=d,mapping=aes(x=value))+ggtitle("Density of Untransformed and log-Transformed Data")+labs(x="Value",y="Density")+
    theme(legend.position="right",legend.text=element_text(size=20),legend.title=element_text(size=20))+
    geom_density(data=subset(d,group=="Group 1"),mapping=aes(fill=group,lty=group),size=1,alpha=1)+
    geom_density(data=subset(d,group=="Group 2"),mapping=aes(fill=group,lty=group),size=1,alpha=.75)+
    facet_wrap(~trans,scales="free") +
    geom_vline(data=aggregate(value~group+trans,data=d,median),mapping=aes(xintercept=value,linetype=group),size=1)+scale_fill_manual(values=c("#999999", "#56B4E9"))
ggsave("dens.pdf",d_plts,height=4,width=10)

d = rbind(data.frame(pc1=Y%*%vY[,1],pc2=Y%*%vY[,2],trans="Untransformed",Group=as.factor(fctr)),
          data.frame(pc1=logY%*%vlogY[,1],pc2=logY%*%vlogY[,2],trans="log",Group=as.factor(fctr)))
pcs_plt = ggplot(data=d,mapping=aes(x=pc1,y=pc2,lty=Group))+geom_point()+facet_wrap(~trans,scales="free")+ggtitle("Density of Untransformed and log-Transformed Data")+labs(x="PC1",y="PC2")+scale_colour_brewer(palette = "Accent")
ggsave("pcs.pdf",pcs_plt,height=5,width=10)

ccY = sapply(1:200,function(k)cancor(fctr,uY[,1:k,drop=FALSE])$cor^2)
cclogY = sapply(1:200,function(k)cancor(fctr,ulogY[,1:k,drop=FALSE])$cor^2)

d = rbind(data.frame(index=1:200,cc=ccY,Transformation="Untransformed"),
          data.frame(index=1:200,cc=cclogY,Transformation="log"))
cca_plt = ggplot(data=d,mapping=aes(x=index,y=cc,lty=Transformation,group=Transformation))+geom_hline(yintercept=1,lty=3)+geom_line(lwd=2)+labs(x="Num. PCs (k)",y="Mean Squared Can. Cor. \n between Batch and First PCs")+ggtitle("Recovery of Group by PCs")+ylim(0,1)+xlim(1,5)+scale_linetype_manual(values = c(1,3))
ggsave("cca.pdf",cca_plt)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

dual_plts = plot_grid(pct_plt+guides(lty=FALSE),cca_plt+guides(lty=FALSE),labels=c("A","B"))

glist = list(ggplotGrob(dual_plts),g_legend(cca_plt+theme(legend.position="bottom")))
lmx = matrix(c(rep(1,3*10),NA,2,NA),nrow=11,byrow=TRUE)
gbs = arrangeGrob(grobs=glist,layout_matrix=lmx)

pdf("pca.pdf",height=4,width=10)
grid.draw(gbs)
dev.off()
