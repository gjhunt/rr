library('ggplot2')

save_plots=function(pl,dir,subdir=NULL,name,uid=NULL,
                    height=10,width=15){
    graphics.off()
    dir = paste0(dir,subdir)
    dir.create(dir,showWarnings=FALSE,recursive=TRUE)
    fname = paste0(dir,name,uid,".pdf")
    cat(paste0("Saving ",fname,"\n"))
    pdf(fname,width=width,height=height)
    if("ggplot"%in%class(pl)){
        print(pl)
    } else {
        for(i in seq_along(pl)){
            tryCatch({
                p = pl[[i]]
                grid.draw(p)
                if(("gtable" %in% class(p)) && (i != length(pl)))
                    grid.newpage()
            },error=function(e){
                print(e)
            })
        }
    }
    dev.off()    
}

proj = function(A){
    return(A%*%MASS::ginv(base::t(A)%*%A)%*%base::t(A))
}

resid = function(A){
    P = proj(A)
    return(diag(rep(1,dim(P)[1]))-P)
}

center=function(x)scale(x,scale=FALSE)
centert=function(x)t(scale(t(x),scale=FALSE))

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

meta_title = function(p,title,size=14){
    title <- ggdraw() + draw_label(title, fontface='bold',size=size)
    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
}

ex_arrange=function(ex_plts){
    xlab = ex_plts[[1]]$labels$x
    ylab = ex_plts[[1]]$labels$y
    leg = g_legend(ex_plts[[1]]+theme(legend.position="bottom"))
    ex_plts = lapply(seq_along(ex_plts),function(i)ex_plts[[i]]+labs(title=paste0("(",LETTERS[i],") ",names(ex_plts)[i])))
    ex_plts = lapply(ex_plts,function(x)x+guides(color=FALSE)+labs(x=NULL,y=NULL))
    lout = t(array(c(rep(1,20),2),c(1,21)))
    ag = arrangeGrob(grobs=ex_plts,ncol=2)
    ag = arrangeGrob(grobs=list(ag,textGrob(xlab)),
                     layout_matrix=lout)
    ag = arrangeGrob(grobs=list(ag,textGrob(ylab,rot=90)),
                     layout_matrix=t(rev(lout)))
    leglout = array(c(rep(1,15),rep(2,3)),c(18,1))
    ag = arrangeGrob(grobs=list(ag,leg),
                     layout_matrix=leglout)
    return(ag)
}

nostrip = theme_dark()+theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                strip.text.y = element_blank(),
                panel.spacing.x=unit(0,"lines"),
                panel.spacing.y=unit(0,"lines"),
                panel.spacing=unit(0,"lines"),
                plot.margin=unit(rep(0,4),"cm"))

ggplotColors <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
