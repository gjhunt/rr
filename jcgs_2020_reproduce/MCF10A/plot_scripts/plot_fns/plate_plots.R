save_plate_plots = function(scaledY,uid=""){

    make_plate_plts = function(sdn,nm=""){
        cat(nm,"\n")
        row_mtx = as.numeric(fctrs$row$mtx)
        col_mtx = as.numeric(fctrs$col$mtx)
        plate_mtx = fctrs$plate$mtx
        well_mtx = fctrs$well$mtx

        trans = TRANSFORMATIONS
        spotwise = function(trans){
            cat(trans,"\n")
            Ytd = melt(as.matrix(sdn[[trans]]))
            Ytd$row=melt(row_mtx)$value
            Ytd$col=melt(col_mtx)$value
            Ytd$plate=melt(plate_mtx)$value
            Ytd$well=melt(well_mtx)$value
            Ytd$trans = trans
            return(Ytd[,-(1:2)])
        }
        d = Reduce("rbind",lapply(trans,spotwise))
        d = data.table(d)
        d$well_letter = factor(gsub("[[:digit:]]","",d$well))
        d$well_number = factor(gsub("[^[:digit:]]","",d$well))

        brks_mtx = d[,.(min = min(value,na.rm=TRUE),
                        max= max(value,na.rm=TRUE),
                        mdn=min(value,na.rm=TRUE)+.5*diff(range(value,na.rm=TRUE))),
                     by=.(trans)]

        lay = Reduce("rbind",
            rep(list(c(rep(1,,16),2)),2))

        wells_heat = function(plte,transform,
                              wells=c("A01","A02","A03","A04","B01","B02","B03","B04"),
                              legend=FALSE,
                              title=FALSE){
            dta = d[plate==plte & trans==transform & well %in% wells,]
            dta = data.table(dta)

            if(length(wells)>1){
                totally_empty = dta[,.(is_na=all(is.na(value))),by=.(well)]
                totally_empty = totally_empty$well[totally_empty$is_na]
                dta = dta[!(well %in% totally_empty),]
                if(nrow(dta)==0)
                    return(NULL)
            }
            
            p = ggplot(data=dta,
                       mapping=aes(x=col,y=row,fill=value))
            p = p + geom_tile() + geom_tile(data = subset(dta, !is.na(value)), aes(fill = value)) +
                #geom_tile(data = subset(dta,  is.na(value)), fill="green") +
                facet_grid(well_letter~well_number)

            
            if(title)
                p = p+ggtitle(plte)
            if(!legend)
                p = p + guides(fill=FALSE)
            p = p+nostrip
            brks = unlist(brks_mtx[trans==transform,c("min","mdn","max")])
            #brks = sapply(brks,round,5)
            p = p+scale_fill_gradient2(
                      mid='white',
                      na.value = "green",
                      midpoint=brks[2],
                      breaks=brks,
                      labels=round(brks,5),limits=range(brks)
                  ) +
                theme(panel.spacing = unit(.5, "lines"))
            p = p + scale_x_continuous(expand=c(0,0)) + 
                scale_y_reverse(expand=c(0,0))
            if(length(wells)>1)
                p = p + theme(plot.margin=unit(rep(.05,4),"cm"))

            return(p)
        }

        by_plate_plts = function(trans){
            cat("plotting...",trans,"\n")
            ws = c("A01","A02","A03","A04","B01","B02","B03","B04")
            wells_Yt_byw=lapply(unique(d$plate),function(plate){
                L = lapply(ws,function(w)wells_heat(plate,trans,w,title=FALSE))
                names(L) = ws
                return(L)
            })
            names(wells_Yt_byw) = unique(d$plate)
            wells_Yt = plot_grid(plotlist=lapply(unique(d$plate),wells_heat,trans),nrow=3)
            wells_Yt = meta_title(wells_Yt,paste0(trans),size=35)
            leg = g_legend(wells_heat(unlist(d$plate[1]),trans,legend=TRUE))
            grd = arrangeGrob(grobs=list(wells_Yt,leg),layout_matrix=lay)
            return(list(plot=grd,wells=wells_Yt_byw))
        }

        pps_byw = lapply(TRANSFORMATIONS,by_plate_plts)
        names(pps_byw) = TRANSFORMATIONS
        pps = lapply(pps_byw,"[[","plot")
        names(pps) = TRANSFORMATIONS
        names(pps) = gsub("[[:punct:]]| ","_",tolower(TRANSFORMATIONS))
        save_plots(pps,dir=plotdir,subdir="obj_dens/plates_trans/",name=nm,uid=uid,
                   height=12,width=25)
        return(list(plot=pps,plates=pps_byw))
    }
    
    plts_byp = lapply(seq_along(scaledY),function(i)make_plate_plts(scaledY[[i]],names(scaledY)[i]))
    names(plts_byp) = names(scaledY)
    plts = lapply(plts_byp,"[[","plot")
    names(plts) = names(scaledY)
    names(plts) = names(scaledY)
    #for(tns in names(plts[[1]])){
    #   sub_plts = lapply(plts,"[[",tns)
    #    save_plots(sub_plts,dir=plotdir,subdir="obj_dens/",name=paste0("plate_",tns),uid=uid,
    #               height=12,width=25)
    #}
    return(list(plot=plts,list=plts_byp))
}
