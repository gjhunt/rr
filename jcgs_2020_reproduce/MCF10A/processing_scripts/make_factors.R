make_factors=function(file_in,file_out=""){
    ww = readRDS(file_in)
    ECMp_mtx = as.matrix(ww$ECMp)
    uECMp = unique(as.vector(ECMp_mtx))
    uECMp <- uECMp[!is.na(uECMp)]
    for(ecmp in uECMp){
        tmp <- ECMp_mtx
        tmp[tmp!=ecmp] <- "Other"
        ww[[ecmp]] <- tmp
    }

    make_factors = function(nme,factor_across=1){
        mtx = ww[[nme]]
        mtx = t(apply(mtx,factor_across,function(x)rep(unique(x[!is.na(x)]),length(x))))
        if(factor_across==2)
            mtx = t(mtx)
        nme=tolower(nme)
        factor = factor(apply(mtx,factor_across,function(x)unique(x[!is.na(x)])),ordered=TRUE)
        uni = unique(factor)
        basis = model.matrix(~0+factor)
        return(list(name=nme,mtx=mtx,factor=factor,unique=uni,basis=basis))
    }

    well_f = make_factors("Well")
    ligand_f = make_factors("Ligand")
    stain_f = make_factors("StainingSet")
    stain_f$name = "Staining Set"
    plate_f = make_factors("Barcode")
    plate_f$name = "Plate"

    row_f = make_factors("PrintRow",factor_across=2)
    col_f = make_factors("PrintColumn",factor_across=2)

    ecmp_f = make_factors("ECMp",factor_across=2)

    fctrs = list(well=well_f,ligand=ligand_f,stain=stain_f,plate=plate_f,ecmp=ecmp_f,col=col_f,row=row_f)

    for(ecmp in uECMp){
        fctrs[[ecmp]] <- make_factors(ecmp,factor_across=2)
    }
    
    saveRDS(fctrs,paste0("../data/fctrs",file_out,".rds"))
}

make_factors("../data/all_no_nid_eln_thbs_fn","_no_nid_eln_thbs_fn")
make_factors("../data/all_no_nid_eln","_no_nid_eln")
make_factors("../data/all")
