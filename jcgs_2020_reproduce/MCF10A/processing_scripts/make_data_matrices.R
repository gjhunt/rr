library('dplyr')
library('data.table')
library('tidyr')

# Read in data
data_dir = "../mcf10a_data/"
files = dir(data_dir, full.names=TRUE)
read_files = lapply(files,read.csv,sep="\t",header=TRUE)

# Bind together rows of plates
d = bind_rows(read_files)
d$PrintRow = floor((d$PrintSpot-1)/20)+1
d$PrintColumn = (d$PrintSpot-1)%%20+1

#Remove punctuation from ECMP names
d$ECMp = gsub('[[:punct:]]+','_',d$ECMp)

feature_cols = grep("PA|CP",colnames(d),value=TRUE)
feature_cols = feature_cols[!grepl("_SE",feature_cols)]

exclude_features = c("Spot_PA_Perimeter",
                     "Nuclei_PA_Cycle_State",
                     "Cytoplasm_PA_Gated_KRTClass",
                     "Spot_PA_ReplicateCount",
                     "Spot_PA_LoessSCC",
                     "Nuclei_PA_Gated_EdUPositive")

feature_cols = setdiff(feature_cols,exclude_features)

meta_cols = colnames(d)
meta_cols = meta_cols[!grepl("PA|CP|_SE|Row|Col",meta_cols)]

for(j in meta_cols){
    d[,j] <- as.factor(d[,j])
}
d = data.table(d,stringsAsFactors=FALSE)

save_data = function(to_rmv,file_name=""){
    d = d[!grepl(to_rmv,d$ECMp),]

    # Make well by spot matrices
    make_v = function(v){
        sd = d[,.(Barcode,PrintSpot,Well,v=get(v))]
        sd = sd %>% spread(PrintSpot,"v")
        rns = apply(sd[,1:2],1,function(x)paste(x,collapse="_"))
        sd = sd[,-(1:2)]
        rownames(sd) = rns
        return(sd)
    }

    v_names = colnames(d)
    d_list = lapply(v_names,make_v)
    names(d_list) = v_names

    # Remove metadata covariates and only includes cell features
    Ys = d_list[feature_cols]

    dir.create('../data/',showWarnings=FALSE)
    saveRDS(d_list,paste0('../data/all',file_name))
    saveRDS(Ys,paste0('../data/Ys',file_name))
}

save_data("fiducial|Fiducial|gelatin|blank|air|PBS|ELN|NID|THBS|FN",'_no_nid_eln_thbs_fn')
save_data("fiducial|Fiducial|gelatin|blank|air|PBS|ELN|NID",'_no_nid_eln')
save_data("fiducial|Fiducial|gelatin|blank|air|PBS")
