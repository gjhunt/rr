source('plot_fns/cca.R')

## LEFT SINGULAR VECTORS

btchs = c("well","stain","plate","ligand")
# SVD FULL
ccas = lapply(fctrs[btchs],function(batch)make_cca(batch,SVDs=SVD_no_center))

plts = lapply(btchs,
              function(batch){
                  print(batch)
                  b = fctrs[[batch]]
                  save_cca_plots(ccas[[batch]]$cca,NULL,batch=b,SVDs=SVD_no_center,uid=gsub(" ","_",b$name),subdir='all_svd/')
              })

# SVD (only examples)
ex_SVDs = SVD_no_center[ex_ss]
ccas_ex = lapply(ccas,function(x)x[[1]][ex_ss])

ex_plts = lapply(btchs,
              function(batch){
                  print(batch)
                  b = fctrs[[batch]]
                  save_cca_plots(ccas_ex[[batch]],NULL,batch=b,SVDs=ex_SVDs,uid=gsub(" ","_",b$name),subdir='svd/')
              })
names(ex_plts) = btchs

lapply(btchs,function(batch)ggsave(paste0(plotdir,"cca/svd/",batch,"_meta_cca.pdf"),meta_cca_plot(ex_plts[[batch]]$data,ex_SVDs,batch),height=4,width=10))

# ASVD
asvd_SVD_list = list("ASVD"=DAPIASVD_no_center[[1]],"HL ASVD"=DAPIASVD_hl[[1]])
asvd_ccas = lapply(fctrs[btchs],function(batch)make_cca(batch,SVDs=asvd_SVD_list))

asvd_plts = lapply(btchs,
              function(batch){
                  print(batch)
                  b = fctrs[[batch]]
                  save_cca_plots(asvd_ccas[[batch]]$cca,NULL,batch=b,SVDs=asvd_SVD_list,uid=gsub(" ","_",b$name),subdir='asvd/')
              })
names(asvd_plts) = btchs

lapply(btchs,function(batch)ggsave(paste0(plotdir,"cca/asvd/",batch,"_meta_cca.pdf"),meta_cca_plot(asvd_plts[[batch]]$data,asvd_SVD_list,batch),height=4,width=10))

## RIGHT SINGULAR VECTORS

# Example Features
#btchs = as.vector(fctrs[["ecmp"]]$unique)
#ex_plts = lapply(fctrs[btchs],
#                 function(batch)save_cca_plots(batch,SVDs=SVD_no_center,uid=gsub(" ","_",batch$name),subdir='right_svd_examples/',chir="v",kmax=2))
