library('rrscale')
library('foreach')
library('doParallel')
library('R.utils')

args = commandArgs(asValues=TRUE)

scale_data = function(data_tag=NULL){

    if(is.null(data_tag))
        data_tag <- ""
    
    cat(paste(Sys.time()),"\n")
    cat("Read data.\n")
    Ys = readRDS(paste0("../data/Ys",data_tag))

    outdir = "./out/"
    dir.create(outdir,showWarnings=FALSE)

    rrscale_i <- function (i) {
        capture.output({        
            cat(i,"\n")
            Y = Ys[[i]]
            ret=tryCatch({
                rrscale::rrscale(Y=Y,trans_list = list(bcn=rrscale::box_cox_negative, asinh=rrscale::asinh),
                                     lims_list=list(c(-100,100),c(0,100)),run_parallel=FALSE,
                                     verbose=FALSE)
            },error=function(e){
                print(e)
                print("Continuing...")
                return(NULL)
            })
        },file=paste0(outdir,"out",i),split=TRUE)
        return(ret)
    }

    cat(paste(Sys.time()),"\n")
    cat("Scaling.\n")

    no_cores <- floor(detectCores()/2)
    cl<-makeCluster(no_cores)
    registerDoParallel(cl)
    scaledY = foreach(i=1:length(Ys)) %dopar% rrscale_i(i)
    stopCluster(cl)

    names(scaledY) = names(Ys)

    scaledY = scaledY[lengths(scaledY)>0]

    dir.create("../scaledY",showWarnings=FALSE)
    saveRDS(scaledY,paste0("../scaledY/scaledY",data_tag,".rds"))

    cat(paste(Sys.time()),"\n")
    cat("Scaled.\n")
}

print(args$tag)
scale_data(data_tag=args$tag)
