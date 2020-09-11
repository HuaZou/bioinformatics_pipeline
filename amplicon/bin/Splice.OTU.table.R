args <- commandArgs(T)
if(length(args)!=3){
    stop(
"
Usage:
\tRscript Splice.OTU.table.R [OTU] [num column] [out dir]
\t[OTU] : OTU table with header and rowname
\t[num columns] : number of columns to 1 file
\t[out dir] : out directory
")
}

otu <- read.table(args[1],row.names=1,header=T,sep="\t",check.names=F)
num <- as.numeric(args[2])
out_dir <- args[3]
n <- 0
for(i in seq(1,ncol(otu),num)){
    n_end <- i+num-1
    if(n_end > ncol(otu)){
        n_end <- ncol(otu)
    }
    aa <- otu[,i:n_end,drop=F]
    n <- n+1
    write.table(aa,paste0(out_dir,"/",n,".txt"),row.names=T,col.names=NA,sep="\t",quote=F)
}

