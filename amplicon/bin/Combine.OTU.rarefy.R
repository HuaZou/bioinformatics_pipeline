args <- commandArgs(T)

Usage <- function(){
    stop(
"Usage:\n\tRscrpts [scripts] [file type] [file] [out]
Description:
\tcombine OTU tables by a given directory or OTU file list
\t[file type]: L for OTU list, D for directory 
\t[out]      : name of output file
Examples:
\tRscripts combine.OTU.R L xx.list.txt OTU.combined.txt
\tRscripts combine.OTU.R D xxxx/ OTU.combined.txt"
        )
}

# args validation
n_args = 3
if (length(args) != n_args){Usage()}

# args parser
file_type = args[1]
dat_input = args[2]
output = args[3]

# data lists
if (file_type == "L"){
    print ("debug")
    dat_lst = read.table(dat_input,header=F)
    dat_lst = as.character(dat_lst[,1])
}else if (file_type == "D"){
    dat_lst = list.files(dat_input,recursive =T, full.names=T)
    print (dat_lst)
}else{
    Usage()
}

# data reading
ReadFiles <- function(file_path){
    tmp <- read.table(file_path,header=T,sep="\t")
    colnames(tmp)[1] <- "OTU_ID"
    #tmp <- tmp[,-ncol(tmp)]
    
    #cn <- unlist(strsplit(file_path,"/"))
    #cn <- cn[length(cn)-1]
    #colnames(tmp) <- c("OTU_ID",cn)
    return(tmp)
}

dat <- ReadFiles(dat_lst[1])
for (i in 2:length(dat_lst)){
    tmp <- ReadFiles(dat_lst[i])
    dat <- merge(dat,tmp,by="OTU_ID",all=T)
}

dat[is.na(dat)] <- 0
# output
write.table(dat,output,row.names=F,sep="\t",quote=F)
