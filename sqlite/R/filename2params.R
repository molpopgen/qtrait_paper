library(stringr)

getparams=function(f)
{
    s=str_split(f,'_')
    return(list(opt=as.double(s[[1]][4]),mu=as.double(s[[1]][6])))
}


