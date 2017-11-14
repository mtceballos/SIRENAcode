library(FITSio)

readKeyword <- function(fitsFileName,keyword){
  #fitsFileName
  #keyword
  
  zz <- file(description = fitsFileName, open = "rb")
  header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
  header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
  close(zz)
  
  rowWhereKeywwordIs <- header[grep(keyword,header)]
  start <- regexpr("=",rowWhereKeywwordIs)[[1]]+1
  end <- regexpr("/",rowWhereKeywwordIs)[[1]]-1
  
  return(as.numeric(substring(rowWhereKeywwordIs,start,end)))
  
}