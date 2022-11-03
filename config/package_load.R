# Check if packages are installed
listOfAllPackages = c("Biobase","dplyr","QDNAseq","QDNAseqmod","GenomicRanges")

for(thisPackage in listOfAllPackages) {
  
  if(thisPackage %in% rownames(installed.packages()) == FALSE) {
	cat(paste("[install_env] Package", thisPackage, "needs installing.\n"))
	install.packages(thisPackage, repos="https://cran.ma.imperial.ac.uk/")
  } else {
	cat(paste("[install_env] Package", thisPackage, "is installed.\n"))
  }
}
a <- lapply(listOfAllPackages,FUN = function(x){suppressPackageStartupMessages(require(x, character.only = TRUE))})

if(!is.null(warnings())){
	print(warnings())
}
