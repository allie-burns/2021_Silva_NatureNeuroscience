#################################################################################
########################  Binning Raw Photometry data  ##########################
##
## Align photometry times and average fiber raw photometry readings into 10ms
## time bins to reduce file sizes for next steps in Igor.
##
## Input (all files to process should be in same directory set by path)
##  - Raw Photometry Input (/data/RawPhotometry_Input/)
##    - Output from 1-site 2-color Fiber Photometry System (Doric Lenses, Canada)
##    - Columns include Time(s), 405nm and 465nM signal readings
##    - One file per mouse per phase (recall, extinction, SR, etc)
##    - Currently no rule for file names
##
## Output:
## - New directory (/data/RawPhotometry_Input/binned_output/) with all new binned data
##   - Includes one output file per input file
##   - File names: InputFileName_binned.csv
##   - File includes:
##     - Times (s)
##       - exclude negative signal readings (before wire connection to mouse)
##       - averaged into 10ms time bins
##
################################################################################

#########################################################
## STEP 1 : set parameters
#########################################################
library(readr)
library(dplyr)

## Get file list from Folder
path <- "../data/RawPhotometry_Input"
files <- list.files(path,full.names=TRUE)

## Set up file for saving (saving will occur wherever you run the script)
if(!dir.exists(paste("../data/RawPhotometry_Input/binned_output",sep=""))) { 
    dir.create(paste("../data/RawPhotometry_Input/binned_output",sep=""))
}

######################################################
## STEP 2: Run the binning script
######################################################
BinningByFile <- function(i,files) {
    ## Select one file
    filePath <- files[[i]]
    cat(paste("Load and process raw file:", basename(filePath),"\n",sep=" "))
    
    ## Load File and remove extra column 
    x <- read_csv(filePath)
    x <- x[,1:6]
    
    ## Remove rows where Aln-3 < 0 for aligning times and fix for time 0
    FC <- x[x$'AIn-3' > 0,]
    FC$`Time(s)`<- FC$`Time(s)`-min(FC$`Time(s)`)
    
    ## Remove duplicated rows
    FC <- FC[!duplicated(FC),]
    
    ## Set bin regions by the 0.01 (s)
    binsize <- 0.01
    lbin <- seq(min(FC$`Time(s)`),max(FC$`Time(s)`),by=binsize)
    hbin <- seq(binsize,max(FC$`Time(s)`)+binsize,by=binsize)

    ## Bin
    cat(paste("Binning:", basename(filePath),"\n",sep=" "))
    binSecs <- function(i,lbin,hbin,FC){
        low <- lbin[i]
        high <- hbin[i]
        
        sVals <- FC[FC$`Time(s)` >= low & FC$`Time(s)`< high,]  
        
        averages <- tbl_df(do.call(cbind,lapply(sVals,mean)))
        averages$`Time(s)`<- low
        averages
    }

    binFC <- lapply(seq(1:length(lbin)),binSecs,lbin,hbin,FC)
    binFC <- do.call(rbind, binFC)
    
    ## Output CSV file
    fl <- substr(basename(filePath),1,nchar(basename(filePath)) -4)
    fl <- paste(fl,"binned.csv",sep="_")
    write.csv(binFC,
              file=paste("../data/RawPhotometry_Input/binned_output/",fl,sep=""),
              sep=',',
              row.names=FALSE)
    
}


lapply(seq(1:length(files)),BinningByFile,files)
