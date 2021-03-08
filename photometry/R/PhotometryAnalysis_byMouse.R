################################################################################
########################  Photometry-Freezing Analysis  ########################
##
## Calculate the photometry patterns during freezing for each animal in each phase
##
## Input (should all be in directory set below):
##  - Photometry Files after Igor processing (IgorOutput)
##     - One file per mouse per phase (recall, extinction, SR, etc)
##     - previously binned using BinTimePoints.R and processed using Igor pipeline
##     - File contains time_ext, df_405, df_465 columns
##     - File name should be:
##        - mouse_phase#_phaseName_binned_binningvalue_dF.txt
##        - EX: 10005_phase6_ext4_6_binned_0.01_dF.txt
##  - Freezing file (TSE_freezingOutput)
##     - One file represents many mice and many phases
##     - File output from TSE software
##     - File converted to tab delimited Txt on a windows computer(using TSE software)
##     - File name should be:
##        - freezing_photoMouse#_listofPhaseNames_freezingcutoff.txt
##        - EX: freezing_photo95540_Ext3_Ext4_Ext5_Ext6_Ext7_cutoff0.5.txt
##
## Output:
##  - CSV file with photometry readings for every freezing window by mouse, by phase
##
################################################################################

#########################################################
## Set parameters
#########################################################
## Loading Libraries
library(plotly)
library(plotrix)
library(RColorBrewer)
library(wesanderson)

directory <- "../data/Photometry_Output/"
dir.create(paste(directory,"byMouse/",sep=''))
setwd(directory)

## Build Comparison Info
fls<- list.files("../Photometry_Input/Igor_output/",pattern="^[0-9]")
comp <- data.frame(
    PHASE_NAME = unlist(tolower(lapply(fls, function(x) {strsplit(x,split="_")[[1]][[3]][[1]]}))),
    ANIMAL = unlist(lapply(fls, function(x) {strsplit(x,split="_")[[1]][[1]][[1]]})),
    PHASE_NUMBER = unlist(tolower(lapply(fls, function(x) {strsplit(x,split="_")[[1]][[2]][[1]]}))))

phases <- unique(comp$PHASE_NAME)

## Set up for saving
date <- gsub('-','',Sys.Date())

## Set Analysis parameters
tcutoff   <- 60   ## Remove first 60 seconds of experiment
minfreeze <- 1.5  ## set minimum freezing length 
maxfreeze <- Inf  ## set maximum freezing length
window    <- 2    ## set  window around end of freezing

###########################################################
## Setup Mouse averaging Function
###########################################################
AverageMouse <- function(i,phases,comp){
    myphase <- comp[grep(phases[[i]],comp$PHASE_NAME,ignore.case=TRUE),]
    cat(paste("Analyzing average freezing behavior for", phases[[i]],"\n",sep =' '))
    
    traceData <-  function(j,myphase) {
        ## ###############################################
        ## Set Data parameters and load a process
        ## ###############################################
        phasename <- myphase$PHASE_NAME[j]
        animal <- myphase$ANIMAL[j]
        phase <- tolower(myphase$PHASE_NUMBER[j])
        ph <- gsub('[a-z]+',"",phase)
        cat(paste(animal, "-", phasename, "-", ph, "\n", sep=" "))
        
        ## Load and process Freezing Data
        freezeTable <-list.files("../Photometry_Input/TSE_freezingOutput",
                                 pattern="freez", full.names = TRUE)
        freezeTable <- freezeTable[grep(phasename,freezeTable,ignore.case=TRUE)]

        rawFreeze <- do.call(rbind,
                             lapply(freezeTable, function(freezeTable) {
                                 read.table(freezeTable,skip=3,dec=",")
                             }))
        columns <- read.table(freezeTable[1],
                              nrow=1,skip=1,sep="\t",
                              stringsAsFactors=FALSE,
                              header=T,check.names=T)
        
        colnames(rawFreeze) <- colnames(columns)
        
        freeze <- rawFreeze[rawFreeze$Animal == animal,] ##filter by animal
        freeze <- freeze[freeze$Phase == ph,] ## filter for phase
        freeze <- freeze[freeze$Time.from..s. >= tcutoff,] ## Remove first minute of data
        freeze_cut <- freeze[freeze$Dur...s. >= minfreeze & freeze$Dur...s. < maxfreeze,]
        
        if(nrow(freeze_cut) == 0){
            ##return ("No freezing fits criteria")
            return(data.frame(animal=animal,
                              phasename=phasename,
                              seconds=2.1,
                              average=0))
        }
        
        freezeRange <- data.frame(pre2 = freeze_cut$Time.to..s. - window,
                                  post2 = freeze_cut$Time.to..s. + window)
        
        ## Load and process photometry data
        photometryTable <- list.files("../Photometry_Input/Igor_output",
                                      pattern="binned", full.names = TRUE)
        photometryTable <- photometryTable[grep(phasename,photometryTable, ignore.case=TRUE)]

        photometryTable <- photometryTable[grep(animal,photometryTable)]
        photometryTable <- photometryTable[grep(phase,photometryTable)]
        
        rawPhot <- read.table(photometryTable,sep='\t',header=TRUE)
        
        getWindows <- function(i){
            start <- freezeRange[i,]$pre2
            end <- freezeRange[i,]$post2
            rawPhot[rawPhot$time_ext>= start & rawPhot$time_ext<= end,]$dF_465
        }
        
        photWindows <- lapply(seq(1:nrow(freezeRange)),getWindows)
        photWindows <- photWindows[unlist(lapply(photWindows,function(x) {length(x) == 401}))]
        photWindows <- do.call(cbind,photWindows)
        colnames(photWindows) <- paste(animal,phasename,seq(1:ncol(photWindows)),sep="_")
        
        ## #########################################################
        ## Freezing Stats by animal
        ## #########################################################
        freezingStats <- data.frame(animal = animal,
                                    phasename = phasename,
                                    seconds = seq(-window,window, length.out = nrow(photWindows)),
                                    average = apply(photWindows,1,mean),
                                    SE = apply(photWindows,1,std.error),
                                    SD = apply(photWindows,1,sd),
                                    percfreeze = round(sum(freeze$Dur...s.)/(180-60) *100,1))
        freezingStats
    }
    
    freezeStats <- lapply(seq(1:nrow(myphase)),traceData,myphase)
    
    phaseAvs <- do.call(cbind,lapply(freezeStats,function(x) {x$average}))
    colnames(phaseAvs) <- paste("M",myphase$ANIMAL,sep= "_")
    
    vals <- data.frame(seconds=freezeStats[[1]]$seconds,
                       phaseAvs,
                       average = rowMeans(phaseAvs),
                       SE = apply(phaseAvs,1,std.error))
    
    write.csv(vals,
              file = paste("./byMouse/",
                           date,"_",phases[[i]],"_",minfreeze,"-",maxfreeze,
                           "_FreezingPhotdata.csv",sep=""),
              row.names=FALSE
              ) 
    
}

freezeVals <- lapply(seq(1:length(phases)),AverageMouse,phases,comp)
