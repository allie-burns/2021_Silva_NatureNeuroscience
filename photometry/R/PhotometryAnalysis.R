################################################################################
########################  Photometry-Freezing Analysis  ########################
##
## Analysis of photometry data aligned to freezing information
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
##  - HTML file of freezing mapped onto photometry data
##    - One file per mouse per phase
##    - File name: date_mouse_phase#_phaseName.html
##  - PDF file of average photometry pattern at freezing exit
##    - One file per mouse per phase
##    - File name: date_mouse_phase#_phaseName.pdf
##  - CSV file of photometry and freezing information
##    - Oe file per mouse per phase
##    - File name: mouse_phase#_phaseName_PhotFreezeTable.csv
##    - File contains following columns:
##       - time_ext : binned times from photometry files
##       - df_405
##       - df_465
##       - freeze : boolean (0 = no freeze; 1 = freeze) (freezing = TSE cutoff (0.5s)
##       - freeze_cut : boolean (as above) (freezing = minfreeze (set below))
##  - CSV file of photo signal during each freezing window(freezing end+-2s) for freezings longer than minfreeze(normally set as 1.5s)
## 
################################################################################

#########################################################
## Set parameters
#########################################################
## Loading Libraries
library(plotly)
library(plotrix)
library(RColorBrewer)
library(reshape)
library(gridExtra)
library(scales)

## Set up for saving
date <- gsub('-','',Sys.Date())
## Set location of files and where your output folder will be saved
directory <- "../data/"

## make output directory
outDir <- paste(directory,"Photometry_Output",sep="")
if(!dir.exists(outDir)) { dir.create(outDir) }
setwd(outDir)

## Build Comparison Info
fls <- list.files("../Photometry_Input/Igor_output/",
                  pattern="^[0-9]",
                  include.dirs=FALSE)

comp <- data.frame(
    PHASE_NAME = unlist(lapply(fls, function(x) {
        strsplit(x,split="_")[[1]][[3]][[1]]
    })),
    ANIMAL = unlist(lapply(fls, function(x) {
        strsplit(x,split="_")[[1]][[1]][[1]]
    })),
    PHASE_NUMBER = unlist(lapply(fls, function(x) {
        strsplit(x,split="_")[[1]][[2]][[1]]
    }))
)

write.csv(comp,file=paste("compareInfo_",
                          date,
                          ".csv",
                          sep=''),
          row.names=FALSE)

## Set freezing parameters
tcutoff <- 60   ## Remove first 60 seconds of experiment
minfreeze <- 1.5   ## set minimum freezing length 
maxfreeze <- Inf ## set maximum freezing length 
window <- 2      ## set  window around end of freezing

###################################
## Run code in function
###################################
CalciumPatterns <- function(i,comp) {

    ## Get data parameters for each file
    phasename <- comp$PHASE_NAME[i]
    animal <- comp$ANIMAL[i]
    phase <- comp$PHASE_NUMBER[i]
    ph <- gsub("([a-z])", "",tolower(phase))
        
    cat(paste("Start Analysis for :", animal,"-", phasename,'\n',sep=' '))
    
    ## Load and process Freezing Data
    freezeTable <-list.files("../Photometry_Input/TSE_freezingOutput",pattern="freez", full.names = TRUE)
    freezeTable <- freezeTable[grep(phasename,freezeTable,ignore.case=TRUE)]

    rawFreeze <- do.call(rbind,
                         lapply(freezeTable, function(freezeTable) {
                             read.table(freezeTable,skip=3,dec=",")
                         }))
    columns <- read.table(freezeTable[1],
                          nrow=1,
                          skip=1,
                          sep="\t",
                          stringsAsFactors=FALSE,
                          header=T,
                          check.names=T)
    colnames(rawFreeze) <- colnames(columns)

    ## Edit this one animal to match correct animal number
    rawFreeze$Animal[rawFreeze$Animal == 9554] <- 95540
    
    ## Get freezing information that correlates with each file
    freeze <- rawFreeze[rawFreeze$Animal == animal,] 
    freeze <- freeze[freeze$Phase == ph,] ## filter for recall
    freeze <- freeze[freeze$Time.from..s. >= tcutoff,] 
    freeze_cut <- freeze[freeze$Dur...s. >= minfreeze & freeze$Dur...s. < maxfreeze,]
    
    freezeRange <- data.frame(pre2 = freeze_cut$Time.to..s. - window,
                              post2 = freeze_cut$Time.to..s. + window)
    
    ## Load and process photometry data
    photometryTable <- list.files("../Photometry_Input/Igor_output/", pattern="binned", full.names = TRUE)
    photometryTable <- photometryTable[grep(phasename,photometryTable, ignore.case=TRUE)]

    photometryTable <- photometryTable[grep(animal,photometryTable)] ##error if no matching file name
    photometryTable <- photometryTable[grep(phase,photometryTable)] ##error if no matching phase
    
    cat(paste("Analyzing data for : ", photometryTable, "; Testing :", i, "\n",sep=''))
    
    rawPhot <- read.table(photometryTable,sep='\t',header=TRUE)
        
    ## Add freezing information to photometry table for output
    rawPhot$freeze <- 0
    boolFreeze <- lapply(seq(1:nrow(freeze)),function(i,freeze) {
        x <- freeze[i,]
        fstart <- x$Time.from..s.
        fend <- x$Time.to..s.
        rawPhot$time_ext >= fstart & rawPhot$time_ext <= fend
    },freeze)
    rawPhot$freeze[Reduce("|",boolFreeze)] <- 1
    
    rawPhot$freeze_cut <- 0
    boolFreeze <- lapply(seq(1:nrow(freeze_cut)),function(i,freeze_cut) {
        x <- freeze_cut[i,]
        fstart <- x$Time.from..s.
        fend <- x$Time.to..s.
        rawPhot$time_ext >= fstart & rawPhot$time_ext <= fend
    },freeze_cut)
    rawPhot$freeze_cut[Reduce("|",boolFreeze)] <- 1
    
    write.csv(rawPhot,
              file= paste("./",animal,"_",phase,"_",phasename,"_PhotFreezeTable.csv",sep=""),
              row.names=FALSE)

    if( nrow(freezeRange) >0 ) {
        getWindows <- function(j){
            start <- freezeRange[j,]$pre2
            end <- freezeRange[j,]$post2
            rawPhot[rawPhot$time_ext>= start & rawPhot$time_ext<= end,]$dF_465
        }
        
        photWindows <- lapply(seq(1:nrow(freezeRange)),getWindows)
        photWindows <- photWindows[unlist(lapply(photWindows,function(x) {length(x) == 401}))]##for different binnings not 401
        photWindows <- do.call(cbind,photWindows)
        
        ##if( exp == "photometry" ) {
        get405Windows <- function(i){
            start <- freezeRange[i,]$pre2
            end <- freezeRange[i,]$post2
            rawPhot[rawPhot$time_ext>= start & rawPhot$time_ext<= end,]$dF_405
        }
        
        Windows405 <- lapply(seq(1:nrow(freezeRange)),get405Windows)
        Windows405 <- Windows405[unlist(lapply(Windows405,function(x) {length(x) == 401}))]
        Windows405 <- do.call(cbind,Windows405)
    
    
        ## Get time points for middle of freezing bouts
        MiddleFreezing <- function(i,freeze) {
            dur <- freeze[i,]$Dur...s.
            middleT <- round(freeze[i,]$Time.from..s. + dur/2,digits=2)
            middleT <- data.frame(middleT,
                                  rawPhot[match(middleT,rawPhot$time_ext),],
                                  dur)
        }
        
        midTime <- do.call(rbind,lapply(seq(1:nrow(freeze)),MiddleFreezing,freeze))
        midTime$color <- "grey"
        midTime$color[midTime$freeze_cut == 1] <- "red"
        
        ## #########################################################
        ## Make figure for Ca+2 patterns at end of freezing
        ## #########################################################
        ## get average freezing information for 2 seconds surrounding end of freezing
        freezingStats <- data.frame(seconds = seq(-window,window, length.out = nrow(photWindows)),
                                    photWindows,
                                    average = apply(photWindows,1,mean),
                                    SE = apply(photWindows,1,std.error))
        
        write.csv(freezingStats,
                  paste('./',animal,'_',phase,'_',phasename,'_freezetrace_dF465.csv',sep=''))
        
        ## Rearrange table for plotting
        meltPhot <- melt(freezingStats,id="seconds")
        meltPhot$col <- "cadetblue"
        meltPhot$alpha <- "0.4"
        meltPhot$col[meltPhot$variable == "average"] <- 'black'
        meltPhot$alpha[meltPhot$variable == "average"] <- '1'
        
        ## Output table of freezing stats (for sanity check)
        Stats405 <- data.frame(seconds = seq(-window,window, length.out = nrow(Windows405)),
                               Windows405,
                               average = apply(Windows405,1,mean),
                               SE = apply(Windows405,1,std.error))
        write.csv(freezingStats,
                  paste('./',animal,'_',phase,'_',phasename,'_freezetrace_dF405.csv',sep=''))
    }
    
    ## Plot all traces on top of eachother
    trace <- ggplot(meltPhot, aes(x=seconds,y=value,group=variable)) +
        geom_line(data=subset(meltPhot,variable %in% "average"),color="black",alpha=1)+
        geom_line(data=subset(meltPhot, variable %in% c(as.character(unique(meltPhot$variable[grep("X",meltPhot$variable)])))),color="cadetblue",alpha=0.5) + 
        scale_fill_manual("",values="blue") +
        scale_y_continuous(limits = c(-10, 10)) +
        labs(title = paste("Average dF_465 value for",
                           animal,
                           "during", phasename, "(", phase, ") at the end of freezing", sep=" "),
             y = "Average dF_465 value",
             x = "End of freezing - Time (s)") +
        theme(plot.title = element_text(size=12, face="bold"),
              legend.position = 'none') +
        annotate("text", hjust=0,x=-2,y = 10,##y= max(freezingStats) + 1,
                 label= paste("n =",ncol(photWindows),"\n",
                              minfreeze, " < Freezing Cutoff < ",maxfreeze,sep=" "))
    
    ## Plot traces surrounding freezing as averaged with ribbon plot
    avTrace <- ggplot(data = freezingStats, aes (x = seconds, y = average)) +
        geom_line() +
        geom_ribbon(data = freezingStats,
                    aes(ymin = average - SE, ymax = average + SE, fill = "band"), alpha = 0.3) +
        scale_fill_manual("",values="blue") +
        geom_line(data=Stats405, aes(x = seconds, y = average)) +
        geom_ribbon(data=Stats405,
                    aes(ymin=average - SE, ymax = average + SE,fill = "band"), alpha = 0.3) +
        scale_y_continuous(limits = c(-10, 10)) +
        labs(title = paste("Average dF_465 value for",
                           animal,
                           "during", phasename, "(", phase, ") at the end of freezing", sep=" "),
             y = "Average dF_465 value",
             x = "End of freezing - Time (s)") +
        theme(plot.title = element_text(size=12, face="bold"),
              legend.position = 'none') +
        annotate("text", hjust=0,x=-2, y = 10,
                 ##y= max(freezingStats$average) + max(freezingStats$SE),
                 label= paste("n =",ncol(photWindows),"\n",
                              minfreeze, " < Freezing Cutoff < ",maxfreeze,sep=" "))

    ## Save trace files
    pdf(file   = paste('./',animal,'_',phase,'_',phasename,'_tracingFreezeEnd.pdf',sep=''),
        width  = 15,
        height = 7)
    plot(midTime$dur,midTime$dF_465,
         col=midTime$color, pch=16,
         ylab = "dF_465 at Middle of Freezing",
         xlab = "Duration time for freezing (s)",
         main = paste("Calcium levels at the middle of freezing by duration of freezing \n",
                      animal, "-", phasename))
    abline(lm(data=midTime,formula = dF_465 ~ dur))
    text(max(midTime$dur) - 1,
         max(midTime$dF_465) -1,
         paste("P.value = ",
               scientific(summary(lm(data=midTime,formula = dF_465 ~ dur))$coefficients[8],digits=2),
               "\n", minfreeze, " < Freezing Cutoff < ",maxfreeze,
               sep=''))
    grid.arrange(avTrace,trace,ncol=2)
    ##dev.off()
    
    ## } else {
    ##     pdf(file   = paste(outDir,'/',animal,'_',phase,'_',phasename,'_tracingFreezeEnd.pdf',sep=''),height=7,width=7)
    ##     ggplot() +
    ##         theme(axis.line=element_blank(),axis.text.x=element_blank(),
    ##               axis.text.y=element_blank(),axis.ticks=element_blank(),
    ##               axis.title.x=element_blank(),
    ##               axis.title.y=element_blank(),legend.position="none",
    ##               panel.background=element_blank(),panel.border=element_blank(),
    ##               panel.grid.major=element_blank(),
    ##               panel.grid.minor=element_blank(),plot.background=element_blank()) +
    ##         annotate("text", x=0,y=0,
    ##                  label=paste("No Freezing between ", minfreeze, " and ",maxfreeze," seconds",sep=""))
    ##}
    dev.off()

    ## #########################################################
    ## STEP 3 : Build Figure for Ca+2 patterns and freezing through experiment
    ## #########################################################
    tempindex <- data.frame(xmin = freeze$Time.from..s.,
                            xmax = freeze$Time.to..s.,
                            dur = freeze$Dur...s.)
    
    pal <- c("red", "blue", "green")
    
    if(length(grep("dF_405",colnames(rawPhot))) > 0){
        p <- plot_ly(rawPhot, x = ~time_ext,color="green") %>%
            add_trace(y = ~dF_465, type = 'scatter', mode = 'lines', name="dF_465",color = I("forestgreen")) %>%
            add_trace(y = ~dF_405, type = 'scatter', mode = 'lines', name="dF_405",color = I("navy"))
    } else {
        p <- plot_ly(rawPhot, x = ~time_ext,color="green") %>%
            add_trace(y = ~dF_465, type = 'scatter', mode = 'lines', name="dF_465",color = I("forestgreen"))
    }
    
    ## initiate a line shape object
    rects <- lapply(seq(1:nrow(tempindex)),function(i){
        if(tempindex[i,]$dur <= minfreeze)
            color <- "cadetblue"
        else( color <- "indianred" )
        list(type = "rect",
             fillcolor = color,
             opacity = 0.2,
             xref = "x",
             yref = "y",
             x0 = tempindex[i,]$xmin,
             x1 = tempindex[i,]$xmax,
             y0 = min(rawPhot$dF_465) - 10,
             y1 = max(rawPhot$dF_465) + 1
             )
    })
    
    
    pp <- layout(p,
                 title = paste("Freezing for animal", animal, "during", phasename, "(",phase,")",
                               "\n",minfreeze, " < Freezing Cutoff < ",maxfreeze,sep=" "),
                 shapes = rects)
    
    htmlwidgets::saveWidget(pp,
                            paste('./',animal,'_',phase,'_',phasename,'_tracings.html',sep=''),
                            selfcontained=FALSE)
    
    }
}

lapply(seq(1:nrow(comp)),CalciumPatterns,comp)

writeLines(capture.output(sessionInfo(),paste("Date Run : ",date)),
           paste(".//sessionInfo.txt",sep=""))

