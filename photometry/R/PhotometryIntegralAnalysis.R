################################################################################
######################  Photometry Signal Power Analysis  ######################
##
## Calculate Signal Power for photometry readings for each mouse in each phase
##
## Input (should all be in directory set below):
##  - Photometry Files
##     - One file per mouse per phase (recall, extinction, SR, etc)
##     - previously binned using BinTimePoints.R script
##     - File contains time_ext, df_405, df_465 columns
##     - File name should be:
##        - mouse_phase#_phaseName_binned_binningvalue_dF.txt
##        - EX: 10005_phase6_ext4_6_binned_0.01_dF.txt
##     - If you want to analyse 405 traces files have to be in a folder called "photometry"
##
## Output:
##  - XLSX file that lists each file in directory and the power analysis
##    - one file for full directory
##    - File name: IntegralAnalysis.xlsx
## 
################################################################################

#########################################################
## Set parameters
#########################################################
## Loading Libraries
library(plotly)
library(wesanderson)
library(writexl)

## Set up working directory
directory <- "../data/Photometry_Output/"
##dir.create(paste(directory,"output_testing",sep=''))
setwd(directory)

## Build Comparison Info
fls <- list.files("../Photometry_Input/Igor_output", pattern="^[0-9]")
comp <- data.frame(PHASE_NAME = unlist(
                       lapply(fls, function(x) {strsplit(x,split="_")[[1]][[3]][[1]]})
                   ),
                   ANIMAL = unlist(
                       lapply(fls, function(x) {strsplit(x,split="_")[[1]][[1]][[1]]})
                   ),
                   PHASE_NUMBER = unlist(
                       lapply(fls, function(x) {strsplit(x,split="_")[[1]][[2]][[1]]})
                   ))

###################################
## Run Code
###################################
CalciumPower <- function(i,comp) {
    ## Set Data parameters
    phasename <- comp$PHASE_NAME[i]
    animal <- comp$ANIMAL[i]
    phase <- tolower(comp$PHASE_NUMBER[i])
    
    ph <- gsub("([a-z])", "",tolower(phase))
      
    filename <- paste(animal,phase,phasename,sep="_")
    
    ## #########################################################
    ## STEP 1 : load and process data
    ## #########################################################
    ## Load and process photometry data
    photometryTable <- list.files("../Photometry_Input/Igor_output",
                                  pattern="binned", full.names = TRUE)
    photometryTable <- photometryTable[grep(phasename,photometryTable, ignore.case=TRUE)]
    photometryTable <- photometryTable[grep(animal,photometryTable)]
    photometryTable <- photometryTable[grep(phase,photometryTable)]
    
    cat(paste("Analyzing data for : ", photometryTable, "; Testing :", i, "\n",sep=''))
    
    rawPhot <- read.table(photometryTable,sep='\t',header=TRUE)

    ## Doesn't output but in future can offer to output for checking raw data
    p_raw <- plot_ly(rawPhot, x = ~time_ext,color="green") %>%
        add_trace(y = ~dF_465, type = 'scatter', mode = 'lines', name="dF_465",color = I("forestgreen")) %>%
        add_trace(y = ~dF_405, type = 'scatter', mode = 'lines', name="dF_405",color = I("navy"))
  
    ## Remove first minute of data (grabbing the mouse)
    tcutoff <- 60 ## Remove first 60 seconds of experiment

    mcutoff <- max(rawPhot$time_ext)
    phot <- rawPhot[rawPhot$time_ext >= tcutoff & rawPhot$time_ext <= mcutoff,]

    ## Doesn't output plot but can offer option in future for checking outputs
    p_min <- plot_ly(phot, x = ~time_ext,color="green") %>%
        add_trace(y = ~dF_465, type = 'scatter', mode = 'lines', name="dF_465",color = I("forestgreen")) %>%
        add_trace(y = ~dF_405, type = 'scatter', mode = 'lines', name="dF_405",color = I("navy"))
  
    ## #########################################################
    ## STEP 2 : Normalize values
    ## #########################################################
    ## average all values
    av <- mean(phot$dF_465)
    
    ## Normalize by subtracting average photometry reading
    phot$sub_465  <-  phot$dF_465 - av

    ## #########################################################
    ## STEP 3 : Calculate Power (integral of the curve)
    ## #########################################################
    ## Square each value
    sq465 <- phot$sub_465^2
    
    ## Sum all squares together and output (1 number per trace)
    sum(sq465)
}

## Run function
power.fF <- lapply(seq(1:nrow(comp)),CalciumPower,comp)

## dataframe with sample names
power.df <- data.frame(sample = paste(comp$ANIMAL,comp$PHASE_NUMBER,comp$PHASE_NAME,sep="_"),
                       power = unlist(power.fF))
## save dataframe
write_xlsx(power.df,
           path = "./IntegralAnalysis.xlsx")

