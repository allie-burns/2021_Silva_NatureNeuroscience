guiscript=true

/// = CODE DESCRIPTION =
 // Detect cells on (up to) 3 different channels (A,B and C) and finally on DAPI.
 // Check if the centroid of the cell found in the channels (A,B and C)
 // is contained in the region of interest (roi) found on the DAPI channel.
 // Assigned each cell to one of the possible categories (A, B, C, AB, AC, BC, ABC or N)
 // Export downsample image of the annotation with final classification ( and intermediate detection) 
 // Specify a annotationName to do the analysis on a defined annotation on all images of a project.
 // or annotationName="" so it process all the annotations of an image.
 // 
 // == INPUTS ==
 // QuPath project and one open image with annotation(s)
 // Specify the Channels index and the Cell Detection Parameters 
 // 
 // == OUTPUTS ==
 // Detected cells (as detection object) are classified (A, B, C, AB, AC, BC, ABC or N)
 // Output Image contains :
 //     - the DAPI channel with the classified dectection
 //     - the other channels with their intermediate detection
 // 
 // = DEPENDENCIES =
 // BioFormats & BioFormats Extension for QuPath
 // 
 // = INSTALLATION = 
 // No installation. Drag & Drop the groovy script on QuPath and Run!
 // 
 // = AUTHOR INFORMATION =
 // Code written by Romain Guiet, Olivier Burri, Nicolas Chiaruttini, EPFL - SV -PTECH - BIOP 
 // DATE 2018.09.27
 // 
 // = COPYRIGHT =
 // © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, BioImaging And Optics Platform (BIOP), 2018
 // 
 // Licensed under GNU General Public License (GLP) version 3
 // This program is free software: you can redistribute it and/or modify
 // it under the terms of the GNU General Public License as published by
 // the Free Software Foundation, either version 3 of the License, or
 // (at your option) any later version.

 // This program is distributed in the hope that it will be useful,
 // but WITHOUT ANY WARRANTY; without even the implied warranty of
 // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 // GNU General Public License for more details.
 // 
 // You should have received a copy of the GNU General Public License
 // along with this program.  If not, see <https://www.gnu.org/licenses/>.
 ///


import ij.IJ
import ij.ImagePlus
import qupath.imagej.gui.IJExtension
import qupath.lib.regions.RegionRequest
import ij.gui.Overlay
// require for builder in ChDetector
import groovy.transform.builder.Builder


def qp  = getQuPath()
// get the project & list of images
def project =  qp.getProject()
// Get the image data & server
def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def server    = imageData.getServerPath()
//def px_size = getCurrentImageData().getServer().getPixelWidthMicrons()
def annotations = getAnnotationObjects()


def outputDir = buildFilePath(PROJECT_BASE_DIR, 'Output Images')
mkdirs( outputDir)

// Create a set of QP_Class, with defined colors to highlight simple, double or triple positive
QupathClassBuilder.build()
 
// Ensure the ImageJ user interface is showing
IJExtension.getImageJInstance()
IJ.run("Close All", "");
// Declare a ToImagePlus  
// custom class required to output easily channels with Annotation and Detections as overlay)
def tip = new ToImagePlus()

// Declare the  ChDetector, 
// custom class required to define cell detection parameters and store detection 
def chDAPI = new ChDetector()
def chA = new ChDetector()
def chB = new ChDetector()
def chC = new ChDetector()

//////////////////////////////////////////////////////////////////////////////////////////////
//
//
// PARAMETERS THAT YOU HAVE TO SET
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
 
def output_downsample = 2    
 
chDAPI.channelIndex  = 1

// specify Channels Index for :
chA.channelIndex     = 2 // channel A
chB.channelIndex     = 4 // channel B OR set to null to only do a channel A analysis
chC.channelIndex     = null // channel C OR set to null to only do a channels A&B analysis

def annotationName = "CM" // specify a name , i.e. "ORB" or set to "" so it processes all the annotations.

chDAPI.parameters = DetectorParams.builder().detectionImageFluorescence(chDAPI.channelIndex)
                                            .requestedPixelSizeMicrons(0.5)
                                            .includeNuclei(false)
                                            .backgroundRadiusMicrons(8)
                                            .medianRadiusMicrons(0)
                                            .sigmaMicrons(2.5)
                                            .minAreaMicrons(18)                                            
                                            .maxAreaMicrons(400)
                                            .cellExpansionMicrons(3)
                                            .threshold(20)
                                            .build()


chA.parameters = DetectorParams.builder().detectionImageFluorescence(chA.channelIndex)
                                         .requestedPixelSizeMicrons(0.5)
                                         .includeNuclei(false)
                                         .backgroundRadiusMicrons(8)
                                         .medianRadiusMicrons(0)
                                         .sigmaMicrons(2.5)
                                         .minAreaMicrons(20)
                                         .maxAreaMicrons(400)
                                         .threshold(17)
                                         .build()
                                         
chB.parameters = DetectorParams.builder().detectionImageFluorescence(chB.channelIndex)
                                         .requestedPixelSizeMicrons(0.5)
                                         .includeNuclei(false)
                                         .backgroundRadiusMicrons(8)
                                         .medianRadiusMicrons(0)
                                         .sigmaMicrons(5)
                                         .minAreaMicrons(80)
                                         .maxAreaMicrons(400)
                                         .threshold(400)
                                         .build()
                                         
chC.parameters = DetectorParams.builder().detectionImageFluorescence(chC.channelIndex)
                                         .requestedPixelSizeMicrons(0.5)
                                         .includeNuclei(false)
                                         .backgroundRadiusMicrons(8)
                                         .medianRadiusMicrons(0)
                                         .sigmaMicrons(3.4)
                                         .minAreaMicrons(85)
                                         .maxAreaMicrons(600)
                                         .threshold(500)
                                         .build()   



// chDAPI should be the last one in the list,
// so other detection are performed before
def channelToProcess = [chA, chB, chC, chDAPI]



annotations.each{ current_Annotation ->
    setSelectedObject( current_Annotation )
    if ( ( current_Annotation.getPathClass()) ==~ annotationName) { // =~ <=> contains , ==~ looks <=> match

        //get the roi corresponding to the annotation
        // required to process only the cells that are within the annotation
        def roi_Annotation = current_Annotation.getROI()

        // remove existing childOjbects
        if (current_Annotation.hasChildren()) hierarchy.removeObjects(current_Annotation.getChildObjects(), false)

        channelToProcess.each { chDetect ->
            // if childObjects exist => clear them
            if (current_Annotation.hasChildren()) hierarchy.removeObjects(current_Annotation.getChildObjects(), false)

            // detecif ( current_Annotation.getPathClass() != null ){
            current_Annotation_className = current_Annotation.getPathClass()
            print "Current Annotation is : " + current_Annotation_className
            

            // detect cells using defined parameters
            runPlugin('qupath.imagej.detect.nuclei.WatershedCellDetection', '{    "detectionImageFluorescence": ' + chDetect.channelIndex +
                    ', "requestedPixelSizeMicrons": ' + chDetect.parameters.requestedPixelSizeMicrons +
                    ', "backgroundRadiusMicrons": ' + chDetect.parameters.backgroundRadiusMicrons +
                    ', "medianRadiusMicrons": ' + chDetect.parameters.medianRadiusMicrons +
                    ', "sigmaMicrons": ' + chDetect.parameters.sigmaMicrons +
                    ', "minAreaMicrons": ' + chDetect.parameters.minAreaMicrons +
                    ', "maxAreaMicrons": ' + chDetect.parameters.maxAreaMicrons +
                    ', "threshold": ' + chDetect.parameters.threshold +
                    ', "watershedPostProcess": ' + chDetect.parameters.watershedPostProcess +
                    ', "cellExpansionMicrons": ' + chDetect.parameters.cellExpansionMicrons +
                    ', "includeNuclei": ' + chDetect.parameters.includeNuclei +
                    ', "smoothBoundaries": ' + chDetect.parameters.smoothBoundaries +
                    ', "makeMeasurements": ' + chDetect.parameters.makeMeasurements + '}')

            // store the CellObjects in the custom ChDetector
            chDetect.cellsObjects = getCellObjects()

            // when the other channels have been processed
            // the chDetect cooresponds to DAPI
            if (chDetect.channelIndex == chDAPI.channelIndex) {
                // Iterate through the detected cells
                chDetect.cellsObjects.each { cell ->
                    // get the cell roi and check if it belongs to the current Annotation
                    roi_DAPI = cell.getROI()
                    if (roi_Annotation.contains(roi_DAPI.getCentroidX(), roi_DAPI.getCentroidY())) {
                        cell.setPathClass(getPathClass(null))

                        def class_string = ""
                        // check if the cell roi in the DAPI channel contains
                            // the centroid of a cell roi detected in the ChA
                        cells_ChA_Match = chA.cellsObjects.find{ roi_DAPI.contains (  it.getROI().getCentroidX() , it.getROI().getCentroidY() ) }
                        if( cells_ChA_Match != null) class_string += "A"
                        // the centroid of a cell roi detected in the ChB
                        cells_ChB_Match = chB.cellsObjects.find{ roi_DAPI.contains (  it.getROI().getCentroidX() , it.getROI().getCentroidY() ) }
                        if( cells_ChB_Match != null) class_string += "B"
                        // the centroid of a cell roi detected in the ChC 
                        cells_ChC_Match = chC.cellsObjects.find{ roi_DAPI.contains (  it.getROI().getCentroidX() , it.getROI().getCentroidY() ) }
                        if( cells_ChC_Match != null) class_string += "C"
                        // otherwise it is a negative cell
                        if ( class_string == "") class_string += "N"
    
                        // finally we can set the class of that cell
                        cell.setPathClass( getPathClass(class_string) )

                    }
                }
            }


            // send an output image of the channel with Annotation and Detections as overlay
            //select the active channels you want to export
            //tip.setActiveChannels( [chDetect.channelIndex] )
            // Get the imagePlus by defining the annotation, downsample factor and if you want the annotation as a ROI
            //chDetect.imp = tip.getImagePlus( current_Annotation,  output_downsample, true).flatten()
        }
        //tip.setActiveChannels( channelToProcess.collect{ it.channelIndex } )
        //fireHierarchyUpdate()
         // all the channels have been processed
        // we can now show the image corresponding
        //chDAPI.imp.show()
        //chA.imp.show()
        //chB.imp.show()
        //chC.imp.show()
        // Make a stack of it
        //IJ.run("Images to Stack", "name=Image title=Flat")
        //def last = IJ.getImage()
        //last.show()
        // save the images
        //def imageName = getProjectEntry().getImageName() +"_"+ current_Annotation_className+ '.tiff'
        //imagePath = buildFilePath( outputDir, imageName)
        //IJ.saveAs("Tiff", imagePath )
        //print 'Image exported to ' +    imagePath
        // and the measurement
        def measureName = getProjectEntry().getImageName()+ '.csv'
        measurePath = buildFilePath( outputDir, measureName)
        saveAnnotationMeasurements(  measurePath )
        print 'Results exported to '+ measurePath

    }
}// end     
    

print getProjectEntry().getImageName()+" : DONE ! "
   

 
         


//////////////////////////////////////////////////////////////////////////////////////////////
//
//
// CUSTOM CLASSes
//
//
//////////////////////////////////////////////////////////////////////////////////////////////

// Create a set of Class, with defined colors to highlight simple, double or triple positive
class QupathClassBuilder{
    
    static build(){
        def available = getQuPath().getAvailablePathClasses()
    
        def chAPos = getPathClass('A')
        def chBPos = getPathClass('B')
        def chCPos = getPathClass('C')
        
        def doublePosAB = getPathClass('AB')
        def doublePosAC = getPathClass('AC')
        def doublePosBC = getPathClass('BC')
        
        def triplePos = getPathClass('ABC')
        
        def negative = getPathClass('N')
        
        // Add Classes if they do not exist yet
        if (!( chAPos in available))      available.add( chAPos )
        if (!( chBPos  in available))     available.add( chBPos )
        if (!( chCPos  in available))     available.add( chCPos )
        
        if (!( doublePosAB in available)) available.add( doublePosAB )
        if (!( doublePosAC in available)) available.add( doublePosAC )
        if (!( doublePosBC in available)) available.add( doublePosBC )
      
        if (!( triplePos in available))   available.add( triplePos)   
            
        if (!( negative in available))    available.add( negative)
        
        //Define the colors
        chAPos.setColor(getColorRGB(   0, 255, 127))
        chBPos.setColor(getColorRGB( 255, 127, 0  ))
        chCPos.setColor(getColorRGB( 127,   0, 255))
        
        doublePosAB.setColor(getColorRGB(127  , 255, 0  ))
        doublePosAC.setColor(getColorRGB(0    , 127, 255)) 
        doublePosBC.setColor(getColorRGB(255  ,  0 , 127))   
        
        triplePos.setColor(getColorRGB(255,255, 255))   
        
        negative.setColor(getColorRGB(50, 50, 50))
    }
}

// to store parameters and result from detection for cell detection
class ChDetector {
    // the channel of interest
    def channelIndex
    // all the parameters required for " runPlugin('qupath.imagej.detect.nuclei.WatershedCellDetection',...) "
    DetectorParams parameters
    // use to store the cellObjects from the detection
    def cellsObjects
    // use to store the ImagePlus, containing the channel, and the detected cells as overlay
    def imp


}

// use of @Builder, give access to getter and setter
class DetectorParams{
    def detectionImageFluorescence
    def requestedPixelSizeMicrons
    def backgroundRadiusMicrons
    def medianRadiusMicrons
    def sigmaMicrons
    def minAreaMicrons
    def maxAreaMicrons
    def threshold
    def watershedPostProcess
    def cellExpansionMicrons
    def includeNuclei
    def smoothBoundaries
    def makeMeasurements

    static DetectorParamsBuilder builder() {
        new DetectorParamsBuilder()
    }
}


@Builder(builderStrategy=groovy.transform.builder.ExternalStrategy, forClass=DetectorParams)
class DetectorParamsBuilder {
    DetectorParamsBuilder(){
         detectionImageFluorescence = 1
         requestedPixelSizeMicrons  = 1
         backgroundRadiusMicrons    = 10.0
         medianRadiusMicrons        = 2.5
         sigmaMicrons               = 3
         minAreaMicrons             = 20.0
         maxAreaMicrons             = 400.0
         threshold                  = 100
         watershedPostProcess       = true
         cellExpansionMicrons       = 1
         includeNuclei              = true
         smoothBoundaries           = true
         makeMeasurements           = true
    }
}

// This class handles the conversion from a region to an imagePlus with the chosen channels (no LUT info)
class ToImagePlus {
    def imageData
    def server
    def annot
    def hierarchy
    def viewer
    def channels
    def display

    // Gets some constants that we will need
    public ToImagePlus() {

        imageData = getCurrentImageData()
        server = imageData.getServer()

        viewer = getCurrentViewer()
        display = viewer.getImageDisplay()
        channels = display.getAvailableChannels()
    }

    // Sets the channels active in QuPath, which dictates the ones being exported to IJ
    public void setActiveChannels(activeChannels) {
        this.channels.eachWithIndex{ch, idx ->

            this.display.setChannelSelected(ch, ( (idx+1) in activeChannels ) )
        }
    }

    // Exports the given annotation and adds the detections as the overlay
    public ImagePlus getImagePlus(annot, downsample, getAnnotationRoi) {
        hierarchy = getCurrentHierarchy()
        def request = RegionRequest.createInstance(imageData.getServerPath(), downsample, annot.getROI())
        def pathImage =  IJExtension.extractROIWithOverlay(server, annot, hierarchy, request, getAnnotationRoi, viewer.getOverlayOptions(), display)
        def imp = pathImage.getImage()
        if ( getAnnotationRoi ){
            if ( imp.getOverlay() == null) {
                def ovrl = new Overlay(imp.getRoi())
                imp.setOverlay( ovrl )
                //logger.warn("imp.getRoi() : "+imp.getRoi())
            } else{
                imp.getOverlay().add( imp.getRoi() )
            }      
        }
        return imp
    }
}