# Silva et al., Nature Neuroscience, 2021 
Code used to perform the analysis presented in *A thalamo-amygdalar circuit underlying the extinction of remote fear memories* (Silva et al., Nature Neuroscience, 2021).

## Abstract
Do you think it's worth it to include the abstract (or a bit more information) about the paper here?  If not, delete this section

## Photometry Analysis
Code provided here was used to analyze the photometry data seen in Figure 3, Figure 5, Supplementary Figure 8, and Supplementary Figure 14. Example files are included in the data folder to help clarify some of the input information for code usage.

  1. `BinTimePoints.R`: Align time points and reduce raw output from fiber photometry apparatus into 10ms bins for downstream processing in Igor.
  2. Igor processing: update file
  3. Photometry Analysis
	 - `PhotometryAnalysis.R`:  Run full photometry analysis for each photometry file.
		- Align photometry and freezing data
		- Plot overlay of photometry and freezing data for full duration of experiment
		- Plot average photometry output for 2 seconds surrounding the end of freezing bouts
	 - `PhotometryAnalysis_byMouse.R`: Calculate the average freezing statistics for each animal in each phase
	 - `PhotometryIntegralAnalysis.R`: Calculate signal power for photometry output	 

## Image Analysis
Code provided here was used to count the number of fluorescently labeled cells on microscopy images, and to convert these counts to the cFos+AAVr+/chance ratios presented in Figures 1d-i and Supplementary Figure 1c-f.

   1. `qps-multiChannels_PositiveCells_Analysis.groovy`: Perform co-localization analysis of AAV2r-, cFos- and Dapi- positive cells. This script outputs one .csv file for each microscopy image, containing the number of Dapi+, Dapi+AAVr+, Dapi+cFos+, and Dapi+cFos+AAVr for each brain region on the image. To run this script, we used the QuPath v0.1.3 software (https://qupath.github.io/).
   2. `CSV_reader.py`: Convert cell counts found with (1) into cFos+AAVr+/chance ratios for each brain region. Ratios are averaged over all brain slices available for each animal. To run this script, we used Python 3.7. An additional requirement is Pandas, which you can install with the command `conda install pandas`.
   3. `QuPath.py`: Contains helper functions necessary for (2).
