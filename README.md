# Fixel_segmenter
Identifies fixel orientation from .mif FODF files (MRTrix3)

__**Improved Angular Resolution of Neuronal Fiber Segmentation in Diffusion Magnetic Resonance Imaging**__ 

**Project Description**

This repository takes as input a Fiber Orientation Distribution Function (FODF) .mif file, as well as some supporting MRTrix3 filetype files, and outputs a Fixel.mif file. The algorithm leverages the cylindrical symmetry inherent in single fiber FODFs and the property of FODFs wherein a multifixel - FODF is the summation of single fiber FODFs. The algorithm iteratively scales, rotates, and fits a brain-specific model FODF (a FODF modeling a parallelly oriented fiber bundle) to every FODF in the white matter. It does this (potentially) multiple times per target FODF, until a threshold is reached. The algorithm is thresholded to find fixels with max heights no smaller than 0.1 units, and no smaller than 20% of that voxel’s largest found fixel’s max height. This is to prevent noise effects and ringing effects. Found fixels are also passed through a geometric filter wherein a voxel’s fixel model is considered acceptable when each candidate fixel deviates less than 35 degrees from matching fixels in the nearest-neighbor voxels. The results of this algorithm can be found in the POSTER file. A short summation - this algorithm found ~1.5 times as many fixels as current best segmentation algorithms (peak-finding and SIFT), and found multiple fixels in ~71% of white matter voxels (compared to multiple fixels being found in ~ 45% of WM voxels via peak-finding and SIFT protocols). Information on datasets used for evaluation can be found in the POSTER file as well. My Macbook Air M2 took ~ 4 hours to run this algorithm on a 96x96x60 whole brain dataset with FODF degree 8 (45 coefficients per voxel). I am presenting this algorithm and its results at the BMES annual meeting in October of this year (2023). I am currently working on rewriting all functions in python.
					
**Install and Run / Dependencies**

I wrote this algorithm when I had a student license, and as such had access to about as many MATLAB toolboxes as I wanted. A lot of these functions can be replicated without the toolbox. The toolboxes used are as follows.

- The DSP System Toolbox
  - Dependency for Phased Array System Toolbox
- The Image Processing Toolbox
  - boundary_adjacent_wrapper uses imgradient3 to write simulated fixels normal to edge of brain (used in validation of edge voxels)
- Optimization Toolbox
  - Fixel_wrapper used eqnproblem to fit the model FODF to the cylindrically symmetric cap of the target FODF
- Phased Array System Toolbox
  - Roty and rotz are used for rotating FODFs (used in multiple functions)
- Signal Processing Toolbox
  - Dependency for Phase Array System Toolbox

*Note, three of the toolboxes are needed just for Roty and Rotz. The transform matrices used are 3x3 matrices and really very simple and easily replicated. The formulas can be found on the help page for rotz - https://www.mathworks.com/help/phased/ref/rotz.html

The algorithm also uses a few suites of functions written by others
- Anton Semechko (2023). Suite of functions to perform uniform sampling of a sphere (https://github.com/AntonSemechko/S2-Sampling-Toolbox), GitHub. Retrieved August 14, 2023.
  - This suite is used for uniform sampling of a unit sphere over its surface
  - MIT license
  - https://github.com/AntonSemechko/S2-Sampling-Toolbox/blob/master/LICENSE.md
- Javier Montalt Tordera (2023). Spherical Harmonics (https://github.com/jmontalt/harmonicY/releases/tag/v2.0.1), GitHub. Retrieved August 14, 2023.
  - HarmonicY is used for the evaluation of spherical harmonics of degree N order M in a single direction.
  - MIT license
  - https://github.com/jmontalt/harmonicY/blob/master/LICENSE
- Archontis Politis, Microphone array processing for parametric spatial audio techniques, 2016, Doctoral Dissertation, Department of Signal Processing and Acoustics, Aalto University, Finland
  - getSHrotMtx is used to get a rotated set of SH coefficients 
  - BSD 3-clause License 
  - Copyright (c) 2015, Archontis Politis 
  - https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/LICENSE
- MRTrix3 (2023). https://github.com/MRtrix3/mrtrix3/tree/master, Github. Retrieved August 14, 2023
  - Read_mrtrix and write_mrtrix used to read and write .mif files
  - Mozilla Public License 2.0
  - https://github.com/MRtrix3/mrtrix3/blob/master/LICENCE.txt

**How to use it**

The file pipeline_github.m shows a clear example of how this algorithm is used. First, a white matter FODF file, and a fixel.mif file (peak-finding fixel file) is loaded, as well as a file containing index information about the voxels used in white matter response function estimation. MRtrix3 has more information on all these file types. From here, the steps are as follows

- Create mask using fixel.mif file
- only segment WM voxels with fixels found using other fixel finding means (easy way to delineate WM vs CSF and GM voxels, also done as to be able to compare this algorithm to other algorithms)
- Equidistant sample points on a sphere
- Evaluate spherical harmonics of even degrees from degree 0, order 0 all the way through degree 8 order 8 (45 coefficients) with all coefficients equal to 1. This matrix when multiplied with coefficient vectors allows for a quick and easy way to evaluate future spherical harmonics in the same directions
- Use indices of wm rf voxels to calculate single fiber FODF (average of all (normalized) FODFs used in wm rf estimation, with only m=0 values kept to preserve cylindrical symmetry)
- Simulate fixels around the edge of WM as fixels normal to the edge of the WM (as to help validate edge cases.)
- Segment every FODF in wm of the brain into constituent fixels
- Use geometric validation to filter out erroneous fixels
- Create histograms to compare fixel segmentation results to those of peak-finding and SIFT protocols (optional)
- Write fixel .mif files
  - Use transform matrices from WM FODF file for writing index files, and use transform matrices from fixel direction file for writing directions and AFD files.


**Credits**

This research was conducted under the guidance of Dr. Adam Anderson of Vanderbilt University. It was supported by the VUSRP, and conducted with the Vanderbilt University Institute of Imaging Sciences.

**License**

Mozilla Public License 2.0 (found in LICENSE file)
