# fMRI Quality Checker

A set of Matlab and SPM12 scripts to compute several quality metrics and visualizations for fMRI datasets. Influenced by the PCP-QAP and MRIQC.

NOTE: THIS TOOL IS UNDER CONTINUOUS DEVELOPMENT AND HAS NOT BEEN SUFFCIENTLY TESTED

Function to calculate multiple quality control measures for an fMRI time series using SPM12 and Matlab. The goal is to create a tool, for use by fMRI technicians/researchers/clinicians familiar with SPM12 and Matlab, that allows the calculation of measures that can assist in diagnosing quality issues in fMRI data. It is not intended for use as ground truth for quality diagnosis, as these measures are known to be relative and to vary based on scanner site, acquisition time, data format, and more. However, it can provide insight into possible data quality issues originating from the scanner or single subject.

Function steps include preprocessing the acquired fMRI and structural data (coregistering structural image to first functional image, segmenting the
coregistered structural image into tissue types, and reslicing the segments to the functional resolution image grid), and calling functions for a standard set of QC measures based on various sources (mainly including the Quality Assessment Protocol (QAP) from the Preprocessed Connectomes Project (PCP) and MRIQC. This function makes use of SPM12 functions and batch routines. If SPM12 batch parameters are not explicitly set, defaults are assumed. 

CURRENT QC MEASURES (27/06/2018):
- Temporal signal to noise ratio (tSNR) (3D image)
- Mean tSNR of brain, grey matter, white matter, and CSF
- Z-score timeseries and mean
- Standard deviation (3D image)
- Framewise displacement (FD) timeseries, total and mean (mm)
- Non-standardized differential variance (DVARS) timeseries (a.u.)
- The Plot (GM, WM and CSF voxel intensities over time)

INPUT:
- funcional4D_fn: filename of 4D functional timeseries (.nii only)
- structural_fn: filename of T1-weighted structural scan (.nii only)
- fwhm: kernel size for smoothing operations (mm)
- spm_dir: SPM12 directory
- out_dir: output directory (for figures and html logfile)
- subject: subject name/code

OUTPUT:
- Matlab figure: timeseries plots (FD, DVARS, Zscore, The Plot)
- Matlab figures: montage images (tSNR, tSNR brain, stddev, mean EPI)
- Matlab figure: coregistration contour plot
- HTML logfile: HTML logfile with measures and figures and metadata

SOURCES:
- PCP-QAP: http://preprocessed-connectomes-project.org/quality-assessment-protocol/index.html
- MRIQC: https://mriqc.readthedocs.io/en/latest/
- The Plot: https://www.sciencedirect.com/science/article/pii/S1053811916303871?via%3Dihub