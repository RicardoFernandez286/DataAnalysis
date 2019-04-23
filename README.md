# The DataAnalysis GitHub repository
© 2018-2019, Ricardo J. Fernández-Terán (ricardo.fernandez [at] chem.uzh.ch)

# 1. Description and contents
This project contains all the files that constitute two separate main graphical user interfaces (GUIs), which are now merged into a single window.
The interface was reprogrammed from the old GUIDE-based GUI to a new AppDesigner-based interface, which offers significant advantages.

The first one, **DataAnalysis_GUI** is designed to process, analyse and plot time-resolved data from pump-probe experiments
and time-resolved fluorescence spectra (obtained by combining several single-wavelength TCSPC measurements). Support for streak camera and general time-resolved data will be implemented in a future release.

The second one, **InterfDataAnalysis_GUI** is designed to process, phase, plot and analyse coherent two-dimensional IR (2D-IR) spectra acquired by spectral interferometry in the pump-probe geometry. The details of the setup used to acquire the data and the main processing algorithm are described in *J. Opt. Soc. Am. B **28**, 171-178 (2011)* [doi: 10.1364/JOSAB.28.000171](https://doi.org/10.1364/JOSAB.28.000171).


# 2. Basic requirements
The first implementation of the program was done in MATLAB R2016a, but current development will always run in the latest release.
The current version of the GUI **DOES NOT RUN** on earlier MATLAB versions than R2019a.
### 2.1 Current MATLAB version: R2019a

Both GUIs need a file, **"GUIoptions.txt"**, to be present in the same folder as the main App (NewGUI.mlapp).
The file no longer needs to be placed in C: (which was a major drawback of the previous implementation).

See the [About the GUIoptions.txt file](README.md#5-about-the-guioptionstxt-file) section below for more details.

### 2.2 Initial setup
To correctly run the scripts for the first time, download all folders and files and add them (including all subfolders) to the MATLAB path.
**IMPORTANT INFO IF UPDATING:**
If updating from previous versions of the GUI, clean the MATLAB path of the "Data Analysis GUI" subfolders and reload it. Many things have changed since then...

# 3. Pump-probe data analysis routine
The pump-probe data analysis routine can be called by running the `DataAnalysis_GUI` command from the MATLAB command line.
The initial screen contains buttons for all the features (which are disabled by default, and will be enabled according to the loaded data type).

The program has the following features:
- [x] 2D (contour) and 3D (surface) Plotting of the transient absorption/time-resolved emission data.
- [x] Automated background and offset subtraction.
- [x] Shifting of time/wavelength vectors for calibration.
- [x] Plotting of individual (selected) kinetic traces at a given wavelength.
- [x] Plotting of individual (selected) spectral traces at a given delay time.
- [x] Normalisation and plotting of normalised kinetic/spectral traces.
- [x] Plotting of the evolution of kinetic or spectral traces as a function of the scan number (to check signal stability) with binning.
- [x] Plotting of the change in kinetic or spectral traces as a function of the concentration or pump energy.
- [x] Singular value decomposition of the time-resolved data and global fitting of the SVD components.
- [x] Compatibility with the Bruker OPUS files from FTIR spectrometers (to show steady-state data together with transient absorption).
- [x] Fitting of multi-exponential (up to 4) kinetics with (or without) a Gaussian IRF.
- [x] Fitting of up to 4 stretched exponential kinetics.
- [ ] Information about the fits, calculation of the relevant statistical parameters, GooF, etc.
- [x] Removal of time-dependent offset from the data, defined as a spectral average or by averaging over a user-selected region.
- [x] Plotting the kinetic trace of the difference in transient absorption/time-resolved emission between two wavelengths.

The unchecked features are still work in progress, incomplete or missing.

# 4. 2D-IR data analysis routine
The 2D-IR data analysis routine can be called by running the `InterfDataAnalysis_GUI` command from the MATLAB command line.
The initial screen contains all the features available so far.

The program has the following features:
- [x] Automated loading, processing, phasing and plotting of 2D-IR data collected by spectral interferometry in pump-probe geometry.
- [x] Automated signal calculation, with background (scattering) subtraction and pump correction.
- [x] 2D (contour) and 3D (surface) Plotting of the 2D-IR data.
- [x] Automated probe axis calibration from scattering (fitted to a 1st or 2nd degree polynomial).
- [x] Automated (but fully customisable) phasing routine to produce fully absorptive 2D-IR spectra.
- [x] User-selectable apodisation functions (Box and Cos^n with n=1-3).
- [x] Automated zeropadding with user-selectable zeropadding factor.
- [x] Secondary graph showing time-domain data for interferometer ("pixel 0", in V) and for any pixel in the MCT detector (in mOD).
- [x] Plotting of waiting time dependence (kinetics) of any signal in the 2D-IR spectrum.
- [x] Plotting of slices from the 2D-IR data (diagonal, fixed pump or probe WL, integrated pump or probe).
- [x] Plotting and fitting of "integral dynamics" - a unified module to calculate the volume (by integration over a rectangular region), sum, difference and so on of different peaks.
- [x] Processing and plotting of transient 2D-IR data.
- [x] Spectral diffusion module, which allows analysis of lineshape parameters (CLS, IvCLS, NLS, etc.) and fitting of their evolution.
- [x] Volume fitting of 2D-IR peaks by means of a 2D Gaussian lineshape function and extraction of the volume kinetics. [NEW FEATURE]
- [ ] Processing, plotting and calculation of 2D-IR anisotropy data.

The unchecked features are still work in progress, incomplete or missing.

# 5. About the "GUIoptions.txt" file
The contents of this file are defined as follows:
```
defaultdir	<route to the default directory for pump-probe data>
defaultIRdir	<route to the default directory for FTIR data>
default2DIRdir	<route to the default directory for 2D-IR data>
defaultSPECdir	<route to the default directory for spectroelectrochemical data>
defaultPPcolormap	<default colormap for the Pump-Probe data plots, can be: RdOr/Wh/Bl, DkRd/Wh/DkBl, Rd/Wh/Bl, Jet>
default2Dcolormap	<default colormap for the 2D-IR data plots, can be: RdOr/Wh/Bl, DkRd/Wh/DkBl, Rd/Wh/Bl, Jet>
2DIRContourPlotVersion	<1 or 2 / version 2 is an attempt to improve the contour representations and to make nicer contour plots with black lines. v1 is recommended for now.>
SelectKnownDataTypes	<true or false / will keep only datafolders containing "pp" or "2D" in their names in the corresponding tabs>
AxisBreak	<limits to break the axis in a linear|log time scale - format: [start end] of the linear scale (first axis) **- TO BE IMPLEMENTED**> 
BreakRatio	<size of the linear part in a linear|log time graph - format: x, with 0<x<1 **- TO BE IMPLEMENTED**>
LinLogScale	<can be off or log. If set to log, will create a lin|log time graph, otherwise a lin|lin or log|log graph without break is created **- TO BE IMPLEMENTED**>
```

# License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

# Acknowledgments
* Prof. Dr. Peter Hamm (for mentoring and helpful discussions).
* Dr. Jan Helbing (for helpful discussions).
* Andrea Pasti (for using this software, helpful discussions and for finding bugs).
* Gökçen Tek (for using this software, helpful discussions and for finding bugs).
* Dr. Kerstin Oppelt (for using this software, helpful discussions and for finding bugs).

