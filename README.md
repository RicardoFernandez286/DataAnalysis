# The DataAnalysis GitHub repository
© 2018-2022, Ricardo J. Fernández-Terán (Ricardo.Fernandez [at] sheffield.ac.uk; Ricardo.FernandezTeran [at] gmail.com)

# 1. Description and contents
This repository contains all the files that constitute three separate main graphical user interfaces (GUIs), now merged into a single window.
The interface was reprogrammed from the old GUIDE-based GUI to a new AppDesigner-based interface, which offers significant advantages.

The first one, **Pump-Probe** is designed to process, analyse and plot time-resolved data from pump-probe experiments
and time-resolved fluorescence spectra (obtained by combining several single-wavelength TCSPC measurements, *support for fluorescence up-conversion spectroscopy will be implemented in a future release*).

The second one, **Interferometric 2D-IR** is designed to process, phase, plot and analyse coherent two-dimensional IR (2D-IR) spectra acquired by spectral interferometry in the pump-probe geometry. The details of the setup used to acquire the data and the main processing algorithm are described in *J. Opt. Soc. Am. B **28**, 171-178 (2011)* [doi: 10.1364/JOSAB.28.000171](https://doi.org/10.1364/JOSAB.28.000171). Additional loading/processing routines were implemented for 2D-IR using pulse shapers at the University of Zurich and the University of Sheffield. 

Finally, **Spectrometer and Shaper calibration** is designed to help with the fundamental task of calibrating the spectral axis of transient absorption spectrometers (both in the IR and UV-Vis regions), and the PhaseTech QuickShape mid-IR shaper. A blank probe spectrum, the spectrum of a solvent/filter and a reference spectrum (measured with a calibrated UV-Vis or FT-IR instrument) need to be provided. Example files can be found in the /Calibration folder for both UV-Vis (Holmium filter) and FT-IR (dioxane and polystyrene). The users are very welcome to provide/suggest additional reference spectra.

# 2. Basic requirements
The first implementation of the program was done in MATLAB R2016a, but current development will always run in the latest release.
The current version of the GUI **DOES NOT RUN** in MATLAB versions older than R2021a.
There are some unexpected bugs in the MacOS version of MATLAB, so please report them to find a workaround in case something does not work properly.
### 2.1 Current MATLAB version: R2022a

A settings file **"GUIoptions.txt"** needs to be present **in the same folder** as the main App (NewGUI.mlapp), in case the user wants to customise their starting datafolders.

See the [About the GUIoptions.txt file](README.md#5-about-the-guioptionstxt-file) section below for more details.

### 2.2 Initial setup
To correctly run the scripts for the first time, download all folders and files and add them (including **all subfolders**) to the MATLAB path.
**IMPORTANT INFO IF UPDATING:**
If updating from previous versions of the GUI, clean the MATLAB path of the "Data Analysis GUI" subfolders and reload it. Many things have changed since then...

# 3. Pump-probe data analysis routine
The initial screen contains buttons for all the features (which are disabled by default, and will be enabled according to the loaded data type).

The program has the following features:
- [x] 2D (contour) and 3D (surface) plotting of the transient absorption/time-resolved emission data.
- [x] Automated background and offset subtraction.
- [x] Shifting of time zero.
- [x] Masking of a probe region (i.e. to avoid plotting pump scatter or to remove noisy regions from spectral/contour plots).
- [x] Plotting of individual (selected) kinetic traces at a given wavelength.
- [x] Plotting of individual (selected) spectral traces at a given delay time.
- [x] Normalisation and plotting of normalised kinetic/spectral traces.
- [x] Plotting of the evolution of kinetic or spectral traces as a function of the scan number (to check signal stability) **with binning options**.
- [x] Plotting of the change in kinetic or spectral traces as a function of the concentration or pump energy.
- [x] Singular value decomposition of the time-resolved data and global fitting of the SVD components.
- [x] Compatibility with the Bruker OPUS files from FTIR spectrometers (to show steady-state data together with transient absorption).
- [x] Compatibility with general steady-state data (UV-Vis and FTIR), to show it together with transient absorption data [must be saved in ASCII format].
- [x] Singular Value Decomposition (SVD) and Global Fit based on selected SVD components.
- [x] Spectrokinetic fit based on <Ref.>
- [x] Automated and manual chirp correction for visible TA data. A new algorithm has been implemented to be used even with raw data (not limited to solvent measurements!).
- [x] Automated spectrometer and pulse shaper calibration routines (fit stretching factor, offset and baseline to known spectrum, for both UV-Vis and mid-IR).
- [x] Automated and manual chirp correction for mid-IR pulse shaper optimisation [THIS HAS BEEN INCLUDED IN A SEPARATE GUI FOR CONVENIENCE].

The unchecked features are still work in progress, incomplete or missing. Any help or comments would be kindly appreciated
Suggestions for new features are very welcome.

The processed pump-probe data can be exported (using the Export... menu) as PDAT [short for *processed data*] files. These files contain the time and wavelength units in the 1st line, then an array of time and wavelength data (the 1st element is a 0 for spacing). This data can be readily plotted/processed/imported into other software.

# 4. 2D-IR data analysis routine
The initial screen contains all the features available so far (disabled by default, and enabled once a suitable dataset is loaded).

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
- [x] Animation of a series of 2D-IR spectra, which can be saved as a GIF file to be used in presentations, etc.
- [ ] Processing, plotting and calculation of 2D-IR anisotropy data.

The unchecked features are still work in progress, incomplete or missing.

# 5. About the "GUIoptions.txt" file
This file contains settings that make the day-to-day operation of the software easier and simplify the user experience.
The default folders need to be defined for the "Go to today..." button to work properly.

The contents of the settings file are defined as follows:
```
defaultTRIRdir	<route to the default TRIR folder>
defaultTAdir	<route to the default electronic TA folder>
defaultSSdir	<route to the default steady-state folder>
defaultFLUPSdir	<route to the default FLUPS folder>
default2DIRdir	<route to the default 2D-IR folder>
defaultSPECdir	<route to the default SEC folder>
defaultPPtype	<indicate the default (starting) data format selected in the Pump-Probe tab -- must match the name on the data type selector>
defaultPPcolormap	<select default colour map for pump-probe -- must match the name on the colormap selector dropdown list>
default2Dcolormap	<select default colour map for 2D-IR -- must match the name on the colormap selector dropdown list>
2DIRContourPlotVersion	<1 or 2; 2 is recommended (newer way to plot contours)>
SelectKnownDataTypes	<true or false; select true to filter out datasets by "known" names (e.g. 2D, etc.)>
SSpecPlotSizePercent	<size of the steady-state plot in % of the transient plot. Does not change size of transient plot; recommended = 60>
SSpecPlotFillArea	<1 or 0; 1 to fill the steady-state spectrum with an area plot>
SSpecPlotFillAreaCol	<RGB triplet in [0,1]  range; colour of the area plot for steady-state spectrum>
AxisBreak	<limits of the linear part of a lin-log kinetic plot- TO BE IMPLEMENTED>
BreakRatio	<relative size of the linear part of a lin-log kinetic plot- TO BE IMPLEMENTED>
LinLogScale	<lin or log; starting format of the pump-probe time axis>
fixScreen	<1 or 0; an attempt to fix the GUI so it fits in smaller screens - TO BE IMPLEMENTED>
```

# License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

# Acknowledgments
I would like to thank my mentors:
* Prof. Dr. Peter Hamm (UZH)
* Dr. Jan Helbing (UZH)
* Prof. Julia A. Weinstein (UoS)

And I would like to thank the following users for finding bugs, requesting features and playing with the software overall:
* Jeannette Ruf (UZH)
* Dr. Gökçen Tek (UZH)
* Dr. Kerstin Oppelt (UZH)
* Dr. James D. Shipp (UoS)
* Martin V. Appleby (UoS)
* Iona Ivalo (UoS)
* Catheryne Royle (UoS)
* Estefanía Sucre-Rosales (UniGE)

