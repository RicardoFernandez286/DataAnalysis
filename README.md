# The DataAnalysis GitHub repository
© 2018, Ricardo J. Fernández-Terán (ricardo.fernandez@chem.uzh.ch)

# 1. Description and contents
This project contains all the files that constitute two separate main graphical user interfaces (GUIs).

The first one, **DataAnalysis_GUI** is designed to process, analyse and plot time-resolved data from pump-probe experiments
and time-resolved fluorescence spectra (obtained by combining several single-wavelength TCSPC measurements). Support for streak camera and general time-resolved data will be implemented in a future release.

The second one, **InterfDataAnalysis_GUI** is designed to process, phase, plot and analyse coherent two-dimensional IR (2D-IR) spectra acquired by spectral interferometry in the pump-probe geometry. The details of the setup used to acquire the data and the main processing algorithm are described in *J. Opt. Soc. Am. B **28**, 171-178 (2011)* (https://doi.org/10.1364/JOSAB.28.000171).


# 2. Basic requirements
The first implementation of the program was done in MATLAB R2016a, but current development will run in the latest release.
### 2.1 Current MATLAB version: R2018a

Both GUIs need a file, **"GUIoptions.txt"**, to be present in the root of the __C:__ drive.
This could be changed to a more general location in a later release.
See the [About the GUIoptions.txt file](README.md#5-about-the-guioptionstxt-file) section below for more details].

### 2.2 Initial setup
To correctly run the scripts for the first time, download all folders and files and add them (including all subfolders) to the MATLAB path.

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
- [ ] Fitting of multi-exponential (up to 4) kinetics with (or without) a Gaussian IRF.
- [ ] Fitting of up to 4 stretched exponential kinetics.
- [x] Removal of time-dependent offset from the data, defined as a spectral average or by averaging a user-selected region.
- [x] Plotting the kinetic trace of the difference in transient absorption/time-resolved emission from two wavelengths.

The unchecked features are still work in progress, incomplete or missing.

# 4. 2D-IR data analysis routine
The pump-probe data analysis routine can be called by running the `InterfDataAnalysis_GUI` command from the MATLAB command line.
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
- [ ] Plotting and fitting of "integral dynamics" - a unified module to calculate the sum, difference and so on of different peaks.
- [ ] Spectral diffusion module, which allows analysis of lineshape parameters (CLS, NLS, etc.) and fitting of their evolution.

The unchecked features are still work in progress, incomplete or missing.

# 5. About the "GUIoptions.txt" file
The contents of this file are defined as follows:
```
defaultdir	<route to the default directory for pump-probe data>
defaultIRdir	<route to the default directory for FTIR data>
default2DIRdir	<route to the default directory for 2D-IR data>
defaultSPECdir	<route to the default directory for spectroelectrochemical data>
AxisBreak	<limits to break the axis in a linear|log time scale - format: [start end] of the linear scale (first axis)>
BreakRatio	<size of the linear part in a linear|log time graph - format: x, with 0<x<1>
LinLogScale	<lin or log>
BinScans	<number of scans to bin for the scan-dependent kinetic or spectral traces>
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
* Prof. Dr. Peter Hamm (for the project and helpful discussions).
* Dr. Jan Helbing (for helpful discussions).
* Dr. Kerstin Oppelt (for use, discussions and debugging).
