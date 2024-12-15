MIRI/MRS destriping tool to remove background systematics
==========

This repository provides a simple fix to background systematics in MIRI/MRS data cubes (*s3d.fits), as described in Matthews et al. (submitted; pre-print on request), Appendix B.


Installation
------------

First clone the github repository:

```sh
git clone mirimrs_destripe
```

Then, navigate to the cloned directory, and install with pip:

```sh
cd mirimrs_destripe
pip install -e .
```

Note that this has not been beta-tested, and there may be issues, e.g. missing dependencies! Please reach out with any problems.


Examples
------------
While the package is not fully documented, we provide some examples for using the code here.

First, define the location of the `s3d.fits` file you wish to correct, and optionally an `x1d` file with the wavelength calibration (the `x1d` file is not used for the systematics removal, but only for a diagnostic plot showing the extent of the systematics).

```python
file = "data_ifualign/Level3_ch1-long_s3d.fits"
file_wavelengths = "data_ifulaign/Level3_ch1-long_x1d.fits"
```

Then, run the analysis. There are several key steps:
*Mask the source, and dithered copies of the source. Here we use the brightest and faintest points in a collapsed cube to mask sources, but a custom approach may be needed especially for other astrophysical sources.
*Run the destriping algorithm on the `s3d.fits` cubes
*Re-run the stage3 extraction from the normal JWST pipeline (please update with your preferred keywords!)
*Make diagnostic plots, showing the spectrum with and without the systematics correction term.

```
cube_destripe = CubeDestripe(file, file_wavelengths=file_wav, outfolder="example_outputs_w0458")
cube_destripe.mask_source()
cube_destripe.destripe()
cube_destripe.run_stage3_extract(fringe=False, run=True, ifu_rscale=2.)
cube_destripe.plot_spectra_with_flattening()
```

When running the code, several files are created in subdirectories of the current directory:
*`flat_cubesfiles` contains `s3d.fits` files, with the background systematics removed
*`flat_extract_1d` contains extracted spectra (`x1d.fits`) from these systematics-corrected cubes
*`basic_extract_1d` contains identically extracted spectra, but without the systematics correction, for comparison

Several diagnostic plots are provided, in the following subdirectories of the current directory:
*`diagnostic_plots` contains images of the cube (collapsed across wavelength), as well as the average pixel value in each row and column. The diagnostics are created both with and without the correction applied; in the case with correction these traces should be close to 0. The average pixel traces are also created for sub-sets of the wavelength range, to indicate how systematics are changing through the cube. **We recommend to check** in these cubes that the source and any negative copies of the source are suitable masked out in these diagnostic images.
*`plots_speccompare` shows a comparison between spectra extracted with and without the systematics corrections.


Improvements
------------
List proposed improvements:
*Beta testing: does installation work correctly? does the code work on other datasets? So far this 
*Better/additional masking algorithms to remove sources from data. Use astrophysical information for this?
*Output plot with all wavelengths plotted (as well as current wavelength-by-wavelength plots)
*Better documentation
*Explore data-driven methods in 2-D detector images (`*cal.fits`) rather than 3-D cubes (`*s3d.fits`, as used here)

This package is still being actively developed. Please reach out if you're interested in contributing!


Attribution
------------

If you make use of this code in your research, please cite Matthews et al. (submitted).