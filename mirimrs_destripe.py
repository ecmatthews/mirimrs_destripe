import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
import glob
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import matplotlib.ticker as ticker
import re

import jwst
from jwst import datamodels # JWST datamodels
from jwst import pipeline # JWST pipeline
import os
import glob


### some plot settings
plt.style.use('~/tools/elisabeth.mplstyle')

from matplotlib import rcParams
rcParams["axes.grid"] = False
## set label sizes
rcParams["xtick.labelsize"] = 12
rcParams["ytick.labelsize"] = 12
rcParams["axes.labelsize"] = 12


def extract_1d(file, out_dir, input_vars): 
    ''' 
    extract 1d array from cube -- from Helena
    '''
    print('start 1d extraction')
    im = datamodels.IFUCubeModel(file)

    step = pipeline.calwebb_spec3.extract_1d_step.Extract1dStep()
    step.output_dir = out_dir
    
    vars = input_vars['stage3']
    step.save_results = vars['save_results']

    if vars['ifu_autocen']:
        step.ifu_autocen = True #, autocenter of circle where to take the mean
    else:
        step.center_xy = vars['center_x'], vars['center_y'] #24, 29
    
    step.ifu_rfcorr = vars['ifu_rfcorr'] #True  # , residual fringe correction instead in spec2
    step.subtract_background = vars['subtract_background'] #False  # , take ring around as background and subtract, only do this the first time
    step.ifu_rscale = vars['ifu_rscale'] #1  # set number of FWHMs fro radius
    
    print('run step')
    step(im)
    print('done')


class CubeDestripe():
    def __init__(self, file_ifualign, file_wavelengths=None, mask_rad=4, outfolder="."):

        self.outfolder = outfolder

        self.file_ifualign = file_ifualign
        self.file_wavelengths = file_wavelengths
        self.cube_ifualign = fits.open(file_ifualign)
        self.cube_copy = fits.open(file_ifualign)

        ## read in cube, and make a copy -- where we can apply masking
        self.cube_sci = self.cube_ifualign["SCI"].data
        ## slightly hacky solution to all-NAN slices at the end of channel 4
        slice_totals = np.nansum(self.cube_sci, axis=(1,2))
        pos = np.where(slice_totals == 0)[0]
        if len(pos) > 0:
            for p in pos:
                self.cube_sci[p,:,:] = 0

        ## store a copy so that we can mask it
        self.cube_sci_mask = self.cube_sci.copy()

        ## work out what the channel names are
        match = re.search("Level3_ch(\d+)-(\w+)_s3d", file_ifualign)

        self.ch = match.group(1)
        self.chpart = match.group(2)
        chcode = {"short": "A", "medium": "B", "long": "C"}
        self.chpartL = chcode[self.chpart]


        self.diagnostic_plots = True
        if self.diagnostic_plots:
            if not os.path.isdir(os.path.join(self.outfolder,"diagnostic_plots")):
                os.mkdir(os.path.join(self.outfolder,"diagnostic_plots"))

        if not os.path.isdir(os.path.join(self.outfolder,"flat_cubefiles")): ## TODO move this !!!
            os.mkdir(os.path.join(self.outfolder,"flat_cubefiles"))

        ## better outfile name!
        self.out_cube_name = os.path.join(self.outfolder,"flat_cubefiles",file_ifualign.replace("s3d","flatten_s3d").split("/")[-1])

        ## parameters
        self.mask_rad = mask_rad

        ## diagnostics
        if not self.cube_ifualign["SCI"].header["BUNIT"] == "MJy/sr":
            raise ValueError(f"{file_ifualign} is not in MJy/sr")


    # def plot_background_diagnostic(self):

    #     width=0.6
    #     height=width*im_sci.shape[0]/im_sci.shape[1]
    #     barwidth=0.2
    #     gap = 0.03
    #     startx = 0.15
    #     starty = 0.1
    #     axs[0].set_position([startx, starty+gap+barwidth, width, height]) ## main
    #     axs[1].set_position([startx, starty, width, barwidth]) ## bottom
    #     axs[2].set_position([startx+gap+width, starty+gap+barwidth, barwidth, height]) ## right
    #     #
    #     axs[1].yaxis.set_major_locator(plt.MultipleLocator(1_000))
    #     axs[2].xaxis.set_major_locator(plt.MultipleLocator(1_000))
    #     axs[2].tick_params(axis='x', labelrotation=90)
    #     axs[1].set_ylabel("Mean pixel value\n[y direction]")
    #     axs[2].set_xlabel("Mean pixel value\n[x direction]")
    #     #
    #     axs[0].set_xticks([])
    #     axs[0].set_yticks([])
    #     axs[1].set_xticks([])
    #     axs[2].set_yticks([])
    #     axs[3].axis("off")

    #     plt.savefig(f"diagnostic_plots/backgroundvariance_{file[file.index('ch'):file.index('s3d')-1]}_mask{masking}.pdf")


    def mask_source(self):
        """
        Mask the source and its dither subtraction copies in the science cube.
    
        This function identifies the brightest source in the summed image and masks it,
        as well as the two darkest regions (which should correspond to negative images of the source
        from dither subtraction). The masking is applied to the cube_sci_mask attribute.
        """

        ##
        im_sci = np.sum(self.cube_sci, axis=0)

        ## block out the max (target)
        pos = np.array(np.unravel_index(np.nanargmax(im_sci), im_sci.shape))
        self.cube_sci_mask[:, pos[0]-self.mask_rad:pos[0]+self.mask_rad, pos[1]-self.mask_rad:pos[1]+self.mask_rad] = np.nan
        im_sci = np.sum(self.cube_sci_mask, axis=0)

        ## and block out the two minima (dither subtraction copies of target)
        for i in range(2):
            pos = np.array(np.unravel_index(np.nanargmin(im_sci), im_sci.shape))
            print("!!!", pos)
            self.cube_sci_mask[:, pos[0]-self.mask_rad:pos[0]+self.mask_rad, pos[1]-self.mask_rad:pos[1]+self.mask_rad] = np.nan
            im_sci = np.sum(self.cube_sci_mask, axis=0)
        im_sci = np.sum(self.cube_sci_mask, axis=0)

        self.cube_sci_mask_hold = self.cube_sci_mask.copy()

    def destripe(self):
        '''
        Destripe the science cube by subtracting mean values along both axes.

        This method performs destriping on the science cube (self.cube_sci) and its masked version
        (self.cube_sci_mask). It subtracts the mean values along both axes for each slice of the cube.
        The destriped cube (non-masked version) is then saved to a file.
        '''

        for q in range(self.cube_sci_mask.shape[0]):

            ## axis 1
            ymean = np.nanmean(self.cube_sci_mask[q,:,:], axis=1)
            sub1 = np.array(ymean)[:,None]
            sub1 = np.repeat(sub1,self.cube_sci_mask.shape[2],axis=1)
            self.cube_sci_mask[q,:,:] = self.cube_sci_mask[q,:,:] - sub1
            self.cube_sci[q,:,:] = self.cube_sci[q,:,:] - sub1
            
            ## axis 2
            xmean = np.nanmean(self.cube_sci_mask[q,:,:], axis=0)
            sub2 = np.array(xmean)[None,:]
            sub2 = np.repeat(sub2,self.cube_sci_mask.shape[1],axis=0)
            self.cube_sci_mask[q,:,:] = self.cube_sci_mask[q,:,:] - sub2
            self.cube_sci[q,:,:] = self.cube_sci[q,:,:] - sub2


        ## save destriped cube to a file
        self.cube_copy["SCI"].data = self.cube_sci        
        self.cube_copy.writeto(self.out_cube_name,overwrite=True)

        im_sci = np.sum(self.cube_sci, axis=0)
        im_sci_mask = np.sum(self.cube_sci_mask, axis=0)

        if self.diagnostic_plots:

            ext = "_corr"
            for cube in [self.cube_sci_mask, self.cube_sci_mask_hold]:

                im_sci = np.sum(cube, axis=0)

                # Create a figure with a 2x3 grid using gridspec
                fig = plt.figure(figsize=(7, 7))
                gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], width_ratios=[3,1])
                axs = [fig.add_subplot(gs[0, 0])]
                axs.append(fig.add_subplot(gs[1, 0], sharex=axs[0]))
                axs.append(fig.add_subplot(gs[0, 1], sharey=axs[0]))
                axs.append(fig.add_subplot(gs[1, 1]))
                ## 

                im = axs[0].imshow(im_sci, origin="lower", cmap="Greys_r")

                xlin = np.linspace(0, im_sci.shape[1]-1, num=im_sci.shape[1])
                ylin = np.linspace(0, im_sci.shape[0]-1, num=im_sci.shape[0])
                #
                xmean = np.nanmean(im_sci, axis=0)
                ymean = np.nanmean(im_sci, axis=1)
                #
                xstd = np.nanstd(im_sci, axis=0)
                ystd = np.nanstd(im_sci, axis=1)

                # axs[1].plot(xlin, xmean-xstd, "-", color="grey")
                # axs[1].plot(xlin, xmean+xstd, "-", color="grey")
                axs[1].fill_between(xlin, xmean-xstd, xmean+xstd, color="grey")
                axs[1].plot(xlin, np.zeros_like(xlin), "--", color="red", lw=2, zorder=10)
                axs[1].plot(xlin, xmean, "o-", color="black", ms=3)
                #
                axs[2].plot(ymean-ystd, ylin, "-", color="grey")
                axs[2].plot(ymean+ystd, ylin, "-", color="grey")
                axs[2].fill_betweenx(ylin, ymean-ystd, ymean+ystd, color="grey")
                axs[2].plot(np.zeros_like(ylin), ylin, "--", color="red", lw=2, zorder=10)
                axs[2].plot(ymean, ylin, "o-", color="black", ms=3)
                #
                ##
                start = 0
                nparts = 3
                colorpart = ["orange", "limegreen", "blueviolet"]
                
                ## label wavelength range
                xpos = 0.05
                ypos = 0.95

                ## extremely hacky addition of wavelength labels (could be improved! )...
                if not self.file_wavelengths == None:
                    file_s3d = fits.open(self.file_ifualign.replace("s3d","x1d").replace("2301024_ifualign","211123"))
                    wav = np.array(file_s3d["EXTRACT1D"].data["WAVELENGTH"])
                    mywav = wav
                    axs[3].text(xpos, ypos, f"{mywav[0]:.1f}-{mywav[-1]:.1f}$\mu$m", color="black", fontweight="bold", horizontalalignment="left", verticalalignment="top", transform=axs[3].transAxes, fontsize=15)
                    axs[3].axis("off")
                    for n in range(nparts):
                        stop = start + int(cube.shape[0]/nparts)
                        im_sci2 = np.sum(cube[start:stop,:,:], axis=0)
                        xmean2 = np.nanmean(im_sci2, axis=0)
                        xstd2 = np.nanstd(im_sci2, axis=0)
                        ymean2 = np.nanmean(im_sci2, axis=1)
                        ystd2 = np.nanstd(im_sci2, axis=1)
                        axs[1].errorbar(xlin+0.1*(n+1), xmean2, yerr=xstd2, color=colorpart[n], fmt="o-", markersize=3)
                        axs[2].errorbar(ymean2, ylin+0.1*(n+1), xerr=ystd2, color=colorpart[n], fmt="o-", markersize=3)

                        mywav = wav[start:stop]

                        ypos -= 0.22
                        axs[3].text(xpos, ypos, f"{mywav[0]:.1f}-{mywav[-1]:.1f}$\mu$m", color=colorpart[n], horizontalalignment="left", verticalalignment="top", transform=axs[3].transAxes, fontsize=15)

                        # 
                        start = stop


                plt.suptitle(f"{self.ch}{self.chpartL}")

                plt.savefig(os.path.join(self.outfolder,"diagnostic_plots",f"cube_masked_bgval_{self.ch}{self.chpartL}{ext}.pdf"))
                ext = "_nocorr"


    def run_stage3_extract(self, input_vars=None, fringe=False, run=True, ifu_rscale=1.):
        if input_vars is None:
            input_vars = {"stage3": {"ifu_autocen": True,
                            "ifu_rfcorr": False, ## THIS ONE MATTERS FOR THE FRINGING
                            "subtract_background": False,
                            "save_results": True,
                            "coord_system": "ifualign"}}
        input_vars["stage3"]["fringe"] = fringe
        input_vars["stage3"]["ifu_rscale"] = ifu_rscale

        outdir_flat = os.path.join(self.outfolder,"flat_extract1d")
        outdir_basic = os.path.join(self.outfolder,"basic_extract1d")
        if fringe:
            outdir_flat = outdir_flat + "_fringe"
            outdir_basic = outdir_basic + "_fringe"
        else:
            outdir_flat = outdir_flat + "_nofringe"
            outdir_basic = outdir_basic + "_nofringe"
        self.outdir_flat = outdir_flat
        self.outdir_basic = outdir_basic

        #
        if not os.path.isdir(outdir_flat):
            os.mkdir(outdir_flat)
        if not os.path.isdir(outdir_basic):
            os.mkdir(outdir_basic)


        if run:
            ## run the extraction!
            extract_1d(self.out_cube_name, outdir_flat, input_vars)

            ## also run the extraction on the "unflattened" cubes
            extract_1d(self.file_ifualign, outdir_basic, input_vars)

    def plot_spectra_with_flattening(self):
        
        ## note have to run run_stage3_extract before calling this function
        ## (can run it with run=False to just set-up path names)

        if not os.path.isdir(os.path.join(self.outfolder,"plots_speccompare")):
            os.mkdir(os.path.join(self.outfolder,"plots_speccompare"))

        hdul_basic = fits.open(self.outdir_basic+"/"+f"Level3_ch{self.ch}-{self.chpart}_extract1dstep.fits")
        hdul_flat = fits.open(self.outdir_flat+"/"+f"Level3_ch{self.ch}-{self.chpart}_flatten_extract1dstep.fits")

        fig0, ax0 = plt.subplots(2, 1, figsize=(11, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.03})
        chdic = {"1short":0, "1medium":1, "1long":2, "2short":3, "2medium":4, "2long":5, "3short":6, "3medium":7, "3long":8, "4short":6, "4medium":7, "4long":8}

        # fig0, ax0 = plt.subplots(2, 1, figsize=(11, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.03})
        fig1, ax1 = plt.subplots(2, 1, figsize=(11, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.03})

        cc_flat = sns.husl_palette(9, l=0.5)
        cc_basic = sns.husl_palette(9, l=0.7)
        ccix = chdic[f"{self.ch}{self.chpart}"]

        for ax in [ax1]:

            ax[0].errorbar(hdul_basic["EXTRACT1D"].data["WAVELENGTH"], 1000*hdul_basic["EXTRACT1D"].data["FLUX"]+0, yerr=1000*hdul_basic["EXTRACT1D"].data["FLUX_ERROR"], fmt="o-", color=cc_basic[ccix], ms=1)

            ax[0].errorbar(hdul_flat["EXTRACT1D"].data["WAVELENGTH"], 1000*hdul_flat["EXTRACT1D"].data["FLUX"], yerr=1000*hdul_flat["EXTRACT1D"].data["FLUX_ERROR"], fmt="o-", label="ch{}-{}".format(self.ch, self.chpart), color=cc_flat[ccix], ms=1)

            diff = hdul_flat["EXTRACT1D"].data["FLUX"] - hdul_basic["EXTRACT1D"].data["FLUX"]
            ediff = np.sqrt(hdul_flat["EXTRACT1D"].data["FLUX_ERROR"]**2 + hdul_basic["EXTRACT1D"].data["FLUX_ERROR"]**2)

            ax[1].errorbar(hdul_flat["EXTRACT1D"].data["WAVELENGTH"], 1000*diff, yerr=1000*ediff, fmt="o-", color=cc_flat[ccix], ms=1)

        ccix += 1

        ax1[1].set_xlabel("Wavelength [um]")
        ax1[0].set_ylabel("Flux [mJy]")
        ax1[1].set_ylabel("Difference [mJy]")
        w0, w1 = [np.min(hdul_basic["EXTRACT1D"].data["WAVELENGTH"]), np.max(hdul_basic["EXTRACT1D"].data["WAVELENGTH"])]
        ax1[1].set_xlim(w0-(w1-w0)*0.02, w1+(w1-w0)*0.02)

        ax1[1].text(0.01, 0.97, f'Channel {self.ch}{self.chpartL}', horizontalalignment='left', verticalalignment='top', transform = ax1[0].transAxes, fontsize=15)
        fig1.align_ylabels(ax1[:])
        fig1.savefig(os.path.join(self.outfolder,"plots_speccompare",f"spec_{self.ch}{self.chpartL}.pdf"))

        ## TODO bring this back somehow?? plot of all spectra overlaid. Broken since I do a single file at a time in the current setup.
        # ax0[1].set_xlabel("Wavelength [um]")
        # ax0[0].set_ylabel("Flux [mJy]")
        # ax0[1].set_ylabel("Difference [mJy]")
        # ax0[1].set_ylim([-0.5,0.5])
        # w0, w1 = [4.900400095357327, 17.978749789996073]
        # ax0[1].set_xlim(w0-(w1-w0)*0.02, w1+(w1-w0)*0.02)
        # ax0[0].text(0.01, 0.97, f'Channels 1-3', horizontalalignment='left', verticalalignment='top', transform = ax0[0].transAxes, fontsize=15)
        # fig0.align_ylabels(ax0[:])
        # fig0.savefig("plots_speccompare/full_spec.pdf")


if __name__ == "__main__":

    filelist = sorted(glob.glob("../W0458/data_reprocessed2301024_ifualign/*s3d*"))
    print(filelist)

    if not os.path.isdir("example_outputs_w0458"):
        os.mkdir("example_outputs_w0458")

    for i_file, file in enumerate(filelist):
        file_wav = file.replace("s3d","x1d").replace("2301024_ifualign","211123")
        cube_destripe = CubeDestripe(file, file_wavelengths=file_wav, outfolder="example_outputs_w0458")
        cube_destripe.mask_source()
        cube_destripe.destripe()
        cube_destripe.run_stage3_extract(fringe=False, run=True, ifu_rscale=2.)
        cube_destripe.plot_spectra_with_flattening()

        # if i_file > 2:
        #     break

