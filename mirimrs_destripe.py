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


filelist = sorted(glob.glob("../W0458/data_reprocessed211123/*s3d*"))

## PARS
sidecut = 6 ## size of box to find max pixel
##
# radlist = [2,3,4,5,6] ## list of n_fwhm to play with !
radlist = [1, 1.5, 2, 2.5, 3]
colorlist = sns.husl_palette(n_colors=len(radlist)) ## associated color list.

if not os.path.isdir("flat_cubefiles"):
    os.mkdir("flat_cubefiles")
if not os.path.isdir("plots"):
    os.mkdir("plots")


## OOP !!!

class CubeDestripe():
    def __init__(self, file_ifualign):

        self.cube_ifualign = fits.open(file_ifualign)
        self.cube_copy = fits.open(file_ifualign)

        ## TODO better out name
        out_cube_name = "flat_cubefiles/" + file_ifualign.replace("s3d","s3d_elisabethfix").split("/")[-1]

        ## diagnostics
        if not self.cube_ifualign["SCI"].header["BUNIT"] == "MJy/sr":
            raise ValueError(f"{file_ifualign} is not in MJy/sr")


    ## function 1: mask pixels!
    def mask_source(self):

        ## DO SOME MASKING on IFUALIGN
        cube_sci = self.cube_ifualign["SCI"].data
        cube_sci_mask = cube_sci.copy()

        mask_rad = 4 ## TODO throughout need to find hardcoding stuff !!! 
        ##
        im_sci = np.sum(cube_sci, axis=0)

        ## block out the max (target)
        pos = np.array(np.unravel_index(np.nanargmax(im_sci), im_sci.shape))
        cube_sci_mask[:, pos[0]-mask_rad:pos[0]+mask_rad, pos[1]-mask_rad:pos[1]+mask_rad] = np.nan
        im_sci = np.sum(cube_sci_mask, axis=0)

        ## and block out the two minima (dither subs)
        for i in range(2):
            pos = np.array(np.unravel_index(np.nanargmin(im_sci), im_sci.shape))
            cube_sci_mask[:, pos[0]-mask_rad:pos[0]+mask_rad, pos[1]-mask_rad:pos[1]+mask_rad] = np.nan
            im_sci = np.sum(cube_sci_mask, axis=0)
            # im_sci[18,24] = 100_000

        self.cube_sci = cube_sci
        self.cube_sci_mask = cube_sci_mask

        ## TODO need a diagnostic plot of cube_sci_mask i thinkkk

    ## function 2: rows, columns, flat
    def destripe(self):

        ## here we destripe

        ## slicebyslice_mask
        for q in range(self.cube_sci_mask.shape[0]):

            ymean = np.nanmean(self.cube_sci_mask[q,:,:], axis=1) ## BLA BLA ax stuff #BLA this is corr (??)
            sub1 = np.array(ymean)[:,None]
            sub1 = np.repeat(sub1,self.cube_sci_mask.shape[2],axis=1)
            self.cube_sci_mask[q,:,:] = self.cube_sci_mask[q,:,:] - sub1
            self.cube_sci[q,:,:] = self.cube_sci[q,:,:] - sub1
            # print(np.nanmedian(cube_sci_mask_slicewise[q,:,:]))
            #
            xmean = np.nanmean(self.cube_sci_mask[q,:,:], axis=0)
            sub2 = np.array(xmean)[None,:]
            sub2 = np.repeat(sub2,self.cube_sci_mask.shape[1],axis=0)
            self.cube_sci_mask[q,:,:] = self.cube_sci_mask[q,:,:] - sub2
            self.cube_sci[q,:,:] = self.cube_sci[q,:,:] - sub2
            # print(np.nanmedian(cube_sci_mask_slicewise[q,:,:]))
            # print("~")

        im_sci = np.sum(self.cube_sci, axis=0) ## TODO this is still a mess !! IS NOT PLOTTING THE IMAGE CORRECTLY <<<<<
        im_sci_mask = np.sum(self.cube_sci_mask, axis=0)

        ## save cube_sci_slicewise!!
        self.cube_copy["SCI"].data = self.cube_sci        
        self.cube_copy.writeto(self.out_cube_name,overwrite=True)



    ## function 3: make the plots
    def diagnostic_plots(self):

        return None ## TODO TODOOOOOO diagnostic plots !!! (split per function ??)

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
        nparts = 3
        start = 0
        # colorpart = ["red", "orange", "green", "blue", "violet", "purple", "pink"]
        colorpart = ["orange", "limegreen", "blueviolet"]
        for n in range(nparts):
            print("!", n)
            print("shape",cube_sci_mask_slicewise.shape)
            stop = start + int(cube_sci_mask_slicewise.shape[0]/nparts)
            print(">",start, stop, "\n")
            im_sci2 = np.sum(cube_sci_mask_slicewise[start:stop,:,:], axis=0)
            xmean2 = np.nanmean(im_sci2, axis=0)
            xstd2 = np.nanstd(im_sci2, axis=0)
            ymean2 = np.nanmean(im_sci2, axis=1)
            ystd2 = np.nanstd(im_sci2, axis=1)
            axs[1].errorbar(xlin+0.1*(n+1), xmean2, yerr=xstd2, color=colorpart[n], fmt="o-", markersize=3)
            axs[2].errorbar(ymean2, ylin+0.1*(n+1), xerr=ystd2, color=colorpart[n], fmt="o-", markersize=3)
            start = stop
        # this adds a scalebar for the images
        # divider = make_axes_locatable(axs[0])
        # cax = divider.append_axes('left', size='5%', pad=0.6)
        # fig.colorbar(im, cax=cax, orientation="vertical", ticklocation="left")
        # cax.set_ylabel("Counts")
        plt.suptitle(file[file.index("ch"):file.index("s3d")-1])

        ## TODO 
        ##
        width=0.6
        height=width*im_sci.shape[0]/im_sci.shape[1]
        barwidth=0.2
        gap = 0.03
        startx = 0.15
        starty = 0.1
        axs[0].set_position([startx, starty+gap+barwidth, width, height]) ## main
        axs[1].set_position([startx, starty, width, barwidth]) ## bottom
        axs[2].set_position([startx+gap+width, starty+gap+barwidth, barwidth, height]) ## right
        #
        axs[1].yaxis.set_major_locator(plt.MultipleLocator(1_000))
        axs[2].xaxis.set_major_locator(plt.MultipleLocator(1_000))
        axs[2].tick_params(axis='x', labelrotation=90)
        axs[1].set_ylabel("Mean pixel value\n[y direction]")
        axs[2].set_xlabel("Mean pixel value\n[x direction]")
        #
        axs[0].set_xticks([])
        axs[0].set_yticks([])
        axs[1].set_xticks([])
        axs[2].set_yticks([])
        axs[3].axis("off")

        plt.savefig(f"plots/backgroundvariance_{file[file.index('ch'):file.index('s3d')-1]}_mask{masking}.pdf")

        # plt.show()
        plt.close("all")



filelist = sorted(glob.glob("../W0458/data_reprocessed211123/*s3d*"))
print(filelist)

for i_file, file in enumerate(filelist):
    cube_destripe = CubeDestripe(file)
    cube_destripe.mask_sources()
    cube_destripe.destripe() ## I have to pass the correct cube here somehow there's some interdependency ... but check this even works, first. TODO !

    if i_file > 3:
        break