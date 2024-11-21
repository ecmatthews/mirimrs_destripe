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


## TODO check with HK, PP what exact radius is being used

for i_file, file in enumerate(filelist):
    print(file)
    print(file.replace("s3d","x1d"))
    print("")

    cube = fits.open(file)
    spec = fits.open(file.replace("s3d","x1d"))
    spec_2r = fits.open(file.replace("s3d","x1d").replace("211123","170124_ifur2"))
    cube_ifualign = fits.open(file.replace("211123","2301024_ifualign"))
    cube_ifualign_copy = fits.open(file.replace("211123","2301024_ifualign"))

    out_cube_name = "flat_cubefiles/" + file.replace("s3d","s3d_elisabethfix").split("/")[-1]

    if not cube["SCI"].header["BUNIT"] == "MJy/sr":
        raise ValueError(f"{file} is not in MJy/sr")

    if file == filelist[0]:
        cube.info()
        spec.info()
        spec_2r.info()
        print(spec["EXTRACT1D"].columns)
        print(spec_2r["EXTRACT1D"].columns)


    for masking in [2]:

        ## DO SOME MASKING on IFUALIGN
        cube_sci = cube_ifualign["SCI"].data
        cube_sci_mask = cube_sci.copy()
        cube_sci_nomask = cube_sci.copy()
        cube_sci_slicewise = cube_sci.copy()

        mask_rad = 4
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

        cube_sci_mask_slicewise = cube_sci_mask.copy()

        ## here we destripe

        ## first do columns
        ymean = np.nanmean(im_sci, axis=1)
        sub1 = np.array(ymean/cube_sci.shape[0])[None,:,None]
        sub1 = np.repeat(sub1,cube_sci.shape[0],axis=0)
        sub1 = np.repeat(sub1,cube_sci.shape[2],axis=2) ## TODO can those three lines be slicker code?

        cube_sci_mask = cube_sci_mask - sub1
        cube_sci_nomask = cube_sci_nomask - sub1
        # plt.imshow(np.sum(sub, axis=0), origin="lower")
        im_sci = np.sum(cube_sci_mask, axis=0)

        ## then rows
        xmean = np.nanmean(im_sci, axis=0)
        sub2 = np.array(xmean/cube_sci.shape[0])[None,None,:]
        sub2 = np.repeat(sub2,cube_sci.shape[0],axis=0)
        sub2 = np.repeat(sub2,cube_sci.shape[1],axis=1) ## TODO can those three lines be slicker code?

        cube_sci_mask = cube_sci_mask - sub2
        cube_sci_nomask = cube_sci_nomask - sub2
        im_sci = np.sum(cube_sci_mask, axis=0)
        # im_sci[18,24] = 100_000

        ## slicebyslice_mask
        for q in range(cube_sci.shape[0]):

            ymean = np.nanmean(cube_sci_mask_slicewise[q,:,:], axis=1) ## BLA BLA ax stuff #BLA this is corr (??)
            sub1 = np.array(ymean)[:,None]
            sub1 = np.repeat(sub1,cube_sci.shape[2],axis=1)
            cube_sci_mask_slicewise[q,:,:] = cube_sci_mask_slicewise[q,:,:] - sub1
            cube_sci_slicewise[q,:,:] = cube_sci_slicewise[q,:,:] - sub1
            # print(np.nanmedian(cube_sci_mask_slicewise[q,:,:]))
            #
            xmean = np.nanmean(cube_sci_mask_slicewise[q,:,:], axis=0)
            sub2 = np.array(xmean)[None,:]
            sub2 = np.repeat(sub2,cube_sci.shape[1],axis=0)
            cube_sci_mask_slicewise[q,:,:] = cube_sci_mask_slicewise[q,:,:] - sub2
            cube_sci_slicewise[q,:,:] = cube_sci_slicewise[q,:,:] - sub2
            # print(np.nanmedian(cube_sci_mask_slicewise[q,:,:]))
            # print("~")

        im_sci = np.sum(cube_sci_mask_slicewise, axis=0) ## TODO this is still a mess !! IS NOT PLOTTING THE IMAGE CORRECTLY <<<<<

        ## save cube_sci_slicewise!!
        cube_ifualign_copy["SCI"].data = cube_sci_slicewise         
        cube_ifualign_copy.writeto(out_cube_name,overwrite=True)


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

    if i_file >= 2:
        break
    break



## TODO add uncert

# TODO HK PP is the script you use to go from s3d to x1d shareable? I want to think more about this <<<
    
# TODO HK cuts of spectra at 0 ... that's weird (!) what do we make of that?
    
# TODO save figures !!
    
# TODO a bunch with different radii !!
    
## The radius is quite instrumental to the spectrum ... that *must* indicate that something is going wrong with the background correction.'




## TODO NEXT you have to save the spectrum with some (["ERR"]) uncerts for a ~fair comparison !
## we're getting there, need to plot the uncerts in all the specy-plots. YOUAREHERE