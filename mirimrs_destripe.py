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


filelist = sorted(glob.glob("../data_reprocessed211123/*s3d*"))

## PARS
sidecut = 6 ## size of box to find max pixel
##
# radlist = [2,3,4,5,6] ## list of n_fwhm to play with !
radlist = [1, 1.5, 2, 2.5, 3]
colorlist = sns.husl_palette(n_colors=len(radlist)) ## associated color list.

if not os.path.isdir("flat_cubefiles"):
    os.mkdir("flat_cubefiles")
if not os.path.isdir("flat_specfiles"):
    os.mkdir("flat_specfiles")
if not os.path.isdir("plots"):
    os.mkdir("plots")


## TODO check with HK, PP what exact radius is being used

for i_file, file in enumerate(filelist):
    print(file)
    print(file.replace("s3d","x1d"))
    print("")

    # if i_file < 5:
    #     continue


    ## ch4 long is all NAN
    if "ch4-long" in file:
        continue

    cube = fits.open(file)
    spec = fits.open(file.replace("s3d","x1d"))
    spec_2r = fits.open(file.replace("s3d","x1d").replace("211123","170124_ifur2"))
    cube_ifualign = fits.open(file.replace("211123","2301024_ifualign"))
    cube_ifualign_copy = fits.open(file.replace("211123","2301024_ifualign"))

    out_cube_name = "flat_cubefiles/" + file.replace("s3d","s3d_elisabethfix").split("/")[-1]
    out_spec_name = "flat_specfiles/" + file.replace("x1d", "x1d_elisabethfix").replace(".fits",".dat").split("/")[-1]

    if not cube["SCI"].header["BUNIT"] == "MJy/sr":
        raise ValueError(f"{file} is not in MJy/sr")

    if file == filelist[0]:
        cube.info()
        spec.info()
        spec_2r.info()
        print(spec["EXTRACT1D"].columns)
        print(spec_2r["EXTRACT1D"].columns)


    for masking in [0, 1, 2]:

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

        if masking >= 1:
            ## in this case try the autocorrect

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

        if masking == 2:
            im_sci = np.sum(cube_sci_mask_slicewise, axis=0) ## TODO this is still a mess !! IS NOT PLOTTING THE IMAGE CORRECTLY <<<<<

        ## save cube_sci_slicewise!!
        if masking == 2:
            # hdu = fits.PrimaryHDU(data=cube_sci_slicewise, header=cube["SCI"].header)
            # hdu.writeto(out_cube_name, overwrite=False)
            #
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
    plt.clf()

    # continue
    # break

    for i_cube, cubehdul in enumerate([cube, cube_ifualign, cube_ifualign]):
    # for cube_sci in [cube_ifualign["SCI"].data]:

        cube_sci = cubehdul["SCI"].data

        if i_cube == 0:
            title2 = "  SCI"
        if i_cube == 1:
            title2 = "  CORR(CUBE)"
            continue
        if i_cube == 2:
            title2 = "  CORR(SLICE)"

        ## this is a big mess because it's just a uniform correction and the striping *changes*. but ... ???
        if i_cube == 2:
            # cube_sci = cube_sci - sub1
            # cube_sci = cube_sci - sub2 ## WHY DOES THIS NOT WORK i'm doing something jazzy with the sub cubes.
            cube_sci = cube_sci_slicewise

        ## get the cube
        # cube_sci = cube["SCI"].data
        # print(cube["SCI"].header)
        # print(spec["EXTRACT1D"].header)

        ## find the maximum of a small cutout in the center
        s = cube_sci.shape
        box = np.nansum(cube_sci[:,s[1]//2-sidecut:s[1]//2+sidecut, s[2]//2-sidecut:s[2]//2+sidecut], axis=0)
        pos = np.array(np.unravel_index(np.nanargmax(box), (2*sidecut, 2*sidecut)))
        pos_full = pos+np.array((s[1]//2-sidecut, s[2]//2-sidecut))

        ## show the main box and the cutout
        fig, ax = plt.subplots(1,2)
        ax[0].imshow(np.sum(cube_sci, axis=0), origin="lower",cmap="Greys")
        ax[0].plot(pos_full[1], pos_full[0], "x", color="red")
        ax[1].imshow(box, origin="lower",cmap="Greys")
        ax[1].plot(pos[1], pos[0], "x", color="red")

        ## calculate aperture radius - slightly hacky to use the other fits? TODO check with HK, PP
        # wav = cube["WMAP"].data
        # plt.imshow(np.sum(wav, axis=0))
        # plt.show()
        wav_med = np.median(spec["EXTRACT1D"].data["WAVELENGTH"])
        res = 1.22*wav_med/(6.5e6) ## this in radians
        res = res * 180/np.pi ## this in deg
        res = res / cube["SCI"].header["CDELT1"] ## this in pix

        #### VERSION TWO: Law+2023 eq. 1: fwhm = 0.033 * (lambda/micron) + 0.106 arcsec
        res_law2023 = 0.033*(wav_med)+0.106 ## this in arcsec
        res_law2023 = res_law2023 / 60 / 60 ## this in deg
        res_law2023 = res_law2023 / cube["SCI"].header["CDELT1"] ## this in pix
        ## ADOPT THIS ONE:
        res = res_law2023


        ## mark this on the images
        for irad, radval in enumerate(radlist):
            r = [radval*res]*1000
            phi = np.linspace(0, 2*np.pi, 1000)        
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            ax[0].plot(pos_full[1]+x, pos_full[0]+y, "-", color=colorlist[irad])


        ## aperture photometry, slice by slice
        positions = [(pos_full[1], pos_full[0])]
        aper1 = CircularAperture(positions, r=res)
        aper2 = CircularAperture(positions, r=2*res)

        my_spec_1, my_spec_2 = [], []
        my_spec_1_e, my_spec_2_e = [], []


        for z in range(s[0]):
            phot_table = aperture_photometry(cube_sci[z,:,:], aper1, error=cubehdul["ERR"].data[z,:,:])
            my_spec_1.append(phot_table["aperture_sum"].value[0])
            my_spec_1_e.append(phot_table["aperture_sum_err"].value[0])
            phot_table = aperture_photometry(cube_sci[z,:,:], aper2, error=cubehdul["ERR"].data[z,:,:])
            val = phot_table["aperture_sum"].value[0]
            my_spec_2.append(phot_table["aperture_sum"].value[0])
            my_spec_2_e.append(phot_table["aperture_sum_err"].value[0])

        ## conversion factor
        ## area of one pixel
        A_px = abs(float(cube["SCI"].header['CDELT1']) * float(cube["SCI"].header['CDELT2']))
        ## area of one sr
        A_sr = (180. / np.pi) ** 2
        ## conversion
        conversionFactor = (A_px / A_sr) * 1e6

        ## and convert!
        my_spec_1 = np.array(my_spec_1) * conversionFactor
        my_spec_2 = np.array(my_spec_2) * conversionFactor
        my_spec_1_e = np.array(my_spec_1_e) * conversionFactor
        my_spec_2_e = np.array(my_spec_2_e) * conversionFactor

        ## set up xticks, now that you're ready to plot
        xtickval = np.r_[0:cube_sci.shape[0]:80]
        xtickval_label = [f'{spec["EXTRACT1D"].data["WAVELENGTH"][iv]:.2f}' for iv in xtickval]

        ## show both the spectra
        fig, ax = plt.subplots(2,1, sharex=True)
        # ax[0].plot(my_spec_1,"o-",color="magenta",ms=2, label="Elisabeth (1FWHM)", alpha=0.7)
        # ax[0].fill_between(np.linspace(0,len(my_spec_1)-1,num=len(my_spec_1)), my_spec_1-my_spec_1_e, my_spec_1+my_spec_1_e, color="magenta", alpha=0.3)
        ax[0].plot(my_spec_2,"o-",color="red",ms=2, label="Elisabeth (2FWHM)", alpha=0.7)
        # ax[0].fill_between(np.linspace(0,len(my_spec_2)-1,num=len(my_spec_2)), my_spec_2-my_spec_2_e, my_spec_2+my_spec_2_e, color="red", alpha=0.3)
        # ax[0].plot(spec["EXTRACT1D"].data["FLUX"], "o-",color="black",ms=2, label="Helena (1FWHM)", alpha=0.7)
        # ax[0].fill_between(np.linspace(0,len(my_spec_1)-1,num=len(my_spec_1)), spec["EXTRACT1D"].data["FLUX"]-spec["EXTRACT1D"].data["FLUX_ERROR"], spec["EXTRACT1D"].data["FLUX"]+spec["EXTRACT1D"].data["FLUX_ERROR"], color="black", alpha=0.3)
        ax[0].plot(spec_2r["EXTRACT1D"].data["FLUX"], "o-",color="black",ms=2, label="Helena (2FWHM)", alpha=0.7)
        # ax[0].fill_between(np.linspace(0,len(my_spec_1)-1,num=len(my_spec_1)), spec_2r["EXTRACT1D"].data["FLUX"]-spec_2r["EXTRACT1D"].data["FLUX_ERROR"], spec_2r["EXTRACT1D"].data["FLUX"]+spec_2r["EXTRACT1D"].data["FLUX_ERROR"], color="black", alpha=0.3)
        # ax[1].plot(np.array(my_spec_1)/np.array(spec["EXTRACT1D"].data["FLUX"]),"o-",color="blue",ms=2, label="Elisabeth/Helena (1FWHM)")
        ax[1].plot(np.array(my_spec_2)/np.array(spec_2r["EXTRACT1D"].data["FLUX"]),"o-",color="forestgreen",ms=2, label="Elisabeth/Helena (2FWHM)")

        ax[0].set_ylabel("Flux [Jy]")
        ax[0].legend(loc="upper right", fontsize=6)
        ax[1].legend(loc="upper right", fontsize=6)
        ax[1].set_ylabel("ratio Elisabeth_flux/Helena_flux")
        ax[1].set_xlabel("wavelength [um]")
        ax[1].set_xticks(xtickval, labels=xtickval_label, rotation=90)
        ax[1].set_ylim([0.55, 0.85])
        plt.suptitle(file[file.index("ch"):file.index("s3d")-1]+title2)

        fig.savefig(f"plots/comparison_flat_noflat_{file[file.index('ch'):file.index('s3d')-1]}_icube{i_cube}.pdf")


        ## also play with radius
        fig, ax = plt.subplots(2, 1,sharex=True)
        for irad, radval in enumerate(radlist):
            positions = [(pos_full[1], pos_full[0])]
            aperture = CircularAperture(positions, r=res*radval)
            my_spec, my_spec_e = [], []
            for z in range(s[0]):
                phot_table = aperture_photometry(cube_sci[z,:,:], aperture, error=cubehdul["ERR"].data[z,:,:])
                my_spec.append(phot_table["aperture_sum"].value[0])
                my_spec_e.append(phot_table["aperture_sum_err"].value[0])
            my_spec = np.array(my_spec) * conversionFactor
            my_spec_e = np.array(my_spec_e) * conversionFactor

            ax[0].plot(my_spec, "o-", color=colorlist[irad], label=f"rad={radval}*fwhm", ms=2, alpha=0.7)

            if irad == 0:
                spec_ref = np.array(my_spec)
            else:
                ax[1].plot(my_spec/spec_ref, "o-", color=colorlist[irad], ms=2, alpha=0.7)
            ax[1].set_xlabel("wavelength [um]")
            ax[1].set_xticks(xtickval, labels=xtickval_label)

            ## SAVE A SPECTRUM
            if (radval == 2) and (i_cube == 2):
                print("SAVING A SPECTRUM !", out_spec_name)
                f_out = open(out_spec_name, "w")
                f_out.write("lbda    flux[Jy]     e_flux[Jy]\n")
                for lbda, val, val_e in zip(spec["EXTRACT1D"].data["WAVELENGTH"], my_spec, my_spec_e):
                    f_out.write(f"{lbda}   {val}   {val_e}\n")
                    ## TODO NEED WAVELENGTHS AND UNCERTS !!!!
                f_out.close()

        # ax[0].plot(spec["EXTRACT1D"].data["FLUX"], "o-", color="black", label="Helena (1FWHM)", ms=2)
        # ax[0].plot(spec_2r["EXTRACT1D"].data["FLUX"], "o-", color="blue", label="Helena (2FWHM)", ms=2)
        ax[0].legend(loc="upper right", fontsize=8)
        ax[1].set_ylim([-1,3])

        ## also compare basic to ifu_align 0-> this is the 0point for column by column correction
        if i_cube == 0:
            keep_1fwhm = my_spec_1
        else:
            fig, ax = plt.subplots(2, 1, sharex=True)
            ax[0].plot(keep_1fwhm, "o-", color="red", label="1FWHM, DEROT", ms=2)
            ax[0].plot(my_spec_1, "o-", color="black", label="1FWHM, IFUALIGN", ms=2)
            ax[1].plot(keep_1fwhm - my_spec_1, "o-", color="blue", ms=2)
            ax[1].set_ylabel("DEROT - IFUALIGN")
            ax[0].set_ylabel("Flux")
            ax[1].set_xlabel("wavelength [um]")
            ax[1].set_xticks(xtickval, labels=xtickval_label)
            ax[0].legend(loc="upper right", fontsize=6)
            if i_cube == 1:
                keep_1fwhm_2 = my_spec_1
            break
            if i_cube == 2:
                fig, ax = plt.subplots(2, 1, sharex=True)
                ax[0].plot(keep_1fwhm_2, "o-", color="red", label="1FWHM, ifualign", ms=2)
                ax[0].plot(my_spec_1, "o-", color="black", label="1FWHM, ifualign + background_sub", ms=2)
                ax[1].plot(keep_1fwhm_2 - my_spec_1, "o-", color="blue", ms=2)
                ax[1].set_ylabel("(ifualign) - (ifualign+bg_sub)")
                ax[0].set_ylabel("Flux")
                ax[1].set_xlabel("wavelength [um]")
                ax[1].set_xticks(xtickval, labels=xtickval_label)
                ax[0].legend(loc="upper right", fontsize=6)





    ## voila
    # plt.show()
    plt.clf()
    ## TODO NEXT SAVE ALL THE FIG :)


    # if i_file == 5:
    #     break
    # break



## X Y LIM
plt.show()





    ## TODO add uncert







# TODO HK PP is the script you use to go from s3d to x1d shareable? I want to think more about this <<<
    
# TODO HK cuts of spectra at 0 ... that's weird (!) what do we make of that?
    
# TODO save figures !!
    
# TODO a bunch with different radii !!
    
## The radius is quite instrumental to the spectrum ... that *must* indicate that something is going wrong with the background correction.'




## TODO NEXT you have to save the spectrum with some (["ERR"]) uncerts for a ~fair comparison !
## we're getting there, need to plot the uncerts in all the specy-plots. YOUAREHERE