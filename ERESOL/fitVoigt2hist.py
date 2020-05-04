#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:21:45 2020

@author: ceballos
"""

from astropy.modeling import models, fitting
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from RxLines import RxLines
from checkLine import checkLine
from ecdf import ecdf
from cdf_interp import cdf_interp

def fitVoigt2hist(data=None, lines=None, nbins=None, ax0=None):
    """
    Fit Voigt profiles (according to description in lines to data histogram) or
    compute curve of residuals vs. number of bins
    Histogram is created and plotted with matplotlib.pyplot.hist in Density
    Voigt functions from astropy.fitting module (fitting by LevMarLSQFitter)

    data: (array) data for 1st histogram
    lines: (RxLines) lines complex for histo
    nbins: (array) number of bins for histogram.
        If single value, perform fit and plot for this number of bins
        If array: do fit for each binning, calculate CDF for ech fit and
                  residuals. Plot comparative of CDF and residuals plot
    ax0: (object) ax object from figure to plot results
    outfig: file to save final figure

    returns:
        fwhm_G: fwhm of Gaussian broadening (if run in single model (scalar nbins))

    Voigt relative intensities keep tied
    Lorentzian FWHMs are kept fixed
    Gaussian broadening are the same for all lines
    """

    # define functions and function-builders to tie parameters
    # Tie FWHM of Gaussians:
    def tie_gauss(model):
        return model.fwhm_G_0

    # Tie amplitudes of Voigt
    def tie_amplitude_builder(line, i):
        def tie_ampl(model):
            rel0 = line.rel_amplitudes[0]
            reli = line.rel_amplitudes[i]
            return model.amplitude_L_0 / rel0 * reli
        return tie_ampl

    # Tie line centres
    def tie_x0_builder(line, i):
        def tie_x0(model):
            return model.x_0_0 / line.energies_eV[0] * line.energies_eV[i]
        return tie_x0


    def fitVoigtLines(nbins=0, density=False, plotFit=False):
        """Do the fitting for the given parameters:

        nbins: number of bins for histogram.
        density: (logical) for histogram
        plotFit: if fit is to be plot

        returns:
            fitted model (astropy.fitting Voigt profiles object model)

        Voigt relative intensities keep tied
        Lorentzian FWHMs are kept fixed
        Gaussian broadening are the same for all lines

        Free parameters: amplitude of MnKa11
                         line_centre of MnKa11
                         gaussian_broadening
        """

        # create histogram for data
        #bin_heights, bin_borders = np.histogram(data, bins=nbins,  density=True)
        bin_heights, bin_borders = np.histogram(data, bins=nbins,  density=density)
        bin_width = bin_borders[1]-bin_borders[0]
        bin_width = '{:0.2f}'.format(bin_borders[1]-bin_borders[0])
        if plotFit:
            ax0.hist(data, bins=bin_borders, density=density, alpha=0.5, label="Histogram")
        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)


        nlines = lines.getNumber()
        # create dictionary of functions for L_amplitude ties
        dict_tie_ampl = {}
        for i in range(0, nlines):
            func_name = "tie_ampl_" + str(i)
            dict_tie_ampl[func_name] = tie_amplitude_builder(lines, i)

        # create dict of functions for centres ties
        dict_tie_x0 = {}
        for i in range(0, nlines):
            func_name = "tie_x0_" + str(i)
            dict_tie_x0[func_name] = tie_x0_builder(lines, i)

        # do first fit with initial values set to those of calibration table
        # -------------------------------------------------------------------
        # define models (initial values, fixed params and ties)
        Vmods = list()
        Vmods.append(models.Voigt1D(x_0=lines.energies_eV[0], amplitude_L=lines.rel_amplitudes[0],
                                    fwhm_L=lines.fwhms_eV[0], fwhm_G=2., fixed={'fwhm_L': True}))
        #Vmods[0].fwhm_G.min = 1.e-6
        sumVoigt = Vmods[0]
        # print("Model=", sumVoigt)
        for i in range(1, nlines):
            Vmods.append(models.Voigt1D(x_0=lines.energies_eV[i], amplitude_L=lines.rel_amplitudes[i],
                                        fwhm_L=lines.fwhms_eV[i], fwhm_G=2., fixed={'fwhm_L': True}))
            Vmods[i].fwhm_G.tied = tie_gauss
            Vmods[i].amplitude_L.tied = dict_tie_ampl["tie_ampl_" + str(i)]
            Vmods[i].x_0.tied = dict_tie_x0["tie_x0_" + str(i)]

            sumVoigt += Vmods[i]
            # print("Adding Model=", Vmods[i])

        # fit nlines1 Voigts to density-histogram
        fitter = fitting.LevMarLSQFitter()
        vv_fit = fitter(sumVoigt, bin_centers, bin_heights, maxiter=300)

        # check there are no errors in fitting
        if not fitter.fit_info['ierr'] in [1,2,3,4]:
            print("Solution not found for nbins=", nbins)
            print("Message (", lines.complabel, ")=", fitter.fit_info['message'])

        # check fitting (ratios ,etc)
        ilabels_fit = [s + "_fit" for s in lines.ilabels]
        energies_eV_fit = np.zeros(nlines, dtype=np.float64)
        fwhms_eV_fit = np.zeros(nlines, dtype=np.float64)
        rel_amplitudes_fit = np.zeros(nlines, dtype=np.float64)

        # plot individual line and global model
        for i in range(nlines):
            # print("Line:", lines1.ilabels[i])
            energies_eV_fit[i] = vv_fit.param_sets[i*4][0]
            rel_amplitudes_fit[i] = vv_fit.param_sets[i*4+1][0]
            fwhms_eV_fit[i] = vv_fit.param_sets[i*4+2][0]
            fwhm_G = vv_fit.param_sets[i*4+3][0]
            # print("x0=", energies_eV_fit[i])
            # print("Ampl=", rel_amplitudes_fit[i])
            # print("FWHM_L=", fwhms_eV_fit[i])
            # print("FWHM_G=", fwhm_G)
            v = models.Voigt1D(x_0=energies_eV_fit[i], amplitude_L=rel_amplitudes_fit[i],
                           fwhm_L=fwhms_eV_fit[i], fwhm_G=fwhm_G)
            if plotFit:
                colorStr = "C" + str(int(i)+1)
                ax0.plot(x_interval_for_fit, v(x_interval_for_fit), color=colorStr, label="Fit " + lines.ilabels[i])
                ax0.axvline(lines.energies_eV[i], ls="--", color=colorStr, label="CalibVoigt " + lines.ilabels[i])

        # calculate model-data residuals
        min_residuals = np.sum((bin_heights-vv_fit(bin_centers))**2)

        param_cov = fitter.fit_info['param_cov']
        if (param_cov is None):
            print("Matrix for ", lines.complabel, " is singular for nbins=", nbins)
            err_fwhm_G = 0.

        else:
            npars = param_cov.shape[0]
            err_centre0 = np.sqrt(param_cov[0, 0])
            err_ampl0 = np.sqrt(param_cov[0, 0])
            err_fwhm_G = np.sqrt(param_cov[(npars-1), (npars-1)])

            iter = 0
            # check fitting in different initial conditions: -3sigma, +3sigma for each free param
            for x0init in (energies_eV_fit[0]-3*err_centre0, energies_eV_fit[0]+3*err_centre0):
                for aLinit in (rel_amplitudes_fit[0]-3*err_ampl0, rel_amplitudes_fit[0]+3*err_ampl0):
                    for fwhmGinit in (fwhm_G-3*err_fwhm_G, fwhm_G+3*err_fwhm_G):
                        iter += 1
                        Vmods = list()
                        Vmods.append(models.Voigt1D(x_0=x0init, amplitude_L=aLinit, fwhm_L=lines.fwhms_eV[0],
                                                    fwhm_G=fwhmGinit, fixed={'fwhm_L': True}))

                        sumVoigt = Vmods[0]
                        for i in range(1, nlines):
                            x0_i = x0init / lines.energies_eV[0] * lines.energies_eV[i]
                            aL_i = aLinit / lines.rel_amplitudes[0] * lines.rel_amplitudes[i]
                            Vmods.append(models.Voigt1D(x_0=x0_i, amplitude_L=aL_i, fwhm_L=lines.fwhms_eV[i],
                                                        fwhm_G=fwhmGinit, fixed={'fwhm_L': True}))
                            Vmods[i].fwhm_G.tied = tie_gauss
                            Vmods[i].amplitude_L.tied = dict_tie_ampl["tie_ampl_" + str(i)]
                            Vmods[i].x_0.tied = dict_tie_x0["tie_x0_" + str(i)]
                            sumVoigt += Vmods[i]

                        fitter = fitting.LevMarLSQFitter()
                        vv_fit_tmp = fitter(sumVoigt, bin_centers, bin_heights, maxiter=300)
                        resid2 = np.sum((bin_heights-vv_fit_tmp(bin_centers))**2)
                        #print("Min residuals=", min_residuals)
                        #print("New residuals=", resid2)

                        if resid2 < min_residuals:
                            min_residuals = resid2
                            vv_fit = vv_fit_tmp
                            #get fit parameters
                            for i in range(nlines):
                                energies_eV_fit[i] = vv_fit.param_sets[i*4][0]
                                rel_amplitudes_fit[i] = vv_fit.param_sets[i*4+1][0]
                                fwhms_eV_fit[i] = vv_fit.param_sets[i*4+2][0]
                                fwhm_G = vv_fit.param_sets[0*4+3][0]

                            #print("New lower residuals for nbins=", nbins, "and iteration=", iter)
                            param_cov = fitter.fit_info['param_cov']
                            if (param_cov is None):
                                print("Matrix for ", lines.complabel, " is singular for nbins=", nbins)
                                err_fwhm_G = 0.
                            else:
                                npars = param_cov.shape[0]
                                err_fwhm_G = np.sqrt(param_cov[(npars-1), (npars-1)])

        # check fitting consistency
        lines_fit = RxLines(complabel=lines.complabel + "_fit", ilabels=ilabels_fit, energies_eV=energies_eV_fit,
                         fwhms_eV=fwhms_eV_fit, rel_amplitudes=rel_amplitudes_fit)
        status = checkLine(lines, lines_fit)
        if status:
            raise RuntimeError
            print("Line consistency check status:", status, "(0 = OK)")

        if plotFit:
            ax0.plot(x_interval_for_fit, vv_fit(x_interval_for_fit), label='Voigt fit', color='black')
            #maxy = max(bin_heights)
            #xtxt = min(data)
            #ax0.text(xtxt, maxy-0.01, "FWHM_G=" + '{:0.2f}'.format(fwhm_G) + "+/-" + '{:0.2f}'.format(err_fwhm_G) + "eV")
            if density:
                ylabel = "Density"
            else:
                ylabel = "counts/" + bin_width + "eV bin"
            ax0.set_ylabel(ylabel, fontsize='x-large')
            ax0.tick_params(axis='both', which='major', labelsize=14)

            ax0_divider = make_axes_locatable(ax0)
            axb = ax0_divider.append_axes("bottom", size="20%", pad=0.5)
            axb.plot(bin_centers, vv_fit(bin_centers)-bin_heights, marker='.', color="black", ls="")
            axb.axhline(0., color="gray", ls="--")
            axb.set_xlabel("Energy (eV)", fontsize='x-large')
            axb.set_ylabel("Residual", fontsize='x-large')
            axb.tick_params(axis='both', which='major', labelsize=14)

            #print("Total counts=", numberCounts)
            ax0.legend()

        return (fwhm_G, err_fwhm_G, vv_fit)


    if isinstance(nbins, int):  # if just a single value for nbins

        fwhm_G, err_fwhm_G, vv_model_fit = fitVoigtLines(nbins=nbins, density=False, plotFit=True)
        return (fwhm_G, err_fwhm_G, vv_model_fit)

    else:  # calculate residuals in diff of CDFs for different number of bins
        # calculate ECDF of data
        x_ecdf, y_ecdf = ecdf(data)
        residuals = np.zeros(len(nbins))

        for ib in range(0, len(nbins)):
            # calculate CDF of model fit
            fwhm_G, err_fwhm_G, vv_model_fit = fitVoigtLines(nbins=nbins[ib], density=True, plotFit=False)
            vv_model_cdf_int = cdf_interp(vv_model_fit, min(data)-30, max(data)+30, 10000)
            model_cdf = vv_model_cdf_int(x_ecdf)
            #ax.plot(x_ecdf, y_ecdf)
            #ax.plot(x_ecdf, model_cdf)
            residuals[ib] = np.sum((y_ecdf-model_cdf)**2)
        ax0.plot(nbins, residuals)
        ax0.set_xlabel("nbins")
        ax0.set_ylabel("residuals: SUM(data_ECDF-model_CDF)²")
