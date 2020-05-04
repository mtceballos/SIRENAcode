#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:19:30 2020

@author: ceballos
"""


def checkLine(lineComplex, lineComplex_fit):
    """ check consistency of theoretical and fitted line
        In particular:
            * ratio of intensities
            * ratio of line centres
            * widths of Lorentzians (fixed)
        Returns:
            1 : if problems found
            0 : if successful
    """
    # check intensity ratios
    nlines = lineComplex.getNumber()
    nlines_fit = lineComplex_fit.getNumber()
    if nlines != nlines_fit:
        print("Inconsistent number of lines:", nlines, nlines_fit)
        return 1
    for i in range(1, nlines):
        amp_0 = lineComplex.rel_amplitudes[0]
        amp_i = lineComplex.rel_amplitudes[i]
        ratioAmp = amp_i/amp_0
        amp_fit_0 = lineComplex_fit.rel_amplitudes[0]
        amp_fit_i = lineComplex_fit.rel_amplitudes[i]
        ratioAmp_fit = amp_fit_i/amp_fit_0
        if abs(ratioAmp-ratioAmp_fit)/ratioAmp > 1e-3:
            print("Inconsistent ratio of amplitudes for", lineComplex.ilabel[i], ":", ratioAmp, ratioAmp_fit)
            return 1

        x0 = lineComplex.energies_eV[0]
        xi = lineComplex.energies_eV[i]
        ratio_x0 = xi/x0
        x0_fit = lineComplex_fit.energies_eV[0]
        xi_fit = lineComplex_fit.energies_eV[i]
        ratio_x0_fit = xi_fit/x0_fit
        # print("Ratio of line centres for", lineComplex.ilabels[i], ":", ratio_x0, ratio_x0_fit)
        # print("Line:", lineComplex.energies_eV[i],lineComplex.energies_eV[0])
        # print("Line fit:", lineComplex_fit.energies_eV[i], lineComplex_fit.energies_eV[0])

        if abs(ratio_x0-ratio_x0_fit)/ratio_x0 > 1e-5:
            print("Inconsistent ratio of line centres for", lineComplex.ilabels[i], ":", ratio_x0, ratio_x0_fit)
            return 1

        if lineComplex.fwhms_eV[i] != lineComplex_fit.fwhms_eV[i]:
            print("Inconsistent fwhm_L for line", lineComplex.ilabels[i], ":",
                  lineComplex.fwhms_eV[i], lineComplex_fit.fwhms_eV[i])
            return 1
    return 0
