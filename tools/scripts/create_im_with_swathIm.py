#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-.

This software is released under a three-clause BSD license:
 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of any author or any participating institution
   may be used to endorse or promote products derived from this software
   without specific prior written permission.
For a full list of authors, refer to the file AUTHORS.
--------------------------------------------------------------------------
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------
$Maintainer: Hannes Roest, Joshua Charkow $
$Authors: Hannes Roest, Joshua Charkow $
--------------------------------------------------------------------------
"""
from __future__ import division
from __future__ import print_function

# Create simulated ion mobility scans for testing

from pyopenms import *
import pyopenms

exp = MSExperiment()


##################################### PARAMETERS ##########################################
NR_RT_SAMPLES = 50
IM_BINS = (300, 700)
NR_PEAKS = 5

# create MS2 precursor across IM for a single spectrum
# allmz, allint, allim = lists storing all of the information for the spectrum
# imCenterIdx = idx IM value, will be divided by 500
# imWidth = width of IM in idx units
# base_int = base intensity for the spectrum
# baseFragmentMz = lowest fragmentMz
# nr_peaks = number of fragment peaks to create
# imWinLim = window limits of IM, do not plot anything outside of this window
def createMS2Precursor(allmz, allint, allim, imCenterIdx, imWidth, base_int, baseFragmentMz = 100, nr_peaks = NR_PEAKS, im_bins=IM_BINS, imWinLim=None):
    if imWinLim is None:
        for im_idx in range(*im_bins):
            for i in range(nr_peaks):
                if (imCenterIdx - imWidth) < im_idx < (imCenterIdx + imWidth):
                    apex_dist = abs( imCenterIdx - im_idx)
                    p = Peak1D()
                    p.setMZ(baseFragmentMz+i)
                    p.setIntensity(base_int * (i + 1) - base_int * (i + 1) * apex_dist / imWidth)
                    allmz.append(p.getMZ())
                    allint.append(p.getIntensity())
                    allim.append( im_idx / 500.0)
    else: # have to check if in range of IM upper limit
        for im_idx in range(*im_bins):
            for i in range(nr_peaks):
                if (imCenterIdx - imWidth) < im_idx < (imCenterIdx + imWidth) and (imWinLim[0] < im_idx < imWinLim[1] ):
                    apex_dist = abs( imCenterIdx - im_idx)
                    p = Peak1D()
                    p.setMZ(baseFragmentMz+i)
                    p.setIntensity(base_int * (i + 1) - base_int * (i + 1) * apex_dist / imWidth)
                    allmz.append(p.getMZ())
                    allint.append(p.getIntensity())
                    allim.append( im_idx / 500.0)


# add values to spectrum across IM corresponding with a precursor
# sp = input spectrum
# fda = input float data array
# mz = mz of precursor
# im = im of precursor (in idx unit will be divided by 500)
# base_intensity = intensity to base other intensity calculations off of
# im_offset_idx = deviation in ion mobility (each side) in idx units
# nr_peaks = number of peaks to create
def addMS1Precursor(sp, fda, mz, im_idx, base_intensity, im_offset_idx=2, nr_peaks=10):
    for i in range(nr_peaks):
        p = Peak1D()
        p.setMZ(mz)
        p.setIntensity(base_intensity+ (i-5))
        sp.push_back(p)

        fda.push_back( (im_idx - im_offset_idx + i*im_offset_idx)/500)


########### 1. Create MS1 spectra #####################


## create peaks from 100m/z to 200m/z
for rt_idx in range(NR_RT_SAMPLES):

    # compute base intensity
    base_int = 100 - abs(20 - rt_idx)*(100)/25.0  # shift the precursor slightly
    base_int = max(5, base_int)


    sp = MSSpectrum()
    sp.setRT(rt_idx)
    sp.setMSLevel(1)
    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(100)

    for i in range(100):
        p = Peak1D()
        p.setMZ(100+i)
        p.setIntensity(base_int+i)
        sp.push_back(p)
        fda[i] = 300/500
    
    # add precursors that are in the library 
    addMS1Precursor(sp, fda, 412.502, 350, base_int)
    addMS1Precursor(sp, fda, 417.502, 450, base_int, im_offset_idx=5)
    addMS1Precursor(sp, fda, 422.502, 550, base_int)

    sp.setFloatDataArrays([fda])
    exp.addSpectrum(sp)


###### CREATE MS2 SPECTRA FOR FIRST SWATH ##########
rt_offset = 0.2
for rt_idx in range(NR_RT_SAMPLES):

    base_int = 100 - abs(25 - rt_idx)*(100)/25.0
    base_int_second = 100 - abs(10 - rt_idx)*(100)/40.0

    allmz = []
    allint = []
    allim = []
    
    # peaks of a precursor at 412.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
    # and ion mobility 100
    createMS2Precursor(allmz, allint, allim, imCenterIdx=350, imWidth=10, base_int=base_int, baseFragmentMz=100)
    # peaks of a precursor at 417.5 m/z : 100.01, 101.01, 102.01, .. 100 + NR_PEAKS. Slight offset so chromatograms can be sorted unambiguously 
    # and ion mobility 150
    createMS2Precursor(allmz, allint, allim, imCenterIdx=450, imWidth=20, base_int=base_int, imWinLim=(300, 450), baseFragmentMz=100.01)

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(len(allmz))
    for k,val in enumerate(allim):
        fda[k] = val
    
    ### create spectrum and set metaData
    sframe = pyopenms.MSSpectrum()
    sframe.setMSLevel(2)
    sframe.setRT(rt_idx + rt_offset)

    ## set IM window boundaries # precursor m/z 417.5 is in this window but not centered, thus should not be extracted from this window
    sframe.setMetaValue('ion mobility lower limit', 300/500)
    sframe.setMetaValue('ion mobility upper limit', 450/500)

    p = pyopenms.Precursor()
    center = 412.5
    width = 25
    p.setMZ(center)
    p.setIsolationWindowUpperOffset(width / 2.0)
    p.setIsolationWindowLowerOffset(width / 2.0)
    sframe.setPrecursors([p])

    # store data in spectrum
    sframe.set_peaks( (allmz, allint) )
    sframe.setFloatDataArrays([fda])
    sframe.sortByPosition()
    exp.addSpectrum(sframe)

############# CREATE MS2 SPECTRA FOR 2ND SWATH #########################
rt_offset = 0.3
for rt_idx in range(NR_RT_SAMPLES):

    base_int = 100 - abs(25 - rt_idx)*(100)/25.0
    base_int_second = 100 - abs(10 - rt_idx)*(100)/40.0

    allmz = []
    allint = []
    allim = []
    
    # peaks of a precursor at 412.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
    # and ion mobility 100
    createMS2Precursor(allmz, allint, allim, imCenterIdx=350, imWidth=10, base_int=base_int, imWinLim=(350, 600))

    # peaks of a precursor at 417.5 m/z : 100.01, 101.01, 102.01, .. 100.01 + NR_PEAKS
    # and ion mobility 150

    createMS2Precursor(allmz, allint, allim, imCenterIdx=450, imWidth=20, base_int=base_int, imWinLim=(350, 600), baseFragmentMz=100.01)

    # peaks of a precursor at 422.5 m/z : 100.02, 101.02, 102.02, .. 100.02 + NR_PEAKS
    # and ion mobility 150

    createMS2Precursor(allmz, allint, allim, imCenterIdx=550, imWidth=10, base_int=base_int, imWinLim=(350, 600), baseFragmentMz=100.02)

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(len(allmz))

    for k,val in enumerate(allim):
        fda[k] = val

    

    
    ### create spectrum and set metaData
    sframe = pyopenms.MSSpectrum()
    sframe.setMSLevel(2)
    sframe.setRT(rt_idx + rt_offset)

    ## set IM window boundaries # precursor m/z 417.5 is in this window but not centered, thus should not be extracted from this window
    sframe.setMetaValue('ion mobility lower limit', 350/500)
    sframe.setMetaValue('ion mobility upper limit', 600/500)

    p = pyopenms.Precursor()
    center = 412.5
    width = 25
    p.setMZ(center)
    p.setIsolationWindowUpperOffset(width / 2.0)
    p.setIsolationWindowLowerOffset(width / 2.0)
    sframe.setPrecursors([p])

    # store data in spectrum
    sframe.set_peaks( (allmz, allint) )
    sframe.setFloatDataArrays([fda])


    sframe.sortByPosition()

    exp.addSpectrum(sframe)



############# EXPOERT MZML FILE ##########################

f = MzMLFile()
pf = f.getOptions()
pf.setCompression(True)
f.setOptions(pf)
exp.sortSpectra()
f.store('OpenSwathWorkflow_23_input.mzML', exp)
