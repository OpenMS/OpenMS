#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright OpenMS Inc. -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-present.

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
$Maintainer: Hannes Roest $
$Authors: Hannes Roest $
--------------------------------------------------------------------------
"""
from __future__ import division
from __future__ import print_function

# Create simulated ion mobility scans for testing

from pyopenms import *

exp = MSExperiment()
# print(dir(exp))

NR_RT_SAMPLES = 50
NR_IM_BINS = 300
NR_PEAKS = 5

# Create MS1 spectra
for rt_idx in range(NR_RT_SAMPLES):
    sp = MSSpectrum()
    sp.setRT(rt_idx)
    sp.setMSLevel(1)

    base_int = 100 - abs(20 - rt_idx)*(100)/25.0  # shift the precursor slightly
    base_int = max(5, base_int)

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(100)

    for i in range(100):
        p = Peak1D()
        p.setMZ(100+i)
        p.setIntensity(base_int+i)
        sp.push_back(p)
        fda[i] = 10

    for i in range(10):
        p = Peak1D()
        p.setMZ(412.502)
        p.setIntensity(base_int+ (i-5))
        sp.push_back(p)
        fda.push_back( 99 + (i-5) )

    for i in range(10):
        p = Peak1D()
        p.setMZ(417.502)
        p.setIntensity(base_int+ (i-5))
        sp.push_back(p)
        fda.push_back( 152 + (i-5) )

    sp.setFloatDataArrays([fda])
    exp.addSpectrum(sp)

# Create MS2 spectra
for rt_idx in range(NR_RT_SAMPLES):
    # Base intensity is a single peak at 25 RT with 100 *i intensity spread across ion mobility scans
    # Second base intensity is a single peak at 10 RT with 100 *i intensity spread across ion mobility scans
    base_int = 100 - abs(25 - rt_idx)*(100)/25.0
    base_int_second = 100 - abs(10 - rt_idx)*(100)/40.0
    print("base int", base_int, abs(25 - rt_idx)*(100)/25.0  )

    allmz = []
    allint = []
    allim = []

    for im_idx in range(NR_IM_BINS):
        sp = MSSpectrum()
        p = pyopenms.Precursor()
        p.setIsolationWindowLowerOffset(12)
        p.setIsolationWindowUpperOffset(12)
        target_mz = im_idx * 2.5 + 400
        p.setMZ(target_mz)
        sp.setPrecursors([p])

        sp.setRT(rt_idx+0.2)
        sp.setMSLevel(2)

        # peaks of a precursor at 412.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
        # and ion mobility 100
        for i in range(NR_PEAKS):
            if im_idx > 90 and im_idx < 90 + 20:
                apex_dist = abs( 100 - im_idx)
                p = Peak1D()
                p.setMZ(100+i)
                p.setIntensity(base_int * (i + 1) - base_int * (i + 1) * apex_dist / 10.0)
                allmz.append(p.getMZ())
                allint.append(p.getIntensity())
                allim.append( im_idx / 500.0)

        # peaks of a precursor at 417.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
        # and ion mobility 150
        for i in range(NR_PEAKS):
            if im_idx > 130 and im_idx < 130 + 40:
                apex_dist = abs( 150 - im_idx)
                p = Peak1D()
                p.setMZ(100+i)
                p.setIntensity(base_int * (i + 1) - base_int * (i + 1) * apex_dist / 20.0)
                allmz.append(p.getMZ())
                allint.append(p.getIntensity())
                allim.append( im_idx / 500.0)

    mz = allmz
    intens = allint
    ims = allim

    # print(mz, intens)

    fda = pyopenms.FloatDataArray()
    fda.setName("Ion Mobility")
    fda.resize(len(mz))
    for k,val in enumerate(ims):
        fda[k] = val

    sframe = pyopenms.MSSpectrum()
    sframe.setMSLevel(2)
    sframe.setRT(rt_idx+0.2)
    sframe.setFloatDataArrays([fda])
    p = pyopenms.Precursor()
    if True:
        center = 412.5
        width = 25
        p.setMZ(center)
        p.setIsolationWindowUpperOffset(width / 2.0)
        p.setIsolationWindowLowerOffset(width / 2.0)
    sframe.setPrecursors([p])
    sframe.set_peaks( (mz, intens) )
    sframe.sortByPosition()
    exp.addSpectrum(sframe)


f = MzMLFile()
pf = f.getOptions()
pf.setCompression(True)
f.setOptions(pf)
exp.sortSpectra()
f.store('output.mzML', exp)

