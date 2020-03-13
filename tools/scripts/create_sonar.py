#!/usr/bin/env python
# -*- coding: utf-8  -*-
"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-2020.

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

# Create simulated sonar scans for testing

from pyopenms import *

exp = MSExperiment()
# print(dir(exp))

# Create MS1 spectra
for rt_idx in range(50):
    sp = MSSpectrum()
    sp.setRT(rt_idx)
    sp.setMSLevel(1)
    for i in range(100):
        p = Peak1D()
        p.setMZ(100+i)
        p.setIntensity(100+i)
        sp.push_back(p)
    exp.addSpectrum(sp)


NR_RT_SAMPLES = 50
NR_SONAR_SP = 200
NR_SONAR_SP = 50
NR_PEAKS = 5

# Create MS2 spectra
for rt_idx in range(NR_RT_SAMPLES):
    # Base intensity is a single peak at 25 RT with 100 *i intensity spread across 9 SONAR scans
    # Second base intensity is a single peak at 10 RT with 100 *i intensity spread across 9 SONAR scans
    base_int = 100 - abs(25 - rt_idx)*(100)/25.0
    base_int_second = 100 - abs(10 - rt_idx)*(100)/40.0
    print("base int", base_int, abs(25 - rt_idx)*(100)/25.0  )

    for sonar_idx in range(NR_SONAR_SP):
        print("======================================= sonar", sonar_idx)
        sp = MSSpectrum()
        p = pyopenms.Precursor()
        p.setIsolationWindowLowerOffset(12)
        p.setIsolationWindowUpperOffset(12)
        target_mz = sonar_idx * 2.5 + 400
        p.setMZ(target_mz)
        sp.setPrecursors([p])

        sp.setRT(rt_idx)
        sp.setMSLevel(2)

        # peaks of a precursor at 412.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
        #  range from window 0 to 10
        for i in range(NR_PEAKS):
            if 412.5 > target_mz - 12 and 412.5 < target_mz + 12:
                p = Peak1D()
                p.setMZ(100+i)
                p.setIntensity(base_int * (i + 1) + sonar_idx)
                sp.push_back(p)
            else:
                # add noise data (6x less)
                p = Peak1D()
                p.setMZ(100+i)
                p.setIntensity(base_int * (i + 1)/ 6.0)
                sp.push_back(p)

        # peaks of a precursor at 462.5 m/z : 100, 101, 102, .. 100 + NR_PEAKS
        #  range from window 20 to 30
        for i in range(NR_PEAKS):
            if 462.5 > target_mz - 12 and 462.5 < target_mz + 12:
                p = Peak1D()
                p.setMZ(300+i)
                p.setIntensity(base_int_second * (i + 1))
                sp.push_back(p)
            else:
                # add noise data (6x less)
                p = Peak1D()
                p.setMZ(300+i)
                p.setIntensity(base_int_second * (i + 1)/ 6.0)
                sp.push_back(p)

        exp.addSpectrum(sp)
    # For debug:
    # break


f = MzMLFile()
pf = f.getOptions()
pf.setCompression(True)
f.setOptions(pf)
f.store('sonar.mzML', exp)


