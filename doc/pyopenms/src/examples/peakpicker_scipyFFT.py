"""
--------------------------------------------------------------------------
                  OpenMS -- Open-Source Mass Spectrometry
--------------------------------------------------------------------------
Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
ETH Zurich, and Freie Universitaet Berlin 2002-2013.

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
"""
###########################################################################
# Create peak-picking figure
###########################################################################

import pylab
import pyopenms
import scipy

filename = "data/NoiseFilterSGolay_1_input.mzML"


def fft_filter_spec(spec, fft_low_cutoff, fft_high_cutoff):
    # get raw data
    peaks = spec.get_peaks()
    mz_values, intensities = peaks.T

    # fft filter on intensities
    fft = scipy.fft(intensities)
    fft[: fft_low_cutoff + 1] = 0
    fft[fft_high_cutoff:] = 0
    filtered_intensities = scipy.ifft(fft).real
    peaks[:, 1] = filtered_intensities

    # build pyopenms.MSSpectrum from filtered raw data
    new_spec = pyopenms.MSSpectrum()
    new_spec.set_peaks(peaks)
    return new_spec


def fft_pick(spec, fft_low_cutoff, fft_high_cutoff):
    # filter spec
    filtered_spec = fft_filter_spec(spec, fft_low_cutoff, fft_high_cutoff)
    # apply pyopenms peakpicker on filtered spec
    picker = pyopenms.PeakPickerHiRes()
    newspec_out = pyopenms.MSSpectrum()
    picker.pick(filtered_spec, newspec_out)
    return newspec_out


def plot_spec_and_picked_spec(spectrum, picked_spectrum):
    original_peaks = spectrum.get_peaks()
    mz_values, intensities = original_peaks.T

    # plots line:
    pylab.plot(mz_values, intensities)

    pylab.xlabel("m/z")
    pylab.ylabel("Intensity")

    # scale y-axis max by 1.1
    axis = pylab.gca()
    ymin, ymax = axis.get_ylim()
    axis.set_ylim([0, ymax * 1.1])

    picked_peaks = picked_spectrum.get_peaks()
    picked_mz_values, picked_intensities = picked_peaks.T
    # plots red dots:
    pylab.plot(picked_mz_values, picked_intensities, "ro")


mse = pyopenms.MSExperiment()
fh = pyopenms.FileHandler()
fh.loadExperiment(filename, mse)

spec = mse[0]
picked_spectrum = fft_pick(spec, -1, 25)
plot_spec_and_picked_spec(spec, picked_spectrum)

pylab.savefig("plot.png")
pylab.show()
