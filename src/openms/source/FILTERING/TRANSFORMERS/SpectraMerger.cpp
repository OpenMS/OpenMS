// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>

using namespace std;
namespace OpenMS
{

  SpectraMerger::SpectraMerger() :
    DefaultParamHandler("SpectraMerger")
  {
    // common
    defaults_.setValue("mz_binning_width", 5.0, "minimum m/z distance for two data points (profile data) or peaks (centroided data) to be considered distinct. Closer data points or peaks will be merged.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("mz_binning_width", 0.0);

    defaults_.setValue("mz_binning_width_unit", "ppm", "Unit in which the distance between two data points or peaks is given.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("mz_binning_width_unit", ListUtils::create<String>("Da,ppm"));

    defaults_.setValue("sort_blocks", "RT_ascending", "Sort blocks by <?> before merging them (useful for precursor order)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("sort_blocks", ListUtils::create<String>("RT_ascending, RT_descending"));

    // Gaussian average
    defaults_.setValue("average_gaussian:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
    defaults_.setValidStrings("average_gaussian:spectrum_type", ListUtils::create<String>("profile,centroid,automatic"));
    defaults_.setValue("average_gaussian:ms_level", 1, "Average spectra of this level. All other spectra remain unchanged.");
    defaults_.setMinInt("average_gaussian:ms_level", 1);
    defaults_.setValue("average_gaussian:rt_FWHM", 5.0, "FWHM of Gauss curve in seconds to be averaged over.");
    defaults_.setMinFloat("average_gaussian:rt_FWHM", 0.0);
    defaults_.setMaxFloat("average_gaussian:rt_FWHM", 10e10);
    defaults_.setValue("average_gaussian:cutoff", 0.01, "Intensity cutoff for Gaussian. The Gaussian RT profile decreases from 1 at its apex to 0 at infinity. Spectra for which the intensity of the Gaussian drops below the cutoff do not contribute to the average.", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("average_gaussian:cutoff", 0.0);
    defaults_.setMaxFloat("average_gaussian:cutoff", 1.0);

    // top-hat average
    defaults_.setValue("average_tophat:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
    defaults_.setValidStrings("average_tophat:spectrum_type", ListUtils::create<String>("profile,centroid,automatic"));
    defaults_.setValue("average_tophat:ms_level", 1, "Average spectra of this level. All other spectra remain unchanged.");
    defaults_.setMinInt("average_tophat:ms_level", 1);
    defaults_.setValue("average_tophat:rt_range", 5.0, "RT range to be averaged over, i.e. +/-(RT range)/2 from each spectrum.");
    defaults_.setMinFloat("average_tophat:rt_range", 0.0);
    defaults_.setMaxFloat("average_tophat:rt_range", 10e10);
    defaults_.setValue("average_tophat:rt_unit", "scans", "Unit for RT range.");
    defaults_.setValidStrings("average_tophat:rt_unit", ListUtils::create<String>("scans,seconds"));

    // block merging
    defaults_.setValue("block_method:ms_levels", ListUtils::create<Int>("1"), "Merge spectra of this level. All spectra with other MS levels remain untouched.");
    defaults_.setMinInt("block_method:ms_levels", 1);
    defaults_.setValue("block_method:rt_block_size", 5, "Maximum number of scans to be summed up.");
    defaults_.setMinInt("block_method:rt_block_size", 1);

    defaults_.setValue("block_method:rt_max_length", 0.0, "Maximum RT size of the block in seconds (0.0 = no size restriction).");
    defaults_.setMinFloat("block_method:rt_max_length", 0.0);
    defaults_.setMaxFloat("block_method:rt_max_length", 10e10);

    // same precursor MS/MS merging
    defaults_.setValue("precursor_method:mz_tolerance", 10e-5, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].");
    defaults_.setMinFloat("precursor_method:mz_tolerance", 0);
    defaults_.setValue("precursor_method:rt_tolerance", 5.0, "Max RT distance of the precursor entries of two spectra to be merged in [s].");
    defaults_.setMinFloat("precursor_method:rt_tolerance", 0);

    defaultsToParam_();
  }

  SpectraMerger::SpectraMerger(const SpectraMerger & source) :
    DefaultParamHandler(source), ProgressLogger() //we probably want a new ProgressLogger when we copy
  {
  }

  SpectraMerger::~SpectraMerger()
  {
  }

  SpectraMerger & SpectraMerger::operator=(const SpectraMerger & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
