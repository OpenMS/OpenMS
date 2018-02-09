// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DIAPRESCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_DIAPRESCORING_H

#include <algorithm>
#include <iterator>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataFrameWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief Scoring of an spectrum given library intensities of a transition group.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    See also class DIAScoring.

    Simulate theoretical spectrum from library intensities of transition group
    and compute manhattan distance and dotprod score between spectrum intensities
    and simulated spectrum.
  */

  class OPENMS_DLLAPI DiaPrescore :
    public DefaultParamHandler
  {
    double dia_extract_window_; //done
    int nr_isotopes_;
    int nr_charges_;
public:

    DiaPrescore();

    DiaPrescore(double dia_extract_window, int nr_isotopes = 4, int nr_charges = 4);

    void defineDefaults();

    void updateMembers_() override;


    /**
      @brief Score a spectrum given a transition group.

      Simulate theoretical spectrum from library intensities of transition group
      and compute manhattan distance and dotprod score between spectrum intensities
      and simulated spectrum.
    */
    void score(OpenSwath::SpectrumPtr spec,
               const std::vector<OpenSwath::LightTransition>& lt,
               double& dotprod,
               double& manhattan);

    /**
      @brief Compute manhattan and dotprod score for all spectra which can be accessed by
      the SpectrumAccessPtr for all transitions groups in the LightTargetedExperiment.
    */
    void operator()(OpenSwath::SpectrumAccessPtr swath_ptr,
                    OpenSwath::LightTargetedExperiment& transition_exp_used,
                    OpenSwath::IDataFrameWriter* ivw);
  };


}

#endif
