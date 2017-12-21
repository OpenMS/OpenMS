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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORING_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORING_H

// OpenMS includes
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>

// stl includes
#include <vector>

namespace OpenMS
{
  /**
    @brief  run with CompNovoIonScoring

      @htmlinclude OpenMS_CompNovoIonScoring.parameters

      @ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI CompNovoIonScoring :
    public CompNovoIonScoringBase
  {

public:

    typedef CompNovoIonScoringBase::IsotopeType IsotopeType;
    typedef CompNovoIonScoringBase::IonScore IonScore;


    /** @name constructors and destructors
     */
    //@{
    /// default constructor
    CompNovoIonScoring();

    /// copy constructor
    CompNovoIonScoring(const CompNovoIonScoring & source);

    /// destructor
    ~CompNovoIonScoring() override;
    //@}

    /// assignment operator
    CompNovoIonScoring & operator=(const CompNovoIonScoring & source);

    /** @name Accessors
     */
    //@{
    void scoreSpectra(Map<double, IonScore> & CID_ion_scores, PeakSpectrum & CID_spec, PeakSpectrum & ETD_spec, double precursor_weight, Size charge);
    //@}

protected:

    void scoreETDFeatures_(Size charge, double precursor_weight, Map<double, IonScore> & CID_nodes, const PeakSpectrum & CID_orig_spec, const PeakSpectrum & ETD_orig_spec);

    void scoreWitnessSet_(Size charge, double precursor_weight, Map<double, IonScore> & CID_nodes, const PeakSpectrum & CID_orig_spec) override;

  };

}

#endif
