// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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


#pragma once

// OpenMS includes
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>

// stl includes
#include <vector>

namespace OpenMS
{
  /**
    @brief  run with CompNovoIdentification

      @htmlinclude OpenMS_CompNovoIdentification.parameters

      @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI CompNovoIdentification :
    public CompNovoIdentificationBase
  {

public:

    /** @name constructors and destructors
     */
    //@{
    /// default constructor
    CompNovoIdentification();

    /// copy constructor
    CompNovoIdentification(const CompNovoIdentification & source);

    /// destructor
    ~CompNovoIdentification() override;
    //@}

    ///
    CompNovoIdentification & operator=(const CompNovoIdentification & source);

    /** @name Accessors
     */
    //@{
    /// performs an ProteinIdentification run on a PeakMap
    void getIdentifications(std::vector<PeptideIdentification> & ids, const PeakMap & exp) override;

    /// performs an ProteinIdentification run on a PeakSpectrum
    void getIdentification(PeptideIdentification & id, const PeakSpectrum & CID_spec, const PeakSpectrum & ETD_spec);
    //@}


    typedef CompNovoIonScoringBase::IsotopeType IsotopeType;
    typedef CompNovoIonScoringBase::IonScore IonScore;
    typedef CompNovoIdentificationBase::Permut Permut;

protected:

    /// call the DAC algorithm for the subspectrum defined via left and right peaks and fill the set with candidates sequences
    void getDecompositionsDAC_(std::set<String> & sequences, Size left, Size right, double peptide_weight, const PeakSpectrum & CID_orig_spec, const PeakSpectrum & ETD_orig_spec, Map<double, IonScore> & CID_nodes);

    /// reduces the given number of permuts by scoring the permutations to the CID and ETD spec
    void reducePermuts_(std::set<String> & permuts, const PeakSpectrum & CID_orig_spec, const PeakSpectrum & ETD_orig_spec, double prefix, double suffix);

    /// estimates an exact precursor weight of the ETD spectrum, because in most of the cases the precursor is found in the MS/MS spec
    double estimatePrecursorWeight_(const PeakSpectrum & ETD_spec, Size & charge);

  };
}

