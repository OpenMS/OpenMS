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


#ifndef OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORINGBASE_H
#define OPENMS_ANALYSIS_DENOVO_COMPNOVOIONSCORINGBASE_H

// OpenMS includes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

// stl includes
#include <vector>

namespace OpenMS
{
  /**
    @brief  run with CompNovoIonScoringBase

      @ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI CompNovoIonScoringBase :
    public DefaultParamHandler
  {

public:

    enum IsotopeType
    {
      PARENT = 0,
      CHILD = 1,
      LONE = 2
    };

    struct OPENMS_DLLAPI IonScore
    {
      IonScore();

      IonScore(const IonScore & rhs);

      virtual ~IonScore();

      IonScore & operator=(const IonScore & rhs);


      double score;
      double s_bion;
      double s_yion;
      double s_witness;
      double position;
      double s_isotope_pattern_1;   // isotope pattern score charge 1
      Int is_isotope_1_mono;   // 0 means not tested, 1 mean is, -1 is tail of isotopes
      double s_isotope_pattern_2;   // "" charge 2
    };


    /** @name constructors and destructors
     */
    //@{
    /// default constructor
    CompNovoIonScoringBase();

    /// copy constructor
    CompNovoIonScoringBase(const CompNovoIonScoringBase & source);

    /// destructor
    ~CompNovoIonScoringBase() override;
    //@}

    ///
    CompNovoIonScoringBase & operator=(const CompNovoIonScoringBase & source);

    /** @name Accessors
     */
    //@{
    double scoreIsotopes(const PeakSpectrum & CID_spec, PeakSpectrum::ConstIterator it, Size charge);
    //@}

protected:

    /// update members method from DefaultParamHandler to update the members
    void updateMembers_() override;


    IsotopeType classifyIsotopes_(const PeakSpectrum & spec, PeakSpectrum::ConstIterator it);

    double scoreIsotopes_(const PeakSpectrum & spec, PeakSpectrum::ConstIterator it, Map<double, IonScore> & CID_nodes, Size charge = 1);

    virtual void scoreWitnessSet_(Size charge, double precursor_weight, Map<double, IonScore> & CID_nodes, const PeakSpectrum & CID_orig_spec) = 0;

    void addSingleChargedIons_(Map<double, IonScore> & ion_scores, PeakSpectrum & CID_spec);

    void initIsotopeDistributions_();

    ///
    Map<Size, std::vector<double> > isotope_distributions_;

    double fragment_mass_tolerance_;

public:

  };

}

#endif
