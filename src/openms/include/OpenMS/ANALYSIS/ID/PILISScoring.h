// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PILISSCORING_H
#define OPENMS_ANALYSIS_ID_PILISSCORING_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief This class actually implements the E-value based scoring of PILIS

      The method which is used to generate the E-values are adapted from

      David Fenyo and Ronald C. Beavis
      Anal. Chem. 2003, 75, 768-774
      A Method for Assessing the Statistical Significance of Mass Spectrometry-Based Protein
      Identifications Using General Scoring Schemes.

      The bases for the calculation are the similarity scores of the simulated
      and experimental spectra. The scores are transformed into a discrete
      score distribution and from this distribution E-values are calculated for
      the peptide hits.

      If more than one spectrum is given two E-values can be calculated, one which gives
      the significance of the peptide hit considering only one spectrum, and the other also
      considering also all other hits of all other spectra. The second type of scoring
      is somewhat more accurate.

      @htmlinclude OpenMS_PILISScoring.parameters

      @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI PILISScoring :
    public DefaultParamHandler
  {

public:

    /** @name constructors and destructors
     */
    //@{
    /// default constructor
    PILISScoring();

    /// copy constructor
    PILISScoring(const PILISScoring & source);

    /// destructor
    virtual ~PILISScoring();
    //@}

    ///
    PILISScoring & operator=(const PILISScoring & source);

    /** @name Accessors
     */
    //@{
    /// performs an ProteinIdentification run on a PeakMap
    void getScores(std::vector<PeptideIdentification> & ids);

    /// performs an ProteinIdentification run on a PeakSpectrum
    void getScore(PeptideIdentification & id);
    //@}

protected:

    ///
    void getFitParameter_(double & slope, double & intercept, const std::vector<double> & scores, double threshold);

    ///
    void getSurvivalFunction_(Map<UInt, double> & points, std::vector<DPosition<2> > & survival_function);

    ///
    void getScore_(PeptideIdentification & id, double global_slope, double global_intercept);

  };
}

#endif
