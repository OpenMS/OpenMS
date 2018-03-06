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
// $Authors: Sandro Andreotti, Andreas Bertsch$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H
#define OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for ion scoring implementation for de novo algorithms



        @ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI DeNovoIonScoring :
    public DefaultParamHandler
  {
public:
    /** @brief IonScore

                  IonScore describes the likelihood of position to be prefix residue masses (in case of
                  an b-ion centric algorithm).
          */
    class IonScore
    {
public:
      /// score
      double score;

      /// position of the ion
      double position;

      /// index of peak in the spectrum, -1 if not in spectrum
      SignedSize index;


      IonScore() :
        score(0),
        position(0.0),
        index(-1)
      {
      }

      IonScore(const IonScore & rhs) :
        score(rhs.score),
        position(rhs.score),
        index(rhs.index)
      {
      }

      virtual ~IonScore()
      {
      }

      IonScore & operator=(const IonScore & rhs)
      {
        if (this != &rhs)
        {
          score = rhs.score;
          position = rhs.position;
          index = rhs.index;
        }
        return *this;
      }

    };


    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    DeNovoIonScoring();

    /// destructor
    ~DeNovoIonScoring() override;

    /// copy constructor
    DeNovoIonScoring(const DeNovoIonScoring & rhs);
    //@}

    /// assignment operator
    DeNovoIonScoring & operator=(const DeNovoIonScoring & rhs);

    virtual void getIonScores(std::vector<IonScore> & ion_scores, const PeakSpectrum & spec) = 0;

    virtual void getIonScores(std::vector<std::vector<IonScore> > & ion_scores, const PeakMap & exp) = 0;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DENOVO_DENOVOIONSCORING_H
