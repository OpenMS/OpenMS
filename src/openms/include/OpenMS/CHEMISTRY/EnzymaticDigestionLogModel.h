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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow, Xiao Liang $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ENZYMATICDIGESTIONLOGMODEL_H
#define OPENMS_CHEMISTRY_ENZYMATICDIGESTIONLOGMODEL_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
       @brief Class for the Log L model of enzymatic digestion of proteins

   An alternative model for tryptic digestion. where the protein is cleaved only at positions where a cleavage model
   trained on real data, exceeds a certain threshold. The model is published in
   Siepen et al. (2007), "Prediction of missed cleavage sites in tryptic peptides aids protein identification in proteomics.", doi: 10.1021/pr060507u
   The model is only available for trypsin and ignores the missed cleavage setting. You should however use setLogThreshold()
   to adjust FP vs FN rates. A higher threshold increases the number of cleavages predicted.

       @ingroup Chemistry
  */
  class OPENMS_DLLAPI EnzymaticDigestionLogModel
  {
public:
    /// Default constructor
    EnzymaticDigestionLogModel();

    /// Copy constructor
    EnzymaticDigestionLogModel(const EnzymaticDigestionLogModel& rhs);

    /// Assignment operator
    EnzymaticDigestionLogModel& operator=(const EnzymaticDigestionLogModel& rhs);

    /// Returns the enzyme for the digestion
    String getEnzymeName() const;

    /// Sets the enzyme for the digestion
    void setEnzyme(const String name);

    /// Performs the enzymatic digestion of a protein.
    void digest(const AASequence& protein, std::vector<AASequence>& output) const;

    /// Returns the number of peptides a digestion of @p protein would yield under the current enzyme and missed cleavage settings.
    Size peptideCount(const AASequence& protein);

    /// Returns the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data)
    double getLogThreshold() const;

    /// Sets the threshold which needs to be exceeded to call a cleavage (only for the trained cleavage model on real data)
    /// Default is 0.25
    void setLogThreshold(double threshold);

protected:
    // define a binding site by position and AA
    struct BindingSite_
    {
      Size position;
      String AAname;

      BindingSite_() :
        position(), AAname() {}

      BindingSite_(const Size& p, const String& name) :
        position(p), AAname(name) {}

      bool operator<(const BindingSite_& rhs) const
      {
        return (position < rhs.position) || ((position == rhs.position) && (AAname < rhs.AAname));
      }

      bool operator==(const BindingSite_& rhs) const
      {
        return position == rhs.position && AAname == rhs.AAname;
      }

    };

    // define the log likelihood for missed and cleavage model
    struct CleavageModel_
    {
      double p_cleave;
      double p_miss;

      CleavageModel_() :
        p_cleave(0), p_miss(0) {}
      CleavageModel_(const double& p_c, const double& p_m) :
        p_cleave(p_c), p_miss(p_m) {}
    };

    /// Moves the iterator @p p behind (i.e., C-term) the next cleavage site of the @p sequence
    void nextCleavageSite_(const AASequence& sequence, AASequence::ConstIterator& p) const;

    /// Tests if position pointed to by @p p (N-term side) is a valid cleavage site
    bool isCleavageSite_(const AASequence& sequence, const AASequence::ConstIterator& p) const;

    /// Used enzyme
    const DigestionEnzyme* enzyme_;

    /// Threshold to decide if position is cleaved or missed (only for the model)
    double log_model_threshold_;
    /// Holds the cleavage model
    Map<BindingSite_, CleavageModel_> model_data_;
  };

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ENZYMATICDIGESTIONLOGMODEL_H
