// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Swenja Wagner, Patricia Scheil $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <map>
#include <vector>
namespace OpenMS
{
  class FeatureMap;
  /**
   * @brief This class is a metric for the QualityControl TOPP Tool.
   *
   * This class counts the number of MissedCleavages per PeptideIdentification given a FeatureMap
   * and returns an agglomeration statistic (observed counts).
   * Additionally the PeptideHits in the FeatureMap are augmented with MetaInformation:
   *  - 'missed_cleavages'
   *  - 'FWHM' (from feature's 'FWHM' or 'model_FWHM')
   *  - 'mass' (experimental mass of peptide)
   */
  class OPENMS_DLLAPI MissedCleavages : public QCBase
  {
  private:
    typedef std::map<UInt32, UInt32> MapU32;
    /// collects number of missed cleavages from PeptideIdentification in a result map (missed cleavages: occurences)
    void get_missed_cleavages_from_peptide_identification_(const ProteaseDigestion& digestor, MapU32& result, const UInt32& max_mc, PeptideIdentification& pep_id);

  public:
    /// constructor
    MissedCleavages() = default;

    /// destructor
    virtual ~MissedCleavages() = default;

    /**
     * @brief Counts the number of missed cleavages per PeptideIdentification.
     *
     * The result is a key/value map: #missed_cleavages --> counts
     * Additionally the first PeptideHit in each PeptideIdentification of the FeatureMap is annotated with metavalue 'missed_cleavages'.
     * The protease and digestion parameters are taken from the first ProteinIdentication (and SearchParameter therein) within the FeatureMap itself.
     *
     * @param fmap FeatureMap with Peptide and ProteinIdentifications
     */
    void compute(FeatureMap& fmap);
    void compute(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids);

    /// returns the name of the metric
    const String& getName() const override;

    /// returns the result as maps of number of missed_cleavages to counts; one map for each call to compute(...)
    const std::vector<std::map<UInt32, UInt32>>& getResults() const;

    /**
     * @brief Returns the input data requirements of the compute(...) function
     * @return Status for POSTFDRFEAT;
     */
    QCBase::Status requirements() const override;

  private:
    /// container that stores results
    std::vector<std::map<UInt32, UInt32>> mc_result_;
  };
} // namespace OpenMS
