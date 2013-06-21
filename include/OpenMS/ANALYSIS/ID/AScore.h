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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef  OPENMS_ANALYSIS_ID_ASCORE_H
#define  OPENMS_ANALYSIS_ID_ASCORE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  class PeptideHit;
  class AASequence;
  struct ProbablePhosphoSites
  {
    Size first;
    Size second;
    Size seq_1;
    Size seq_2;
    Size peak_depth;
    Size AScore;
  };
  /**
      @brief Implementation of the Ascore
      For a given peptide sequence and its MS^2 spectrum it is tried to identify the most probable phosphorylation-site(s).
      For each phosphorylation site a score is calculated, which is an indicator for the probability that this site is phosphorylated.
      The algorithm is implemented according to Beausoleil et al.

  */
  class OPENMS_DLLAPI AScore
  {
public:
    ///Default constructor
    AScore();
    ///Destructor
    ~AScore();

    /**
        @brief Computes the AScore and returns all computed phospho-sites. The saved sequences contain only phospho informations. All other modifications are dropped due to simplicity.

        @param	hit a PeptideHit
        @param real_spectrum spectrum mapped to hit
        @param fmt fragment_mass_tolerance, when comparing real_spectrum to a theoretical spectrum of the amino acid sequence of hit.
        @param number_of_phospho_sites which directs the method to search for this number of phosphorylated sites.

        @note the original sequence is saved in the PeptideHits as MetaValue Search_engine_sequence.
    */
    PeptideHit compute(PeptideHit & hit, RichPeakSpectrum & real_spectrum, DoubleReal fmt, Int number_of_phospho_sites);

    ///Computes the cumulative binomial probabilities.
    DoubleReal computeCumulativeScore(UInt N, UInt n, DoubleReal p);

    /**
        @brief Finds the peptides with the highest PeptideScores and outputs all informations for computing the AScore
        @note This function assumes that there are more permutations than the assumed number of phosphorylations!
    */
    void computeHighestPeptides(std::vector<std::vector<DoubleReal> > & peptide_site_scores, std::vector<ProbablePhosphoSites> & sites, std::vector<std::vector<Size> > & permutations);
    ///Computes the site determining_ions for the given AS and sequences in candidates
    void compute_site_determining_ions(std::vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, Int charge, std::vector<RichPeakSpectrum> & site_determining_ions);
private:
    ///computes number of matched ions between windows and the given spectrum. All spectra have to be sorted by Position!
    Int numberOfMatchedIons_(const RichPeakSpectrum & th, const RichPeakSpectrum & windows, Size depth, DoubleReal fmt);
    ///computes the peptide score according to Beausoleil et al. page 1291
    DoubleReal peptideScore_(std::vector<DoubleReal> & scores);
public:
    ///helperfunction
    std::vector<Size> computeTupel_(AASequence & without_phospho);
    ///helper function
    std::vector<std::vector<Size> > computePermutations_(std::vector<Size> tupel, Int number_of_phospho_sites);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_ASCORE_H
