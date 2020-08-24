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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

using namespace std;


namespace OpenMS
{
  void IDConflictResolverAlgorithm::resolve(FeatureMap & features)
  {
    resolveConflict_(features);
  }
  
  void IDConflictResolverAlgorithm::resolve(ConsensusMap & features)
  {
    resolveConflict_(features);
  }
  
  void IDConflictResolverAlgorithm::resolveBetweenFeatures(FeatureMap & features)
  {
    resolveBetweenFeatures_(features);
  }
  
  void IDConflictResolverAlgorithm::resolveBetweenFeatures(ConsensusMap & features)
  {
    resolveBetweenFeatures_(features);
  }
  
  // static
  void IDConflictResolverAlgorithm::resolveConflict_(
    vector<PeptideIdentification> & peptides, 
    vector<PeptideIdentification> & removed,
    UInt64 uid)
  {
    if (peptides.empty()) { return; }

    for (PeptideIdentification & pep : peptides)
    {
      // sort hits
      pep.sort();

      // remove all but the best hit
      if (!pep.getHits().empty())
      {
        vector<PeptideHit> best_hit(1, pep.getHits()[0]);
        pep.setHits(best_hit);
      }
      // annotate feature id
      pep.setMetaValue("feature_id", String(uid));
    }

    vector<PeptideIdentification>::iterator pos;
    if (peptides[0].isHigherScoreBetter())     // find highest-scoring ID
    {
      pos = max_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }
    else  // find lowest-scoring ID
    {
      pos = min_element(peptides.begin(), peptides.end(), compareIDsSmallerScores_);
    }

    // copy conflicting ones left to best one
    for (auto it = peptides.begin(); it != pos; ++it)
    {
      removed.push_back(*it);
    }
     
    // copy conflicting ones right of best one
    vector<PeptideIdentification>::iterator pos1p = pos + 1;
    for (auto it = pos1p; it != peptides.end(); ++it) // OMS_CODING_TEST_EXCLUDE
    {
      removed.push_back(*it);
    }

    // set best one to first position and shrink vector
    peptides[0] = *pos;
    peptides.resize(1);
  }

  // static
  bool IDConflictResolverAlgorithm::compareIDsSmallerScores_(const PeptideIdentification & left, const PeptideIdentification & right)
  {
    // if any of them is empty, the other is considered "greater"
    // independent of the score in the first hit
    if (left.getHits().empty() || right.getHits().empty()) 
    { // also: for strict weak ordering, comp(x,x) needs to be false
      return left.getHits().size() < right.getHits().size();
    }

    if (left.getHits()[0].getScore() < right.getHits()[0].getScore())
    {
      return true;
    }
    return false;
  }
}

/// @endcond
