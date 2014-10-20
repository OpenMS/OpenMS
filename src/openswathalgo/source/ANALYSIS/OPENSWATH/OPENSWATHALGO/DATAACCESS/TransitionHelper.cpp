// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenSwath
{

  void TransitionHelper::convert(LightTargetedExperiment& lte,
                                 std::map<std::string,
                                 std::vector<OpenSwath::LightTransition> >& transmap)
  {

    typedef std::pair<std::string, std::vector<OpenSwath::LightTransition> > Mpair;
    typedef std::map<std::string, std::vector<OpenSwath::LightTransition> > Mmap;
    std::vector<LightTransition> ltrans = lte.getTransitions();

    /*iterate over transitions*/
    std::vector<LightTransition>::iterator ltrit = ltrans.begin();
    std::vector<LightTransition>::iterator ltrend = ltrans.end();
    for (; ltrit != ltrend; ++ltrit)
    {
      std::string pepref = ltrit->getPeptideRef();

      Mmap::iterator it = transmap.find(pepref);
      if (it == transmap.end())
      {
        std::vector<LightTransition> ltv;
        ltv.push_back(*ltrit);
        transmap.insert(Mpair(pepref, ltv));
      }
      else
      {
        it->second.push_back(*ltrit);
      }
    }
  } //end convert

  bool TransitionHelper::findPeptide(const LightTargetedExperiment& lte,
                                     const std::string& peptideRef,
                                     LightPeptide& pep)
  {
    std::vector<LightPeptide>::const_iterator beg = lte.peptides.begin();
    std::vector<LightPeptide>::const_iterator end = lte.peptides.end();
    for (; beg != end; ++beg)
    {
      //std::cout << beg->id << " " << peptideRef << std::endl;
      if (beg->id.compare(peptideRef) == 0)
      {
        pep = *beg;
        return true;
      }
    }
    return false;
  }

} // namespace OpenSwath
