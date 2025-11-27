// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionHelper.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

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
                                     LightCompound& pep)
  {
    std::vector<LightCompound>::const_iterator beg = lte.compounds.begin();
    std::vector<LightCompound>::const_iterator end = lte.compounds.end();
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
