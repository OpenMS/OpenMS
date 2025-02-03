// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest,  Witold Wolski$
// $Authors:  Hannes Roest, Witold Wolski$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>

namespace OpenMS
{

  bool SimpleOpenMSSpectraFactory::isExperimentCached(const boost::shared_ptr<PeakMap>& exp)
  {
    for (std::size_t i = 0; i < exp->getSpectra().size(); ++i)
    {
      for (std::size_t j = 0; j < exp->getSpectra()[i].getDataProcessing().size(); j++)
      {
        if (exp->getSpectra()[i].getDataProcessing()[j]->metaValueExists("cached_data"))
        {
          return true;
        }
      }
    }
    for (std::size_t i = 0; i < exp->getChromatograms().size(); ++i)
    {
      for (std::size_t j = 0; j < exp->getChromatograms()[i].getDataProcessing().size(); j++)
      {
        if (exp->getChromatograms()[i].getDataProcessing()[j]->metaValueExists("cached_data"))
        {
          return true;
        }
      }
      }
    return false;
  }

  OpenSwath::SpectrumAccessPtr SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(const boost::shared_ptr<PeakMap>& exp)
  {
    bool is_cached = SimpleOpenMSSpectraFactory::isExperimentCached(exp);
    if (is_cached)
    {
      OpenSwath::SpectrumAccessPtr experiment(new OpenMS::SpectrumAccessOpenMSCached(exp->getLoadedFilePath()));
      return experiment;
    }
    else
    {
      OpenSwath::SpectrumAccessPtr experiment(new OpenMS::SpectrumAccessOpenMS(exp));
      return experiment;
    }
  }

}//end Namespace

