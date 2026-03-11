// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MRMFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMChromHandler.h>
#include <cstdio>

namespace OpenMS
{

std::vector<::OpenSwath::SwathMap> MRMFile::loadMzML(const String& file,
                                                   const String& /*tmp*/, 
                                                   std::shared_ptr<ExperimentalSettings>& exp_meta)
{
  OPENMS_LOG_INFO << "Loading SRM/MRM mzML " << file << std::endl;
  std::shared_ptr<PeakMap> full_exp(new PeakMap);
  FileHandler fh;
  // load full chromatograms
  fh.loadExperiment(file, *full_exp, {FileTypes::MZML});

  // Normalize / canonicalize chromatogram precursor/product m/z and nativeIDs
  for (Size i = 0; i < full_exp->getChromatograms().size(); ++i)
  {
    MSChromatogram & chrom = full_exp->getChromatograms()[i];
    MRMChromHandler::normalizeChromatogramMZ(chrom);
  }

  exp_meta = full_exp;
  std::vector<::OpenSwath::SwathMap> swath_maps;
  ::OpenSwath::SwathMap sm;
  sm.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(full_exp);
  sm.lower = -1;
  sm.upper = -1;
  sm.center = -1;
  sm.imLower = -1;
  sm.imUpper = -1;
  sm.ms1 = false;
  swath_maps.push_back(sm);
  return swath_maps;
}

} // namespace OpenMS
