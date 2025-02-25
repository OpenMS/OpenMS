// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  const std::string QCBase::names_of_requires[] = {"fail", "raw.mzML", "postFDR.featureXML", "preFDR.featureXML", "contaminants.fasta", "trafoAlign.trafoXML"};

  const std::string QCBase::names_of_toleranceUnit[] = {"auto", "ppm", "da"};

  QCBase::SpectraMap::SpectraMap(const MSExperiment& exp)
  {
    calculateMap(exp);
  }

  void QCBase::SpectraMap::calculateMap(const MSExperiment& exp)
  {
    nativeid_to_index_.clear();
    for (Size i = 0; i < exp.size(); ++i)
    {
      nativeid_to_index_[exp[i].getNativeID()] = i;
    }
  }

  UInt64 QCBase::SpectraMap::at(const String& identifier) const
  {
    const auto& it = nativeid_to_index_.find(identifier);
    if (it == nativeid_to_index_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("No spectrum with identifier '") + identifier + "' in MSExperiment!");
    }
    return it->second;
  }

  void QCBase::SpectraMap::clear()
  {
    nativeid_to_index_.clear();
  }

  bool QCBase::SpectraMap::empty() const
  {
    return nativeid_to_index_.empty();
  }

  Size QCBase::SpectraMap::size() const
  {
    return nativeid_to_index_.size();
  }

  // function tests if a metric has the required input files
  // gives a warning with the name of the metric that can not be performed
  bool QCBase::isRunnable(const Status& s) const
  {
    if (s.isSuperSetOf(this->requirements()))
    {
      return true;
    }
    for (Size i = 0; i < (UInt64)QCBase::Requires::SIZE_OF_REQUIRES; ++i)
    {
      if (this->requirements().isSuperSetOf(QCBase::Requires(i)) && !s.isSuperSetOf(QCBase::Requires(i)))
      {
        OPENMS_LOG_WARN << "Note: Metric '" << this->getName() << "' cannot run because input data '" << QCBase::names_of_requires[i] << "' is missing!\n";
      }
    }
    return false;
  }

  bool QCBase::isLabeledExperiment(const ConsensusMap& cm)
  {
    bool iso_analyze = true;
    auto cm_dp = cm.getDataProcessing(); // get a copy to avoid calling .begin() and .end() on two different temporaries
    if (all_of(cm_dp.begin(), cm_dp.end(), [](const OpenMS::DataProcessing& dp) { return (dp.getSoftware().getName() != "IsobaricAnalyzer"); }))
    {
      iso_analyze = false;
    }
    return iso_analyze;
  }

} // namespace OpenMS
