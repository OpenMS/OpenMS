// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/DataProcessing.h>

using namespace std;

namespace OpenMS
{
  const std::string DataProcessing::NamesOfProcessingAction[] =
  {
    "Data processing action",
    "Charge deconvolution",
    "Deisotoping",
    "Smoothing",
    "Charge calculation",
    "Precursor recalculation",
    "Baseline reduction",
    "Peak picking",
    "Retention time alignment",
    "Calibration of m/z positions",
    "Intensity normalization",
    "Data filtering",
    "Quantitation",
    "Feature grouping",
    "Identification mapping",
    "File format conversion",
    "Conversion to mzData format",
    "Conversion to mzML format",
    "Conversion to mzXML format",
    "Conversion to DTA format",
    "Identification",
    "Ion mobility binning"
  };

  DataProcessing::~DataProcessing() = default;

  DataProcessing::DataProcessing(DataProcessing&& rhs) noexcept :
    MetaInfoInterface(std::move(rhs)),
    software_(std::move(rhs.software_)),
    processing_actions_(std::move(rhs.processing_actions_)),
    completion_time_(std::move(rhs.completion_time_))
  {
  }

  bool DataProcessing::operator==(const DataProcessing& rhs) const
  {
    return software_ == rhs.software_ &&
           processing_actions_ == rhs.processing_actions_ &&
           completion_time_ == rhs.completion_time_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool DataProcessing::operator!=(const DataProcessing& rhs) const
  {
    return !(operator==(rhs));
  }

  const Software& DataProcessing::getSoftware() const
  {
    return software_;
  }

  Software& DataProcessing::getSoftware()
  {
    return software_;
  }

  void DataProcessing::setSoftware(const Software& software)
  {
    software_ = software;
  }

  const DateTime& DataProcessing::getCompletionTime() const
  {
    return completion_time_;
  }

  void DataProcessing::setCompletionTime(const DateTime& completion_time)
  {
    completion_time_ = completion_time;
  }

  const set<DataProcessing::ProcessingAction>& DataProcessing::getProcessingActions() const
  {
    return processing_actions_;
  }

  set<DataProcessing::ProcessingAction>& DataProcessing::getProcessingActions()
  {
    return processing_actions_;
  }

  void DataProcessing::setProcessingActions(const set<DataProcessing::ProcessingAction>& processing_actions)
  {
    processing_actions_ = processing_actions;
  }

}
