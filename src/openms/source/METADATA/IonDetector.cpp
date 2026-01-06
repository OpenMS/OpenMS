// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>

using namespace std;

namespace OpenMS
{

  const std::string IonDetector::NamesOfType[] = {"Unknown", "Electron multiplier", "Photo multiplier", "Focal plane array", "Faraday cup", "Conversion dynode electron multiplier", "Conversion dynode photo multiplier", "Multi-collector", "Channel electron multiplier", "channeltron", "daly detector", "microchannel plate detector", "array detector", "conversion dynode", "dynode", "focal plane collector", "ion-to-photon detector", "point collector", "postacceleration detector", "photodiode array detector", "inductive detector", "electron multiplier tube"};

  const std::string IonDetector::NamesOfAcquisitionMode[] = {"Unknown", "Pulse counting", "Analog-digital converter", "Time-digital converter", "Transient recorder"};

  IonDetector::IonDetector() :
    MetaInfoInterface(),
    type_(TYPENULL),
    acquisition_mode_(ACQMODENULL),
    resolution_(0.0),
    ADC_sampling_frequency_(0.0),
    order_(0)
  {

  }

  IonDetector::~IonDetector() = default;

  bool IonDetector::operator==(const IonDetector & rhs) const
  {
    return order_ == rhs.order_ &&
           type_ == rhs.type_ &&
           acquisition_mode_ == rhs.acquisition_mode_ &&
           resolution_ == rhs.resolution_ &&
           ADC_sampling_frequency_ == rhs.ADC_sampling_frequency_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool IonDetector::operator!=(const IonDetector & rhs) const
  {
    return !(operator==(rhs));
  }

  IonDetector::Type IonDetector::getType() const
  {
    return type_;
  }

  void IonDetector::setType(IonDetector::Type type)
  {
    type_ = type;
  }

  IonDetector::AcquisitionMode IonDetector::getAcquisitionMode() const
  {
    return acquisition_mode_;
  }

  void IonDetector::setAcquisitionMode(IonDetector::AcquisitionMode acquisition_mode)
  {
    acquisition_mode_ = acquisition_mode;
  }

  double IonDetector::getResolution() const
  {
    return resolution_;
  }

  void IonDetector::setResolution(double resolution)
  {
    resolution_ = resolution;
  }

  double IonDetector::getADCSamplingFrequency() const
  {
    return ADC_sampling_frequency_;
  }

  void IonDetector::setADCSamplingFrequency(double ADC_sampling_frequency)
  {
    ADC_sampling_frequency_ = ADC_sampling_frequency;
  }

  Int IonDetector::getOrder() const
  {
    return order_;
  }

  void IonDetector::setOrder(Int order)
  {
    order_ = order;
  }

  StringList IonDetector::getAllNamesOfType()
  {
    StringList names;
    names.reserve(SIZE_OF_TYPE);
    for (size_t i = 0; i < SIZE_OF_TYPE; ++i)
    {
      names.push_back(NamesOfType[i]);
    }
    return names;
  }

  StringList IonDetector::getAllNamesOfAcquisitionMode()
  {
    StringList names;
    names.reserve(SIZE_OF_ACQUISITIONMODE);
    for (size_t i = 0; i < SIZE_OF_ACQUISITIONMODE; ++i)
    {
      names.push_back(NamesOfAcquisitionMode[i]);
    }
    return names;
  }

  const std::string& IonDetector::typeToString(Type type)
  {
    if (type == SIZE_OF_TYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_TYPE");
    }
    return NamesOfType[static_cast<size_t>(type)];
  }

  IonDetector::Type IonDetector::toType(const std::string& name)
  {
    auto first = &NamesOfType[0];
    auto last = &NamesOfType[SIZE_OF_TYPE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<Type>(it - first);
  }

  const std::string& IonDetector::acquisitionModeToString(AcquisitionMode mode)
  {
    if (mode == SIZE_OF_ACQUISITIONMODE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_ACQUISITIONMODE");
    }
    return NamesOfAcquisitionMode[static_cast<size_t>(mode)];
  }

  IonDetector::AcquisitionMode IonDetector::toAcquisitionMode(const std::string& name)
  {
    auto first = &NamesOfAcquisitionMode[0];
    auto last = &NamesOfAcquisitionMode[SIZE_OF_ACQUISITIONMODE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<AcquisitionMode>(it - first);
  }

}

