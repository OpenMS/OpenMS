// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IonDetector.h>

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

}

