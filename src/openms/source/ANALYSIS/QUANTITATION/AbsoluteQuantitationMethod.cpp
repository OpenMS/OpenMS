// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

#include <tuple>

namespace OpenMS
{

  bool AbsoluteQuantitationMethod::operator==(const AbsoluteQuantitationMethod& other) const
  {
    return
      std::tie(
        component_name_,
        feature_name_,
        IS_name_,
        llod_,
        ulod_,
        lloq_,
        uloq_,
        n_points_,
        correlation_coefficient_,
        concentration_units_,
        transformation_model_,
        transformation_model_params_
      ) == std::tie(
        other.component_name_,
        other.feature_name_,
        other.IS_name_,
        other.llod_,
        other.ulod_,
        other.lloq_,
        other.uloq_,
        other.n_points_,
        other.correlation_coefficient_,
        other.concentration_units_,
        other.transformation_model_,
        other.transformation_model_params_
      );
  }

  bool AbsoluteQuantitationMethod::operator!=(const AbsoluteQuantitationMethod& other) const
  {
    return !(*this == other);
  }
  void AbsoluteQuantitationMethod::setLLOD(const double llod)
  {
    llod_ = llod;
  }

  void AbsoluteQuantitationMethod::setULOD(const double ulod)
  {
    ulod_ = ulod;
  }

  double AbsoluteQuantitationMethod::getLLOD() const
  {
    return llod_;
  }

  double AbsoluteQuantitationMethod::getULOD() const
  {
    return ulod_;
  }

  void AbsoluteQuantitationMethod::setLLOQ(const double lloq)
  {
    lloq_ = lloq;
  }

  void AbsoluteQuantitationMethod::setULOQ(const double uloq)
  {
    uloq_ = uloq;
  }

  double AbsoluteQuantitationMethod::getLLOQ() const
  {
    return lloq_;
  }

  double AbsoluteQuantitationMethod::getULOQ() const
  {
    return uloq_;
  }

  void AbsoluteQuantitationMethod::setFeatureName(const String& feature_name)
  {
    feature_name_ = feature_name;
  }

  String AbsoluteQuantitationMethod::getFeatureName() const
  {
    return feature_name_;
  } 

  void AbsoluteQuantitationMethod::setISName(const String& IS_name)
  {
    IS_name_ = IS_name;
  }

  String AbsoluteQuantitationMethod::getISName() const
  {
    return IS_name_;
  }

  void AbsoluteQuantitationMethod::setComponentName(const String& component_name)
  {
    component_name_ = component_name;
  }

  String AbsoluteQuantitationMethod::getComponentName() const
  {
    return component_name_;
  }

  void AbsoluteQuantitationMethod::setConcentrationUnits(const String& concentration_units)
  {
    concentration_units_ = concentration_units;
  }

  String AbsoluteQuantitationMethod::getConcentrationUnits() const
  {
    return concentration_units_;
  }

  void AbsoluteQuantitationMethod::setTransformationModel(const String& transformation_model)
  {
    transformation_model_ = transformation_model;
  }

  void AbsoluteQuantitationMethod::setTransformationModelParams(const Param& transformation_model_params)
  {
    transformation_model_params_ = transformation_model_params;
  }

  String AbsoluteQuantitationMethod::getTransformationModel() const
  {
    return transformation_model_;
  }

  Param AbsoluteQuantitationMethod::getTransformationModelParams() const
  {
    return transformation_model_params_;
  }

  void AbsoluteQuantitationMethod::setNPoints(const Int n_points)
  {
    n_points_ = n_points;
  }

  void AbsoluteQuantitationMethod::setCorrelationCoefficient(const double correlation_coefficient)
  {
    correlation_coefficient_ = correlation_coefficient;
  }

  Int AbsoluteQuantitationMethod::getNPoints() const
  {
    return n_points_;
  }

  double AbsoluteQuantitationMethod::getCorrelationCoefficient() const
  {
    return correlation_coefficient_;
  }

  bool AbsoluteQuantitationMethod::checkLOD(const double value) const
  {
    return value >= llod_ && value <= ulod_; // is it bracketed or not
  }

  bool AbsoluteQuantitationMethod::checkLOQ(const double value) const
  {
    return value >= lloq_ && value <= uloq_; // is it bracketed or not
  }
} // namespace
