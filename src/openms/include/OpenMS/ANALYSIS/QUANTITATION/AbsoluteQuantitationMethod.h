// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{
  /**
    @brief AbsoluteQuantitationMethod is a class to hold information about the
    quantitation method and for applying and/or generating the quantitation method.

    The quantitation method describes all parameters required to define the
    calibration curve used for absolute quantitation by Isotope Dilution Mass Spectrometry (IDMS).
    The quantitation method also defines the statistics of the fitted calibration curve as well
    as the lower and upper bounds of the calibration for later Quality Control.
  */
  class OPENMS_DLLAPI AbsoluteQuantitationMethod
  {
public:
    bool operator==(const AbsoluteQuantitationMethod& other) const;
    bool operator!=(const AbsoluteQuantitationMethod& other) const;

    void setComponentName(const String& component_name); ///< Component name setter
    String getComponentName() const; ///< Component name getter

    void setFeatureName(const String& feature_name); ///< Feature name setter
    String getFeatureName() const; ///< Feature name getter

    void setISName(const String& IS_name); ///< IS name setter
    String getISName() const; ///< IS_name getter

    void setLLOD(const double llod); ///< LLOD setter
    double getLLOD() const; ///< LLOD getter
    void setULOD(const double ulod); ///< ULOD setter
    double getULOD() const; ///< ULOD getter
    bool checkLOD(const double value) const; ///< This function checks if the value is within the limits of detection (LOD)

    void setLLOQ(const double lloq); ///< LLOQ setter
    double getLLOQ() const; ///< LLOQ getter
    void setULOQ(const double uloq); ///< ULOQ setter
    double getULOQ() const; ///< ULOQ getter
    bool checkLOQ(const double value) const; ///< This function checks if the value is within the limits of quantitation (LOQ)

    void setNPoints(const Int n_points); ///< Set the number of points
    Int getNPoints() const; ///< Get the number of points

    void setCorrelationCoefficient(const double correlation_coefficient); ///< Set the correlation coefficient
    double getCorrelationCoefficient() const; ///< Get the correlation coefficient

    void setConcentrationUnits(const String& concentration_units); ///< Concentration units setter
    String getConcentrationUnits() const; ///< Concentration units getter

    void setTransformationModel(const String& transformation_model); ///< Transformation model setter
    String getTransformationModel() const; ///< Transformation model getter

    void setTransformationModelParams(const Param& transformation_model_params); ///< Transformation model parameters setter
    Param getTransformationModelParams() const; ///< Transformation model parameters getter

private:
    Param transformation_model_params_; ///< transformation model parameters
    String component_name_; ///< id of the component
    String feature_name_; ///< name of the feature (i.e., peak_apex_int or peak_area)
    String IS_name_; ///< the internal standard (IS) name for the transition
    String concentration_units_; ///< concentration units of the component's concentration
    String transformation_model_; ///< transformation model
    double llod_ { 0.0 }; ///< lower limit of detection (LLOD) of the transition
    double ulod_ { 0.0 }; ///< upper limit of detection (ULOD) of the transition
    double lloq_ { 0.0 }; ///< lower limit of quantitation (LLOQ) of the transition
    double uloq_ { 0.0 }; ///< upper limit of quantitation (ULOQ) of the transition
    double correlation_coefficient_ { 0.0 }; ///< the Pearson R value for the correlation coefficient of the calibration curve
    Int n_points_ { 0 }; ///< number of points used in a calibration curve
  };
}

