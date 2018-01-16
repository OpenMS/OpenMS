// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATIONMETHOD_H
#define OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATIONMETHOD_H

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
    AbsoluteQuantitationMethod() = default; ///< Constructor
    ~AbsoluteQuantitationMethod() = default; ///< Destructor

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
    String component_name_; ///< id of the component
    String feature_name_; ///< name of the feature (i.e., peak_apex_int or peak_area)
    String IS_name_; ///< the internal standard (IS) name for the transition
    double llod_; ///< lower limit of detection (LLOD) of the transition
    double ulod_; ///< upper limit of detection (ULOD) of the transition
    double lloq_; ///< lower limit of quantitation (LLOQ) of the transition
    double uloq_; ///< upper limit of quantitation (ULOQ) of the transition
    Int n_points_; ///< number of points used in a calibration curve
    double correlation_coefficient_; ///< the Pearson R value for the correlation coefficient of the calibration curve
    String concentration_units_; ///< concentration units of the component's concentration
    String transformation_model_; ///< transformation model
    Param transformation_model_params_; ///< transformation model parameters
  };
}

#endif // OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATIONMETHOD_H
