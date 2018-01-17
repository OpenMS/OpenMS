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

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

namespace OpenMS
{
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
