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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

//Kernal classes
#include <OpenMS/KERNEL/StandardTypes.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{
    
  AbsoluteQuantitationMethod::AbsoluteQuantitationMethod()
  {
  }
  
  AbsoluteQuantitationMethod::~AbsoluteQuantitationMethod()
  {
  }
  
  //LOD getters and setters
  void AbsoluteQuantitationMethod::setLLOD(const double& llod)
  {
    llod_ = llod;
  }
  
  void AbsoluteQuantitationMethod::setULOD(const double& ulod)
  {
    ulod_ = ulod;
  }
    
  double AbsoluteQuantitationMethod::getLLOD()
  {
    return ulod_;
  }
  
  double AbsoluteQuantitationMethod::getULOD()
  {
    return ulod_;
  }
  
  void AbsoluteQuantitationMethod::setLLOQ(const double& lloq)
  {
    lloq_ = lloq;
  }
  
  void AbsoluteQuantitationMethod::setULOQ(const double& uloq)
  {
    uloq_ = uloq;
  }
  
  double AbsoluteQuantitationMethod::getLLOQ()
  {
    return lloq_;
  }
  
  double AbsoluteQuantitationMethod::getULOQ()
  {
    return uloq_;
  }
  
  //Component, IS, and Feature name setters  
  void AbsoluteQuantitationMethod::setFeatureName(const String& feature_name)
  {
    feature_name_ = feature_name;
  }
  
  String AbsoluteQuantitationMethod::getFeatureName()
  {
    return feature_name_;
  } 
  
  void AbsoluteQuantitationMethod::setISName(const String& IS_name)
  {
    IS_name_ = IS_name;
  }
  
  String AbsoluteQuantitationMethod::getISName()
  {
    return IS_name_;
  }
  
  void AbsoluteQuantitationMethod::setComponentName(const String& component_name)
  {
    component_name_ = component_name;
  }
  
  String AbsoluteQuantitationMethod::getComponentName()
  {
    return component_name_;
  }
  
  //Concentration unit getter and setter
  void AbsoluteQuantitationMethod::setConcentrationUnits(const String& concentration_units)
  {
    concentration_units_ = concentration_units;
  }

  String AbsoluteQuantitationMethod::getConcentrationUnits()
  {
    return concentration_units_;
  }
  
  //Transformation model getters and setters
  void AbsoluteQuantitationMethod::setTransformationModel(const String& transformation_model)
  {
    transformation_model_ = transformation_model;
  }
  void AbsoluteQuantitationMethod::setTransformationModelParams(const Param& transformation_model_params)
  {
    transformation_model_params_ = transformation_model_params;
  }

  String AbsoluteQuantitationMethod::getTransformationModel()
  {
    return transformation_model_;
  }

  Param AbsoluteQuantitationMethod::getTransformationModelParams()
  {
    return transformation_model_params_;
  }
  
  //Actual concentration getter and setter
  void AbsoluteQuantitationMethod::setActualConcentration(const double& actual_concentration)
  {
    actual_concentration_ = actual_concentration;
  }

  double AbsoluteQuantitationMethod::getActualConcentration()
  {
    return actual_concentration_;
  }
  
  //Statistics getters and setters
  void AbsoluteQuantitationMethod::setNPoints(const int& n_points)
  {
    n_points_ = n_points;
  }
  void AbsoluteQuantitationMethod::setCorrelationCoefficient(const double& correlation_coefficient)
  {
    correlation_coefficient_ = correlation_coefficient;
  }
  
  int AbsoluteQuantitationMethod::getNPoints()
  {
    return n_points_;
  }
  double AbsoluteQuantitationMethod::getCorrelationCoefficient()
  {
    return correlation_coefficient_;
  }

  //Non getter/setter methods
  bool AbsoluteQuantitationMethod::checkLOD(const double & value)
  {
    bool bracketted = false;
    if (value <= ulod_ && value >= llod_)
    {
      bracketted = true;
    }
    return bracketted;
  }

  bool AbsoluteQuantitationMethod::checkLOQ(const double & value)
  {
    bool bracketted = false;
    if (value <= uloq_ && value >= lloq_)
    {
      bracketted = true;
    }
    return bracketted;
  }

  Param AbsoluteQuantitationMethod::fitTransformationModel(const String & transformation_model,
    const TransformationModel::DataPoints& data,
    const Param& transformation_model_params)
  {
    Param params;
    if (transformation_model == "TransformationModelLinear")
    {
      TransformationModelLinear tm(data, transformation_model_params);
      params = tm.getParameters();
    }
    else if (transformation_model == "TransformationModelBSpline")
    {
      TransformationModelBSpline tm(data, transformation_model_params);
      params = tm.getParameters();
    }
    else if (transformation_model == "TransformationModelInterpolated")
    {
      TransformationModelInterpolated tm(data, transformation_model_params);
      params = tm.getParameters();
    }
    else if (transformation_model == "TransformationModelLowess")
    {
      TransformationModelLowess tm(data, transformation_model_params);
      params = tm.getParameters();
    }
    else
    {
      LOG_INFO << "TransformationModel " << transformation_model << " is not supported.";
      LOG_INFO << "default TransformationModel will be used.";
      TransformationModel tm(data, transformation_model_params);
      params = tm.getParameters();
    }
    return params;
  }
  
  double AbsoluteQuantitationMethod::evaluateTransformationModel(const String & transformation_model,
    const double& datum,
    const Param& transformation_model_params)
  {
    double result = datum;
    TransformationModel::DataPoints data;
    if (transformation_model == "TransformationModelLinear")
    {
      TransformationModelLinear tm(data, transformation_model_params);
      tm.invert();
      result = tm.evaluate(datum);
    }
    else if (transformation_model == "TransformationModelBSpline")
    {
      TransformationModelBSpline tm(data, transformation_model_params);
      // tm.invert(); // not supported
      result = tm.evaluate(datum);
    }
    else if (transformation_model == "TransformationModelInterpolated")
    {
      TransformationModelInterpolated tm(data, transformation_model_params);
      // tm.invert(); // not supported
      result = tm.evaluate(datum);
    }
    else if (transformation_model == "TransformationModelLowess")
    {
      TransformationModelLowess tm(data, transformation_model_params);
      // tm.invert(); // not supported
      result = tm.evaluate(datum);
    }
    else
    {
      LOG_INFO << "TransformationModel " << transformation_model << " is not supported.";
      LOG_INFO << "The original datum will be returned.";
    }
    return result;
  }

} // namespace

