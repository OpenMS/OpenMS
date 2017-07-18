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
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <math.h>    // std::log

namespace OpenMS
{

  TransformationModel::TransformationModel(const TransformationModel::DataPoints&, const Param&) :
    params_()
  {
  }

  TransformationModel::~TransformationModel()
  {
  }

  double TransformationModel::evaluate(double value) const
  {
    return value;
  }

  const Param& TransformationModel::getParameters() const
  {
    return params_;
  }

  void TransformationModel::getDefaultParameters(Param& params)
  {
    params.clear();
  }
  
  void TransformationModel::unWeightData(TransformationModel::DataPoints& data, const Param& params)
  {
    // unweight x values 
    std::vector<std::string> valid_weights;
    valid_weights = getValidXWeights();
    const std::string x_weight = params.exists("x_weight") ? params.getValue("x_weight") : "";
    if (!x_weight.empty() && checkValidWeight(x_weight, valid_weights) && !data.empty())
    {
      x_weight_ = x_weight;
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].first = unWeightDatum(data[i].first,x_weight_);
      }
    }
    else
    {
      x_weight_ = x_weight;
    }
    // unweight y values
    valid_weights = getValidYWeights();
    const std::string y_weight = params.exists("y_weight") ? params.getValue("y_weight") : "";
    if (!y_weight.empty() && checkValidWeight(y_weight, valid_weights) && !data.empty())
    {
      y_weight_ = y_weight;
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].second = unWeightDatum(data[i].second,y_weight_);
      }
    }   
    else
    {
      y_weight_ = y_weight;
    }
  }

  bool TransformationModel::checkValidWeight(const std::string& weight, const std::vector<std::string>& valid_weights) const
  {    
    bool valid = false;
    if (std::find(valid_weights.begin(), valid_weights.end(), weight) != valid_weights.end())
    {
      valid=true;
    }
    else
    {
      LOG_INFO << "weight " + weight + " is not supported.";
    }
    return valid;
  }

  double TransformationModel::checkDatumRange(const double& datum, const double& datum_min, const double& datum_max)
  {    
    double datum_checked = datum;
    if (datum >= datum_max)
    {
      LOG_INFO << "datum " << datum << " is not out of range.";
      LOG_INFO << "datum will be truncated to " << datum_max << ".";
      datum_checked = datum_max;
    }
    else if (datum <= datum_min)
    {
      LOG_INFO << "datum " << datum << " is not out of range.";
      LOG_INFO << "datum will be truncated to " << datum_min << ".";
      datum_checked = datum_min;
    }
    return datum_checked;
  }
  
  std::vector<std::string> TransformationModel::getValidXWeights() const
  {
    //std::vector<std::string> valid_weights{"1/x","1/x2","ln(x)",""}; C++ 11
    std::vector<std::string> valid_weights;
    valid_weights.push_back("1/x");
    valid_weights.push_back("1/x2");
    valid_weights.push_back("ln(x)");
    valid_weights.push_back("");
    return valid_weights;
  }
  
  std::vector<std::string> TransformationModel::getValidYWeights() const
  {
    // std::vector<std::string> valid_weights{"1/y","1/y2","ln(y)",""}; C++ 11
    std::vector<std::string> valid_weights;
    valid_weights.push_back("1/y");
    valid_weights.push_back("1/y2");
    valid_weights.push_back("ln(y)");
    valid_weights.push_back("");
    return valid_weights;
  }

  double TransformationModel::weightDatum(const double& datum, const std::string& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      if (datum < 10e-5)
      {
        datum_weighted = std::log(10e-5);
      }
      else
      {
        datum_weighted = std::log(datum);
      }
    }
    else if (weight == "ln(y)")
    {
      if (datum < 10e-8)
      {
        datum_weighted = std::log(10e-8);
      }
      else
      {
        datum_weighted = std::log(datum);
      }
    }
    else if (weight == "1/x")
    {
      if (datum < 10e-5)
      {
        datum_weighted = 1/10e-5;
      }
      else
      {
        datum_weighted = 1/std::abs(datum);
      }
    }
    else if (weight == "1/y")
    {
      if (datum < 10e-8)
      {
        datum_weighted = 1/10e-8;
      }
      else
      {
        datum_weighted = 1/std::abs(datum);
      }
    }
    else if (weight == "1/x2")
    {
      if (datum < 10e-5)
      {
        datum_weighted = 1/std::pow(10e-5,2);
      }
      else
      {
        datum_weighted = 1/std::pow(datum,2);
      }
    }
    else if (weight == "1/y2")
    {
      if (datum < 10e-8)
      {
        datum_weighted = 1/std::pow(10e-8,2);
      }
      else
      {
        datum_weighted = 1/std::pow(datum,2);
      }
    }
    else if (weight == "")
    {
      datum_weighted = datum;
    }
    else
    {
      datum_weighted = datum;
      LOG_INFO << "weight " + weight + " not supported.";
      LOG_INFO << "no weighting will be applied.";
    }
    return datum_weighted;
  } 

  double TransformationModel::unWeightDatum(const double& datum, const std::string& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      if (datum > std::log(10e5))
      {
        datum_weighted = 10e5;
      }
      else
      {
        datum_weighted = std::abs(std::exp(datum));
      }
    }
    else if (weight == "ln(y)")
    {
      if (datum > std::log(10e8))
      {
        datum_weighted = 10e8;
      }
      else
      {
        datum_weighted = std::exp(datum);
      }
    }
    else if (weight == "1/x")
    {
      if (datum > 1/std::abs(10e-5))
      {
        datum_weighted = 10e-5;
      }
      else
      {
        datum_weighted = 1/std::abs(datum);
      }
    }
    else if (weight == "1/y")
    {
      if (datum > 1/std::abs(10e-8))
      {
        datum_weighted = 10e-8;
      }
      else
      {
        datum_weighted = 1/std::abs(datum);
      }
    }
    else if (weight == "1/x2")
    {
      if (datum > 1/std::pow(10e-5,2))
      {
        datum_weighted = 10e-5;
      }
      else
      {
        datum_weighted = std::sqrt(1/std::abs(datum));
      }
    }
    else if (weight == "1/y2")
    {
      if (datum >  1/std::pow(10e-8,2))
      {
        datum_weighted = 10e-8;
      }
      else
      {
        datum_weighted = std::sqrt(1/std::abs(datum));
      }
    }
    else if (weight == "")
    {
      datum_weighted = datum;
    }
    else
    {
      datum_weighted = datum;
      LOG_INFO << "weight " + weight + " not supported.";
      LOG_INFO << "no weighting will be applied.";
    }
    return datum_weighted;
  }

}
