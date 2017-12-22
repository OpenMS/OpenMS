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
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <math.h>    // std::log

namespace OpenMS
{

  TransformationModel::TransformationModel(const TransformationModel::DataPoints&, const Param& p) :
    params_(p),
    x_weight_(""),
    x_datum_min_(0),
    x_datum_max_(0),
    y_weight_(""),
    y_datum_min_(0),
    y_datum_max_(0),
    weighting_(false)
  {
    // get x datum ranges
    x_datum_min_ = params_.exists("x_datum_min") ? (double)params_.getValue("x_datum_min") : 1e-15;
    x_datum_max_ = params_.exists("x_datum_max") ? (double)params_.getValue("x_datum_max") : 1e15;

    // get y datum ranges
    y_datum_min_ = params_.exists("y_datum_min") ? (double)params_.getValue("y_datum_min") : 1e-15;
    y_datum_max_ = params_.exists("y_datum_max") ? (double)params_.getValue("y_datum_max") : 1e15;

    y_weight_ = params_.exists("y_weight") ? (String)params_.getValue("y_weight") : "";
    x_weight_ = params_.exists("x_weight") ? (String)params_.getValue("x_weight") : "";

    std::vector<String> valid_x_weights = getValidXWeights();
    std::vector<String> valid_y_weights = getValidYWeights();
    if (!x_weight_.empty() && !checkValidWeight(x_weight_, valid_x_weights))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + x_weight_ + "' is not a valid weight parameter for x values.");
    }
    if (!y_weight_.empty() && !checkValidWeight(y_weight_, valid_y_weights))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + y_weight_ + "' is not a valid weight parameter for y values.");
    }

    // easily remember whether we do weighting or not
    weighting_ = (x_weight_ != "" || y_weight_ != ""); 
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
  
  void TransformationModel::weightData(TransformationModel::DataPoints& data)
  {
    if (!weighting_ ) return;

    // weight x values 
    if (!x_weight_.empty() && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        // check x datum ranges
        data[i].first = checkDatumRange(data[i].first,x_datum_min_,x_datum_max_);
        // weight x datum
        data[i].first = weightDatum(data[i].first,x_weight_);
      }
    }

    // weight y values
    if (!y_weight_.empty() && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        // check y datum ranges
        data[i].second = checkDatumRange(data[i].second,y_datum_min_,y_datum_max_);
        // weight y datum
        data[i].second = weightDatum(data[i].second,y_weight_);
      }
    } 
  }
  
  void TransformationModel::unWeightData(TransformationModel::DataPoints& data)
  {
    if (!weighting_ ) return;

    // unweight x values 
    if (!x_weight_.empty() && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].first = unWeightDatum(data[i].first, x_weight_);
      }
    }
    // unweight y values
    if (!y_weight_.empty() && !data.empty())
    {
      for (size_t i = 0; i < data.size(); ++i)
      {
        data[i].second = unWeightDatum(data[i].second, y_weight_);
      }
    }  
  }

  bool TransformationModel::checkValidWeight(const String& weight, const std::vector<String>& valid_weights) const
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
      LOG_INFO << "datum " << datum << " is out of range.";
      LOG_INFO << "datum will be truncated to " << datum_max << ".";
      datum_checked = datum_max;
    }
    else if (datum <= datum_min)
    {
      LOG_INFO << "datum " << datum << " is out of range.";
      LOG_INFO << "datum will be truncated to " << datum_min << ".";
      datum_checked = datum_min;
    }
    return datum_checked;
  }
  
  std::vector<String> TransformationModel::getValidXWeights() const
  {
    //std::vector<String> valid_weights{"1/x","1/x2","ln(x)",""}; C++ 11
    std::vector<String> valid_weights;
    valid_weights.push_back("1/x");
    valid_weights.push_back("1/x2");
    valid_weights.push_back("ln(x)");
    valid_weights.push_back("");
    return valid_weights;
  }
  
  std::vector<String> TransformationModel::getValidYWeights() const
  {
    // std::vector<String> valid_weights{"1/y","1/y2","ln(y)",""}; C++ 11
    std::vector<String> valid_weights;
    valid_weights.push_back("1/y");
    valid_weights.push_back("1/y2");
    valid_weights.push_back("ln(y)");
    valid_weights.push_back("");
    return valid_weights;
  }

  double TransformationModel::weightDatum(const double& datum, const String& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      datum_weighted = std::log(datum);
    }
    else if (weight == "ln(y)")
    {
      datum_weighted = std::log(datum);
    }
    else if (weight == "1/x")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/y")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/x2")
    {
      datum_weighted = 1/std::pow(datum,2);
    }
    else if (weight == "1/y2")
    {
      datum_weighted = 1/std::pow(datum,2);
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

  double TransformationModel::unWeightDatum(const double& datum, const String& weight) const
  { 
    double datum_weighted = 0;   
    if (weight == "ln(x)")
    {
      datum_weighted = std::exp(datum);
    }
    else if (weight == "ln(y)")
    {
      datum_weighted = std::exp(datum);
    }
    else if (weight == "1/x")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/y")
    {
      datum_weighted = 1/std::abs(datum);
    }
    else if (weight == "1/x2")
    {
      datum_weighted = std::sqrt(1/std::abs(datum));
    }
    else if (weight == "1/y2")
    {
      datum_weighted = std::sqrt(1/std::abs(datum));
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
