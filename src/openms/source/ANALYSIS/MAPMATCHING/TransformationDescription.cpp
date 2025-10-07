// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

#include <iomanip>
#include <iostream>

using namespace std;

namespace OpenMS
{
  TransformationDescription::TransformationDescription() :
    data_(TransformationDescription::DataPoints()),
    model_type_("none"),
    model_(new TransformationModel())
  {
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription::DataPoints& data) :
    data_(data), model_type_("none"),
    model_(new TransformationModel())
  {
  }

  TransformationDescription::~TransformationDescription()
  {
    delete model_;
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription& rhs)
  {
    data_ = rhs.data_;
    model_type_ = "none";
    model_ = nullptr; // initialize this before the "delete" call in "fitModel"!
    Param params = rhs.getModelParameters();
    fitModel(rhs.model_type_, params);
  }

  TransformationDescription& TransformationDescription::operator=(
    const TransformationDescription& rhs)
  {
    if (this == &rhs)
      return *this;

    data_ = rhs.data_;
    model_type_ = "none";
    Param params = rhs.getModelParameters();
    fitModel(rhs.model_type_, params);

    return *this;
  }

  void TransformationDescription::fitModel(const String& model_type,
                                           const Param& params)
  {
    // if (previous) transformation is the identity, don't fit another model:
    if (model_type_ == "identity") return;

    delete model_;
    model_ = nullptr; // avoid segmentation fault in case of exception
    if ((model_type == "none") || (model_type == "identity"))
    {
      model_ = new TransformationModel();
    }
    else if (model_type == "linear")
    {
      model_ = new TransformationModelLinear(data_, params);
      // // debug output:
      // double slope, intercept;
      // TransformationModelLinear* lm = dynamic_cast<TransformationModelLinear*>(model_);
      // lm->getParameters(slope, intercept);
      // cout << "slope: " << slope << ", intercept: " << intercept << endl;
    }
    else if (model_type == "b_spline")
    {
      model_ = new TransformationModelBSpline(data_, params);
    }
    else if (model_type == "lowess")
    {
      model_ = new TransformationModelLowess(data_, params);
    }
    else if (model_type == "interpolated")
    {
      model_ = new TransformationModelInterpolated(data_, params);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unknown model type '" + model_type + "'");
    }
    model_type_ = model_type;
  }

  double TransformationDescription::apply(double value) const
  {
    return model_->evaluate(value);
  }

  const String& TransformationDescription::getModelType() const
  {
    return model_type_;
  }

  void TransformationDescription::getModelTypes(StringList& result)
  {
    result = ListUtils::create<String>("linear,b_spline,interpolated,lowess");
    // "none" and "identity" don't count
  }

  void TransformationDescription::setDataPoints(const DataPoints& data)
  {
    data_ = data;
    model_type_ = "none"; // reset the model even if it was "identity"
    delete model_;
    model_ = new TransformationModel();
  }

  void TransformationDescription::setDataPoints(const vector<pair<double, double> >& data)
  {
    data_.resize(data.size());
    for (Size i = 0; i < data.size(); ++i)
    {
      data_[i] = data[i];
    }
    model_type_ = "none"; // reset the model even if it was "identity"
    delete model_;
    model_ = new TransformationModel();
  }

  const TransformationDescription::DataPoints&
  TransformationDescription::getDataPoints() const
  {
    return data_;
  }

  const Param& TransformationDescription::getModelParameters() const
  {
    return model_->getParameters();
  }

  void TransformationDescription::invert()
  {
    for (TransformationDescription::DataPoints::iterator it = data_.begin();
         it != data_.end(); ++it)
    {
      *it = TransformationDescription::DataPoint(it->second, it->first,
                                                 it->note);
    }
    // ugly hack for linear model with explicit slope/intercept parameters:
    if ((model_type_ == "linear") && data_.empty())
    {
      TransformationModelLinear* lm =
        dynamic_cast<TransformationModelLinear*>(model_);
      lm->invert();
    }
    else
    {
      Param params = getModelParameters();
      fitModel(model_type_, params);
    }
  }

  void TransformationDescription::getDeviations(vector<double>& diffs, 
                                                bool do_apply,
                                                bool do_sort) const
  {
    diffs.clear();
    diffs.reserve(data_.size());
    for (DataPoints::const_iterator it = data_.begin(); it != data_.end(); ++it)
    {
      double x = it->first;
      if (do_apply) x = apply(x);
      diffs.push_back(abs(x - it->second));
    }
    if (do_sort) sort(diffs.begin(), diffs.end());
  }

  TransformationDescription::TransformationStatistics TransformationDescription::getStatistics() const
  {
    TransformationDescription::TransformationStatistics s;

    if (data_.empty()) return s;

    // x/y data ranges:
    double xmin, xmax, ymin, ymax;
    xmin = xmax = data_[0].first;
    ymin = ymax = data_[0].second;
    for (DataPoints::const_iterator it = ++data_.begin(); it != data_.end();
         ++it)
    {
      if (xmin > it->first) xmin = it->first;
      if (xmax < it->first) xmax = it->first;
      if (ymin > it->second) ymin = it->second;
      if (ymax < it->second) ymax = it->second;
    }

    s.xmin = xmin;
    s.xmax = xmax;
    s.ymin = ymin;
    s.ymax = ymax;

    // deviations:
    vector<double> diffs;
    getDeviations(diffs);
    bool no_model = (model_type_ == "none") || (model_type_ == "identity");

    for (const Size p : s.percents)
    {
      Size index = p / 100.0 * diffs.size() - 1;
      s.percentiles_before[p] = diffs[index];
    }

    // if we have a model, calculate deviations after applying the model
    // else set the same values
    if (!no_model) { getDeviations(diffs, true); }
  
    for (const Size p : s.percents)
    {
      Size index = p / 100.0 * diffs.size() - 1;
      s.percentiles_after[p] = diffs[index];
    }

    return s;
  }

  void TransformationDescription::printSummary(ostream& os) const
  {
    const TransformationStatistics s = getStatistics();

    os << "Number of data points (x/y pairs): " << data_.size() << "\n";
    if (data_.empty()) return;

    os << "Data range (x): " << s.xmin << " to " << s.xmax
       << "\nData range (y): " << s.ymin << " to " << s.ymax << "\n";

    // deviations:
    vector<double> diffs;
    getDeviations(diffs);
    bool no_model = (model_type_ == "none") || (model_type_ == "identity");
    os << String("Summary of x/y deviations") +
      (no_model ? "" : " before transformation") + ":\n";

    for (Size p : s.percents)
    {
      os << "- " << setw(3) << p << "% of data points within (+/-)"
         << s.percentiles_before.at(p) << "\n";
    }
    if (no_model)
    {
      os << endl;
      return;
    }
    // else:
    getDeviations(diffs, true);
    os << "Summary of x/y deviations after applying '" << model_type_
       << "' transformation:\n";

    for (Size p : s.percents)
    {
      os << "- " << setw(3) << p << "% of data points within (+/-)"
         << s.percentiles_after.at(p) << "\n";
    }
    os << endl;
  }

} // end of namespace OpenMS
