// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Guillaume Belz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS::Internal
{

  AcqusHandler::AcqusHandler(const String & filename)
  {
    params_.clear();

    std::ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    //temporary variables
    String line;
    std::vector<String> strings(2);

    //read lines
    while (getline(is, line, '\n'))
    {
      if (line.size() < 5)
      {
        continue;                    // minimal string = "##x=x"
      }
      if (line.prefix(2) != String("##"))
      {
        continue;
      }
      if (line.split('=', strings))
      {
        if (strings.size() == 2)
        {
          params_[strings[0].substr(2)] = strings[1].trim();
        }
      }
    }

    // TOF calibration params
    dw_ = params_[String("$DW")].toDouble();
    delay_ = (Size)params_[String("$DELAY")].toInt();
    ml1_ = params_[String("$ML1")].toDouble();
    ml2_ = params_[String("$ML2")].toDouble();
    ml3_ = params_[String("$ML3")].toDouble();
    td_ = (Size) params_[String("$TD")].toInt();

    is.close();
  }

  AcqusHandler::~AcqusHandler()
  {
    params_.clear();
  }

  Size AcqusHandler::getSize() const
  {
    return td_;
  }

  double AcqusHandler::getPosition(const Size index) const
  {
    double sqrt_mz_;
    double tof_ = dw_ * index + delay_;
    double a_ = ml3_;
    double b_ = sqrt(1000000000000.0 / ml1_);
    double c_ = ml2_ - tof_;

    if (ml3_ == 0.0)
    {
      sqrt_mz_ = c_ / b_;
    }
    else
    {
      sqrt_mz_ = (sqrt(b_ * b_ - 4 * a_ * c_) - b_) / (2 * a_);
    }
    return sqrt_mz_ * sqrt_mz_;
  }

  String AcqusHandler::getParam(const String & param)
  {
    if (param == String("mzMax"))
    {
      return String(getPosition(td_ - 1));
    }
    else if (param == String("mzMin"))
    {
      return String(getPosition(0));
    }
    return params_[param];
  }

} // namespace OpenMS    // namespace Internal
