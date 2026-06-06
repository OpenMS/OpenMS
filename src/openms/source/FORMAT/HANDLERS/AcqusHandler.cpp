// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Guillaume Belz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS::Internal
{

  AcqusHandler::AcqusHandler(const std::string & filename)
  {
    params_.clear();

    std::ifstream is(filename.c_str());
    if (!is)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }

    //temporary variables
    std::string line;
    std::vector<std::string> strings(2);

    //read lines
    while (getline(is, line, '\n'))
    {
      if (line.size() < 5)
      {
        continue;                    // minimal string = "##x=x"
      }
      if (StringUtils::prefix(line, 2) !=StringUtils::toStr("##"))
      {
        continue;
      }
      if (StringUtils::split(line, '=', strings))
      {
        if (strings.size() == 2)
        {
          params_[strings[0].substr(2)] = StringUtils::trim(strings[1]);
        }
      }
    }

    // TOF calibration params
    dw_ = StringUtils::toDouble(params_[std::string("$DW")]);
    delay_ = (Size)StringUtils::toInt32(params_[std::string("$DELAY")]);
    ml1_ = StringUtils::toDouble(params_[std::string("$ML1")]);
    ml2_ = StringUtils::toDouble(params_[std::string("$ML2")]);
    ml3_ = StringUtils::toDouble(params_[std::string("$ML3")]);
    td_ = (Size) StringUtils::toInt32(params_[std::string("$TD")]);

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

  std::string AcqusHandler::getParam(const std::string & param)
  {
    if (param ==StringUtils::toStr("mzMax"))
    {
      return StringUtils::toStr(getPosition(td_ - 1));
    }
    else if (param ==StringUtils::toStr("mzMin"))
    {
      return StringUtils::toStr(getPosition(0));
    }
    return params_[param];
  }

} // namespace OpenMS    // namespace Internal
