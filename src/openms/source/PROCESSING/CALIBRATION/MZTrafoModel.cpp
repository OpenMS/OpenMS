// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>

#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/RANSAC/RANSACModelQuadratic.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  MZTrafoModel::MZTrafoModel()
   : coeff_(),
     use_ppm_(true),
     rt_(std::numeric_limits<double>::quiet_NaN())
  {
  }

  MZTrafoModel::MZTrafoModel(bool ppm_model)
    : coeff_(),
      use_ppm_(ppm_model),
      rt_(std::numeric_limits<double>::quiet_NaN())
  {
  }

  const std::string MZTrafoModel::names_of_modeltype[] = {"linear", "linear_weighted", "quadratic", "quadratic_weighted", "size_of_modeltype"};

  Math::RANSACParam* MZTrafoModel::ransac_params_ = nullptr;
  int MZTrafoModel::ransac_seed_ = time(nullptr);
  double MZTrafoModel::limit_offset_ = std::numeric_limits<double>::max(); // no limit by default
  double MZTrafoModel::limit_scale_ = std::numeric_limits<double>::max(); // no limit by default
  double MZTrafoModel::limit_power_ = std::numeric_limits<double>::max(); // no limit by default

  MZTrafoModel::MODELTYPE MZTrafoModel::nameToEnum(const std::string& name)
  {
    const std::string* qb = names_of_modeltype;
    const std::string* qe = qb + (int)SIZE_OF_MODELTYPE;
    const std::string* qm = std::find(qb, qe, name);
    return (MODELTYPE)std::distance(qb, qm);
  }

  const std::string& MZTrafoModel::enumToName(MZTrafoModel::MODELTYPE mt)
  {
    return names_of_modeltype[mt];
  }

  void MZTrafoModel::setRANSACParams(const Math::RANSACParam& p)
  {
    delete ransac_params_;
    ransac_params_ = new Math::RANSACParam(p);
  }

  void MZTrafoModel::setRANSACSeed(int seed)
  {
    ransac_seed_ = seed;
  }

  void MZTrafoModel::setCoefficientLimits(double offset, double scale, double power)
  {
    limit_offset_ = fabs(offset);
    limit_scale_ = fabs(scale);
    limit_power_ = fabs(power);
  }

  bool MZTrafoModel::isValidModel( const MZTrafoModel& trafo )
  {
    if (trafo.coeff_.empty()) return false;

    // go through coefficients and see if they are too extreme
    if (limit_offset_ < fabs(trafo.coeff_[0]))
    {
      return false;
    }
    if (limit_scale_ < fabs(trafo.coeff_[1]))
    {
      return false;
    }
    if (limit_power_ < fabs(trafo.coeff_[2]))
    {
      return false;
    }
    return (true);
  }

  bool MZTrafoModel::isTrained() const
  {
    return !coeff_.empty();
  }

  double MZTrafoModel::getRT() const
  {
    return rt_;
  }

  double MZTrafoModel::predict( double mz ) const
  {
    // mz = a + b * mz + c * mz^2
    double predict =
      coeff_[0] +
      coeff_[1] * mz +
      coeff_[2] * mz * mz;
    if (use_ppm_) // the above prediction is the ppm error
    { // ... so we convert to actual mass diff
      predict = Math::ppmToMass(-predict, mz) + mz;
    }
    else
    {
      predict = (-predict) + mz;
    }
    return predict;
  }

  bool MZTrafoModel::train( const CalibrationData& cd, MODELTYPE md, bool use_RANSAC, double rt_left /*= -std::numeric_limits<double>::max()*/, double rt_right /*= std::numeric_limits<double>::max() */ )
  {
    std::vector<double> obs_mz;
    std::vector<double> theo_mz;
    std::vector<double> weights;
    const CalibrationData* p_cd;
    CalibrationData cdm;
    Size i, ie; // CalibrationData's start-to-end interval
    if (cd.getNrOfGroups() > 0) // we have lock mass traces
    { // this is extra work, since we need to collect peak groups and compute the median
      cdm = cd.median(rt_left, rt_right);
      p_cd = &cdm;
      i = 0;
      ie = cdm.size();
    }
    else
    {
      i = std::distance(cd.begin(), lower_bound(cd.begin(), cd.end(), rt_left, RichPeak2D::RTLess()));
      ie = std::distance(cd.begin(), upper_bound(cd.begin(), cd.end(), rt_right, RichPeak2D::RTLess()));
      p_cd = &cd;
    }
    for (Size j = i; j != ie; ++j)
    {
      obs_mz.push_back(p_cd->getError(j)); // could be ppm or [Th], depending on cd::use_ppm_
      theo_mz.push_back(p_cd->getRefMZ(j));
      weights.push_back(p_cd->getWeight(j));
    }

    this->rt_ = (rt_left + rt_right) / 2;

    return (train(obs_mz, theo_mz, weights, md, use_RANSAC));
  }

  bool MZTrafoModel::train( std::vector<double> obs_mz, std::vector<double> theo_mz, std::vector<double> weights, MODELTYPE md, bool use_RANSAC )
  {
    coeff_.clear();

    if (obs_mz.empty())
    {
      //OPENMS_LOG_ERROR << "Input to calibration model is empty!" << std::endl;
      return false;
    }

    if (use_RANSAC)
    {
      if (ransac_params_ == nullptr)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "TrafoModel::train(): no RANSAC parameters were set before calling train(). Internal error!");
      }
      if (!(md == LINEAR || md == QUADRATIC))
      {
        OPENMS_LOG_ERROR << "RANSAC is implemented for LINEAR and QUADRATIC models only! Please disable RANSAC or choose the LINEAR or QUADRATIC model." << std::endl;
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    }

    try
    {
      if (md == LINEAR)
      {
        if (obs_mz.size() < 2)
        {
          return false;
        }
        if (use_RANSAC && 
          (obs_mz.size() > ransac_params_->n)) // with fewer points, RANSAC will fail
        {
          std::vector<std::pair<double, double> > r, pairs;
          for (Size i = 0; i < obs_mz.size(); ++i)
          {
            pairs.emplace_back(theo_mz[i], obs_mz[i]);
          }
          r = Math::RANSAC<Math::RansacModelLinear>(ransac_seed_).ransac(pairs, *ransac_params_);
          if (r.size() < 2)
          {
            return false; // RANSAC failed
          }
          obs_mz.clear();
          theo_mz.clear();
          for (Size i = 0; i < r.size(); ++i)
          {
            theo_mz.push_back(r[i].first);
            obs_mz.push_back(r[i].second);
          }
        }

        double confidence_interval_P(0.0);
        Math::LinearRegression lr;
        lr.computeRegression(confidence_interval_P, theo_mz.begin(), theo_mz.end(), obs_mz.begin(), false);
        coeff_.push_back(lr.getIntercept());
        coeff_.push_back(lr.getSlope());
        coeff_.push_back(0.0);
      }
      else if (md == LINEAR_WEIGHTED)
      {
        if (obs_mz.size() < 2)
        {
          return false;
        }
        double confidence_interval_P(0.0);
        Math::LinearRegression lr;
        lr.computeRegressionWeighted(confidence_interval_P, theo_mz.begin(), theo_mz.end(), obs_mz.begin(), weights.begin(), false);
        coeff_.push_back(lr.getIntercept());
        coeff_.push_back(lr.getSlope());
        coeff_.push_back(0.0);
      }
      else if (md == QUADRATIC)
      {
        if (obs_mz.size() < 3)
        {
          return false;
        }
        if (use_RANSAC && 
          (obs_mz.size() > ransac_params_->n)) // with fewer points, RANSAC will fail
        {

          std::vector<std::pair<double, double> > r, pairs;
          for (Size i = 0; i < obs_mz.size(); ++i)
          {
            pairs.emplace_back(theo_mz[i], obs_mz[i]);
          }
          r = Math::RANSAC<Math::RansacModelQuadratic>(ransac_seed_).ransac(pairs, *ransac_params_);
          obs_mz.clear();
          theo_mz.clear();
          for (Size i = 0; i < r.size(); ++i)
          {
            theo_mz.push_back(r[i].first);
            obs_mz.push_back(r[i].second);
          }
        }
        // Quadratic fit
        Math::QuadraticRegression qr;
        qr.computeRegression(theo_mz.begin(), theo_mz.end(), obs_mz.begin());
        coeff_.push_back(qr.getA());
        coeff_.push_back(qr.getB());
        coeff_.push_back(qr.getC());
      }
      else if (md == QUADRATIC_WEIGHTED)
      {
        if (obs_mz.size() < 3)
        {
          return false;
        }
        // Quadratic fit (weighted)
        Math::QuadraticRegression qr;
        qr.computeRegressionWeighted(theo_mz.begin(), theo_mz.end(), obs_mz.begin(), weights.begin());
        coeff_.push_back(qr.getA());
        coeff_.push_back(qr.getB());
        coeff_.push_back(qr.getC());
      }

#ifdef DEBUG_TRAFOMODEL
      printf("# mz regression parameters: Y = %3.10f + %3.10f X + %3.10f X^2\n",
        coeff_[0],
        coeff_[1],
        coeff_[2]);

      // print results
      std::cout << "Calibration details:\n\n";
      std::cout << "m/z(theo) m/z(obs) ppm(before) | ppm(after)\n";
      std::vector<double> st_ppm_before, st_ppm_after;
      for (Size i = 0; i < obs_mz.size(); i++)
      {
        if (use_ppm_)
        {
          st_ppm_before.push_back(obs_mz[i]);
        }
        else 
        {
          st_ppm_before.push_back(Math::getPPM(theo_mz[i], obs_mz[i]));
        }
        double obs_mz_v = obs_mz[i];
        if (use_ppm_)
        {
          obs_mz_v = Math::ppmToMass(obs_mz_v, theo_mz[i]) + theo_mz[i];
        }
        st_ppm_after.push_back(Math::getPPM(theo_mz[i], predict(obs_mz_v))); // predict() is ppm-aware itself

        printf("%4.5f  %4.5f  %2.1f | %2.1f\n", theo_mz[i], obs_mz_v, st_ppm_before.back(), st_ppm_after.back());
      }
      // use median and MAD to ignore outliers
      double m = Math::median(st_ppm_before.begin(), st_ppm_before.end());
      std::cout << "ppm before: median = " << m << "  MAD = " << Math::MAD(st_ppm_before.begin(), st_ppm_before.end(), m) << "\n";
      m = Math::median(st_ppm_after.begin(), st_ppm_after.end()); 
      std::cout << "ppm after : median = " << m << "  MAD = " << Math::MAD(st_ppm_after.begin(), st_ppm_after.end(), m) << "\n";
#endif

      return true;
    }
    catch (Exception::BaseException& /*e*/)
    {
      //OPENMS_LOG_ERROR << "Exception during model fitting: " << e.what() << std::endl;
      return false;
    }
  }

  Size MZTrafoModel::findNearest( const std::vector<MZTrafoModel>& tms, double rt )
  {
    // no peak => no search
    if (tms.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one model to determine the nearest model!");
    }
    // search for position for inserting
    std::vector<MZTrafoModel>::const_iterator it = lower_bound(tms.begin(), tms.end(), rt, MZTrafoModel::RTLess());

    // border cases
    if (it == tms.begin())
    {
      return 0;
    }
    if (it == tms.end())
    {
      return tms.size() - 1;
    }
    // the model before or the current model are closest
    std::vector<MZTrafoModel>::const_iterator it2 = it;
    --it2;
    if (std::fabs(it->rt_ - rt) < std::fabs(it2->rt_ - rt))
    {
      return Size(it - tms.begin());
    }
    else
    {
      return Size(it2 - tms.begin());
    }
  }

  void MZTrafoModel::setCoefficients( const MZTrafoModel& rhs )
  {
    coeff_ = rhs.coeff_;
  }

  void MZTrafoModel::setCoefficients( double intercept, double slope, double power )
  {
    coeff_.clear();
    coeff_.push_back(intercept);
    coeff_.push_back(slope);
    coeff_.push_back(power);
  }

  OpenMS::String MZTrafoModel::toString() const
  {
    String s;
    if (coeff_.empty())
    {
      s = "nan, nan, nan";
    }
    else 
    {
      s = ListUtils::concatenate(coeff_, ", ");
    }
    return s;
  }

  void MZTrafoModel::getCoefficients( double& intercept, double& slope, double& power )
  {
    if (!isTrained())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Model is not trained yet.");
    }
    intercept = coeff_[0];
    slope = coeff_[1];
    power = coeff_[2];
  }


}
