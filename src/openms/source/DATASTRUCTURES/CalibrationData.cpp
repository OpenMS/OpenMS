// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/CalibrationData.h>

#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>


namespace OpenMS
{

  CalibrationData::CalibrationData() :
    data_(),
    use_ppm_(true),
    groups_()
  {

  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getMZ( Size i ) const
  {
    return data_[i].getMZ();
  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getRT( Size i ) const
  {
    return data_[i].getRT();
  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getIntensity( Size i ) const
  {
    return data_[i].getIntensity();
  }

  OpenMS::CalibrationData::const_iterator CalibrationData::begin() const
  {
    return data_.begin();
  }

  OpenMS::CalibrationData::const_iterator CalibrationData::end() const
  {
    return data_.end();
  }

  Size CalibrationData::size() const
  {
    return data_.size();
  }

  bool CalibrationData::empty() const
  {
    return data_.empty();
  }

  void CalibrationData::clear()
  {
    data_.clear();
  }

  void CalibrationData::setUsePPM( bool usePPM )
  {
    use_ppm_ = usePPM;
  }

  bool CalibrationData::usePPM() const
  {
    return use_ppm_;
  }

  void CalibrationData::insertCalibrationPoint( CalDataType::CoordinateType rt,
                                                CalDataType::CoordinateType mz_obs,
                                                CalDataType::IntensityType intensity,
                                                CalDataType::CoordinateType mz_ref,
                                                double weight,
                                                int group /*= -1*/ )
  {
    RichPeak2D p(Peak2D::PositionType(rt, mz_obs), intensity);
    p.setMetaValue("mz_ref", mz_ref);
    p.setMetaValue("ppm_error", Math::getPPM(mz_obs, mz_ref));
    p.setMetaValue("weight", weight);

    if (group >= 0)
    {
      p.setMetaValue("peakgroup", group);
      groups_.insert(group);
    }
    data_.push_back(p);
  }

  Size CalibrationData::getNrOfGroups() const
  {
    return groups_.size();
  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getError( Size i ) const
  {
    if (use_ppm_)
    {
      return data_[i].getMetaValue("ppm_error");
    }
    else
    {
      return (data_[i].getMZ() - getRefMZ(i));
    }
  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getRefMZ( Size i ) const
  {
    if (!data_[i].metaValueExists("mz_ref"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "getRefMz() received invalid point without meta data!");
    }

    return data_[i].getMetaValue("mz_ref");
  }

  CalibrationData::CalDataType::CoordinateType CalibrationData::getWeight( Size i ) const
  {
    if (!data_[i].metaValueExists("weight"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "getWeight() received invalid point without meta data!");
    }
    return data_[i].getMetaValue("weight");
  }

  int CalibrationData::getGroup( Size i ) const
  {
    if (!data_[i].metaValueExists("peakgroup"))
    { 
      return -1;
    }
    return data_[i].getMetaValue("peakgroup");
  }

  OpenMS::StringList CalibrationData::getMetaValues()
  {
    return ListUtils::create<String>("mz_ref,ppm_error,weight");
  }

  OpenMS::CalibrationData CalibrationData::median( double rt_left, double rt_right ) const
  {
    CalibrationData cd;
    cd.setUsePPM(this->usePPM());

    Size i = std::distance(data_.begin(), lower_bound(data_.begin(), data_.end(), rt_left, RichPeak2D::PositionLess()));
    Size ie = std::distance(data_.begin(), upper_bound(data_.begin(), data_.end(), rt_right, RichPeak2D::PositionLess()));
    if (i == ie)
    {
      return cd;
    }

    double rt = (rt_left + rt_right) / 2;

    for (std::set<int>::const_iterator it_group = groups_.begin(); it_group!= groups_.end(); ++it_group)
    { 
      std::vector<double> mzs, ints;
      double mz_ref(0);
      for (Size j = i; j < ie; ++j)
      {
        if (getGroup(j) == *it_group)
        {
          mzs.push_back(data_[j].getMZ());
          ints.push_back(data_[j].getIntensity());
          mz_ref = getRefMZ(j);
        }
      }
      if (ints.empty())
      { 
        continue; // no data points for this peak group in this RT range
      }
      double int_median = Math::median(ints.begin(), ints.end());
      cd.insertCalibrationPoint(rt, Math::median(mzs.begin(), mzs.end()), int_median, mz_ref, log(int_median));
    }
    return cd;
  }

  void CalibrationData::sortByRT()
  {
    std::sort(data_.begin(), data_.end(), RichPeak2D::PositionLess());
  }

} // namespace OpenMS
