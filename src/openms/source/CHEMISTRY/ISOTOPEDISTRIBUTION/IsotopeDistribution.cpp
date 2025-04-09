// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>

using namespace std;

namespace OpenMS
{
  
  IsotopeDistribution::IsotopeDistribution()
  {
    distribution_.push_back(Peak1D(0, 1));
  }

  IsotopeDistribution& IsotopeDistribution::operator=(const IsotopeDistribution & iso)
  {
    if (this != &iso)
    {
      distribution_ = iso.distribution_;
    }
    return *this;
  }

  void IsotopeDistribution::set(const ContainerType & distribution)
  {
    distribution_ = distribution;
  }

  void IsotopeDistribution::set(ContainerType && distribution)
  {
    distribution_ = std::move(distribution);
  }

  const IsotopeDistribution::ContainerType& IsotopeDistribution::getContainer() const
  {
    return distribution_;
  }

  Peak1D::CoordinateType IsotopeDistribution::getMax() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return std::max_element(begin(), end(), MassAbundance::MZLess())->getMZ();
  }

  Peak1D::CoordinateType IsotopeDistribution::getMin() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return std::min_element(begin(), end(), MassAbundance::MZLess())->getMZ();
  }

  Peak1D IsotopeDistribution::getMostAbundant() const
  {
    if (distribution_.empty())
    {
        return Peak1D(0, 1);
    }
    return *std::max_element(begin(), end(), MassAbundance::IntensityLess());
  }

  Size IsotopeDistribution::size() const
  {
    return distribution_.size();
  }

  void IsotopeDistribution::clear()
  {
    distribution_.clear();
  }

  void IsotopeDistribution::resize(UInt new_size)
  {
    distribution_.resize(new_size);
  }

  void IsotopeDistribution::trimIntensities(double cutoff)
  {
    distribution_.erase(
      remove_if(distribution_.begin(),
                distribution_.end(),
                [&cutoff](const MassAbundance& sample)
                {
                  return sample.getIntensity() < cutoff;
                }), distribution_.end());
  }

  void IsotopeDistribution::sort_(
    function<bool(const MassAbundance& p1, const MassAbundance& p2)> sorter)
  {
    sort(distribution_.begin(), distribution_.end(), std::move(sorter));
  }

  void IsotopeDistribution::sortByIntensity()
  {
    sort_([](const MassAbundance& p1, const MassAbundance& p2){
        return p1.getIntensity() > p2.getIntensity();
      });
  }

  void IsotopeDistribution::sortByMass()
  {
    sort_([](const MassAbundance& p1, const MassAbundance& p2){
        return p1.getMZ() < p2.getMZ();
      });
  }

  void IsotopeDistribution::transform_(function<void(MassAbundance&)> lambda)
  {
    for_each(distribution_.begin(), distribution_.end(), std::move(lambda));
  }

  bool IsotopeDistribution::operator==(const IsotopeDistribution & isotope_distribution) const
  {
    return distribution_ == isotope_distribution.distribution_;
  }

  bool IsotopeDistribution::operator<(const IsotopeDistribution & rhs) const
  {
    if (distribution_.size() != rhs.distribution_.size()) 
    { 
      return distribution_.size() < rhs.distribution_.size(); 
    }

    // both vectors have same size
    auto it = distribution_.begin();
    auto rhs_it = rhs.distribution_.begin();
    for (; it != distribution_.end(); ++it, ++rhs_it)
    {
      if (*it != *rhs_it) 
      { 
        const double mz = it->getMZ();
        const double in = it->getIntensity();
        const double rhs_mz = rhs_it->getMZ();
        const double rhs_in = rhs_it->getIntensity();
        
        return tie(mz, in) < tie(rhs_mz, rhs_in);
      }
    }

    return false;
  }

  bool IsotopeDistribution::operator!=(const IsotopeDistribution & isotope_distribution) const
  {
    return !(isotope_distribution == *this);
  }

  void IsotopeDistribution::renormalize()
  {
    if (!distribution_.empty())
    {
      double sum(0);
      // loop backwards as most distributions contains a lot of small values at the end
      for (auto it = distribution_.rbegin(); it != distribution_.rend(); ++it)
      {
        sum += it->getIntensity();
      }

      for (Iterator it = distribution_.begin(); it != distribution_.end(); ++it)
      {
        it->setIntensity(it->getIntensity() / sum);
      }
    }
  }

  void IsotopeDistribution::trimRight(double cutoff)
  {
    auto riter = distribution_.rbegin();

    // loop from right to left until an entry is larger than the cutoff
    for (; riter != distribution_.rend(); ++riter)
    {
      if (riter->getIntensity() >= cutoff)
      {
        break;
      }
    }
    // trim the container
    distribution_.resize(riter.base() - distribution_.begin());
  }

  void IsotopeDistribution::trimLeft(double cutoff)
  {
    for (auto iter = distribution_.begin(); iter != distribution_.end(); ++iter)
    {
      if (iter->getIntensity() >= cutoff)
      {
        distribution_.erase(distribution_.begin(), iter);
        break;
      }
    }
  }

  double IsotopeDistribution::averageMass() const
  {
    double prob_sum = accumulate(distribution_.begin(),
                                 distribution_.end(), 
                                 0.0,
                                 [](double total_prob, const Peak1D& iso)
                                 {
                                   return  total_prob + iso.getIntensity();
                                 });

    return accumulate(distribution_.begin(), distribution_.end(), 0.0,
                      [&prob_sum](double average_mass, const Peak1D& iso)
                      {
                        return average_mass + 
                          iso.getMZ() * (iso.getIntensity() / prob_sum);
                      });
  }

  void IsotopeDistribution::merge(double resolution, double min_prob)
  {
    // Sort by mass and trim the tails of the container
    sortByMass();
    trimLeft(min_prob);
    trimRight(min_prob);
    
    ContainerType raw = distribution_;
    double mass_range = (raw.back().getMZ() - raw.front().getMZ());
    UInt output_size = ceil(mass_range / resolution);
    if (output_size > distribution_.size())
    {
      throw Exception::IllegalArgument(__FILE__,
                                       __LINE__,
                                       OPENMS_PRETTY_FUNCTION,
                                       "New Isotope Distribution "
                                       "has more points than the old one.");
    }

    distribution_.clear();
    ContainerType distribution(output_size, Peak1D(0, 0));
    double delta = mass_range / output_size;

    for (auto& p : raw)
    {
      UInt index = round((p.getMZ() - raw.front().getMZ())/resolution);
      if (index >= distribution.size()) {continue;}
      double mass = raw.front().getMZ() + (index * delta);
      distribution[index].setMZ(mass);
      distribution[index].setIntensity(distribution[index].getIntensity() + p.getIntensity());
    }
    distribution_ = distribution;
    trimIntensities(min_prob);
  }


} // namespace OpenMS

