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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//


#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <functional>
#include <numeric>
#include <fstream>

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <OpenMS/KERNEL/Peak1D.h>


using namespace std;

namespace OpenMS
{

  IsotopeDistribution::IsotopeDistribution()
  {}

  IsotopeDistribution::IsotopeDistribution(const IsotopeDistribution & isotope_distribution) :
    distribution_(isotope_distribution.distribution_)
  {}

  IsotopeDistribution::~IsotopeDistribution()
  {}

  

  IsotopeDistribution & IsotopeDistribution::operator=(const IsotopeDistribution & iso)
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

  const IsotopeDistribution::ContainerType & IsotopeDistribution::getContainer() const
  {
    return distribution_;
  }

  Peak1D::CoordinateType IsotopeDistribution::getMax() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[distribution_.size() - 1].getMZ();
  }

  Peak1D::CoordinateType IsotopeDistribution::getMin() const
  {
    if (distribution_.empty())
    {
      return 0;
    }
    return distribution_[0].getMZ();
  }


  Size IsotopeDistribution::size() const
  {
    return distribution_.size();
  }

  void IsotopeDistribution::clear()
  {
    distribution_.clear();
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

  void IsotopeDistribution::sort_(function<bool(const MassAbundance& p1, const MassAbundance& p2)> sorter)
  {
    sort(distribution_.begin(), distribution_.end(), sorter);
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
    for_each(distribution_.begin(), distribution_.end(), lambda);
  }


  bool IsotopeDistribution::operator==(const IsotopeDistribution & isotope_distribution) const
  {
    return distribution_ == isotope_distribution.distribution_;
  }

  bool IsotopeDistribution::operator!=(const IsotopeDistribution & isotope_distribution) const
  {
    return !(isotope_distribution == *this);
  }

  void IsotopeDistribution::renormalize()
  {
    if (distribution_.size() != 0)
    {
      double sum(0);
      // loop backwards as most distributions contains a lot of small values at the end
      for (ContainerType::const_reverse_iterator it = distribution_.rbegin(); it != distribution_.rend(); ++it)
      {
        sum += it->getIntensity();
      }

      for (Iterator it = distribution_.begin(); it != distribution_.end(); ++it)
      {
        it->setIntensity(it->getIntensity() / sum);
      }
    }
    return;
  }

  void IsotopeDistribution::trimRight(double cutoff)
  {
    ContainerType::reverse_iterator riter = distribution_.rbegin();

    // loop from right to left until an entry is larger than the cutoff
    for (; riter != distribution_.rend(); ++riter)
    {
      if (riter->getIntensity() >= cutoff)
        break;
    }
    // trim the container
    distribution_.resize(riter.base() - distribution_.begin());
  }

  void IsotopeDistribution::trimLeft(double cutoff)
  {
    for (ContainerType::iterator iter = distribution_.begin(); iter != distribution_.end(); ++iter)
    {
      if (iter->getIntensity() >= cutoff)
      {
        distribution_.erase(distribution_.begin(), iter);
        break;
      }
    }
  }

  bool IsotopeDistribution::isNormalized() const
  {
    return distribution_.front().getIntensity() == 1.0 && 
      is_sorted(distribution_.begin(), distribution_.end(),[](const MassAbundance& fr1, const MassAbundance& fr2){
          return fr1.getIntensity() > fr2.getIntensity();
        });
  }

  bool IsotopeDistribution::isConvolutionUnit() const
  { 
    return distribution_.size() == 1  && distribution_.front().getMZ() == 0.0;
  }

  double IsotopeDistribution::averageMass() const
  {
    double prob_sum = accumulate(distribution_.begin(), distribution_.end(), 0.0,
               [](double total_prob, const Peak1D& iso)
               {
                 return  total_prob + iso.getIntensity();
               });
    return accumulate(distribution_.begin(), distribution_.end(), 0.0,
                      [&prob_sum](double average_mass, const Peak1D& iso)
                      {
                        return average_mass + iso.getMZ()*(iso.getIntensity()/prob_sum);
                      });
  }


}
