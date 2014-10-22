// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  ConsensusFeature::ConsensusFeature() :
    BaseFeature(), handles_(), ratios_()
  {
  }

  ConsensusFeature::ConsensusFeature(const ConsensusFeature & rhs) :
    BaseFeature(rhs), handles_(rhs.handles_), ratios_()
  {
    ratios_ = rhs.ratios_;
  }

  ConsensusFeature::ConsensusFeature(const BaseFeature & feature) :
    BaseFeature(feature), handles_(), ratios_()
  {
  }

  ConsensusFeature::ConsensusFeature(UInt64 map_index, const Peak2D & element, UInt64 element_index) :
    BaseFeature(element), handles_(), ratios_()
  {
    insert(map_index, element, element_index);
  }

  ConsensusFeature::ConsensusFeature(UInt64 map_index, const BaseFeature & element) :
    BaseFeature(element), handles_(), ratios_()
  {
    insert(FeatureHandle(map_index, element));
  }

  ConsensusFeature & ConsensusFeature::operator=(const ConsensusFeature & rhs)
  {
    if (&rhs == this)
      return *this;

    BaseFeature::operator=(rhs);
    handles_ = rhs.handles_;
    ratios_ = rhs.ratios_;
    return *this;
  }

  ConsensusFeature::~ConsensusFeature()
  {
  }

  void ConsensusFeature::insert(const ConsensusFeature & cf)
  {
    handles_.insert(cf.handles_.begin(), cf.handles_.end());
  }

  void ConsensusFeature::insert(const FeatureHandle & handle)
  {
    if (!(handles_.insert(handle).second))
    {
      String key = String("map") + handle.getMapIndex() + "/feature" + handle.getUniqueId();
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The set already contained an element with this key.", key);
    }
  }

  void ConsensusFeature::insert(const HandleSetType & handle_set)
  {
    for (ConsensusFeature::HandleSetType::const_iterator it = handle_set.begin(); it != handle_set.end(); ++it)
    {
      insert(*it);
    }
  }

  void ConsensusFeature::insert(UInt64 map_index, const Peak2D & element, UInt64 element_index)
  {
    insert(FeatureHandle(map_index, element, element_index));
  }

  void ConsensusFeature::insert(UInt64 map_index, const BaseFeature & element)
  {
    insert(FeatureHandle(map_index, element));
    peptides_.insert(peptides_.end(), element.getPeptideIdentifications().begin(), element.getPeptideIdentifications().end());
  }

  const ConsensusFeature::HandleSetType & ConsensusFeature::getFeatures() const
  {
    return handles_;
  }

  DRange<2> ConsensusFeature::getPositionRange() const
  {
    DPosition<2> min = DPosition<2>::maxPositive();
    DPosition<2> max = DPosition<2>::minPositive();
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      if (it->getRT() < min[0])
        min[0] = it->getRT();
      if (it->getRT() > max[0])
        max[0] = it->getRT();
      if (it->getMZ() < min[1])
        min[1] = it->getMZ();
      if (it->getMZ() > max[1])
        max[1] = it->getMZ();
    }
    return DRange<2>(min, max);
  }

  DRange<1> ConsensusFeature::getIntensityRange() const
  {
    DPosition<1> min = DPosition<1>::maxPositive();
    DPosition<1> max = DPosition<1>::minPositive();
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      if (it->getIntensity() < min[0])
        min[0] = it->getIntensity();
      if (it->getIntensity() > max[0])
        max[0] = it->getIntensity();
    }
    return DRange<1>(min, max);
  }

  void ConsensusFeature::computeConsensus()
  {
    // for computing average position and intensity
    double rt = 0.0;
    double mz = 0.0;
    double intensity = 0.0;

    // The most frequent charge state wins.  Tie breaking prefers smaller charge.
    std::map<Int, UInt> charge_occ;
    Int charge_most_frequent = 0;
    UInt charge_most_frequent_occ = 0;

    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      rt += it->getRT();
      mz += it->getMZ();
      intensity += it->getIntensity();
      const Int it_charge = it->getCharge();
      const UInt it_charge_occ = ++charge_occ[it_charge];
      if (it_charge_occ > charge_most_frequent_occ)
      {
        charge_most_frequent_occ = it_charge_occ;
        charge_most_frequent = it_charge;
      }
      else
      {
        if (it_charge_occ >= charge_most_frequent_occ && abs(it_charge) < abs(charge_most_frequent))
        {
          charge_most_frequent = it_charge;
        }
      }
    }

    // compute the average position and intensity
    setRT(rt / size());
    setMZ(mz / size());
    setIntensity(intensity / size());
    setCharge(charge_most_frequent);
    return;
  }

  void ConsensusFeature::computeMonoisotopicConsensus()
  {
    // for computing average rt position, minimal m/z position and intensity
    double rt = 0.0;
    double mz = std::numeric_limits<double>::max();
    double intensity = 0.0;

    // The most frequent charge state wins.  Tie breaking prefers smaller charge.
    std::map<Int, UInt> charge_occ;
    Int charge_most_frequent = 0;
    UInt charge_most_frequent_occ = 0;

    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      rt += it->getRT();
      if (it->getMZ() < mz)
        mz = it->getMZ();
      intensity += it->getIntensity();
      const Int it_charge = it->getCharge();
      const UInt it_charge_occ = ++charge_occ[it_charge];
      if (it_charge_occ > charge_most_frequent_occ)
      {
        charge_most_frequent_occ = it_charge_occ;
        charge_most_frequent = it_charge;
      }
      else
      {
        if (it_charge_occ >= charge_most_frequent_occ && abs(it_charge) < abs(charge_most_frequent))
        {
          charge_most_frequent = it_charge;
        }
      }
    }

    // compute the position and intensity
    setRT(rt / size());
    setMZ(mz);
    setIntensity(intensity / size());
    setCharge(charge_most_frequent);
    return;
  }

  void ConsensusFeature::computeDechargeConsensus(const FeatureMap & fm, bool intensity_weighted_averaging)
  {
    // for computing average position and intensity
    double rt = 0.0;
    double m = 0.0;
    double intensity = 0.0;

    double proton_mass = Constants::PROTON_MASS_U;

    // intensity sum (for weighting)
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      intensity += it->getIntensity();
    }

    // unweighted averaging by default
    // TODO: add outlier removal
    // TODO: split cluster for each channel (in FD.C)
    double weighting_factor = 1.0 / size();

    // RT and Mass
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      Int q = it->getCharge();
      if (q == 0)
        LOG_WARN << "ConsensusFeature::computeDechargeConsensus() WARNING: Feature's charge is 0! This will lead to M=0!\n";
      double adduct_mass;
      Size index = fm.uniqueIdToIndex(it->getUniqueId());
      if (index > fm.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, fm.size());
      if (fm[index].metaValueExists("dc_charge_adduct_mass"))
      {
        adduct_mass = (double) fm[index].getMetaValue("dc_charge_adduct_mass");
      }
      else
      {
        adduct_mass = q * proton_mass;
      }

      if (intensity_weighted_averaging)
        weighting_factor = it->getIntensity() / intensity;
      rt += it->getRT() * weighting_factor;
      m += (it->getMZ() * q - adduct_mass) * weighting_factor;
    }

    // compute the average position and intensity
    setRT(rt);
    setMZ(m);
    setIntensity(intensity);
    setCharge(0);
    return;
  }

  void ConsensusFeature::addRatio(const ConsensusFeature::Ratio & r)
  {
    ratios_.push_back(r);
  }

  void ConsensusFeature::setRatios(std::vector<ConsensusFeature::Ratio> & rs)
  {
    ratios_ = rs;
  }

  std::vector<ConsensusFeature::Ratio> & ConsensusFeature::getRatios()
  {
    return ratios_;
  }

  std::vector<ConsensusFeature::Ratio> ConsensusFeature::getRatios() const
  {
    return ratios_;
  }

  Size ConsensusFeature::size() const
  {
    return handles_.size();
  }

  ConsensusFeature::const_iterator ConsensusFeature::begin() const
  {
    return handles_.begin();
  }

  ConsensusFeature::iterator ConsensusFeature::begin()
  {
    return handles_.begin();
  }


  ConsensusFeature::const_iterator ConsensusFeature::end() const
  {
    return handles_.end();
  }

  ConsensusFeature::iterator ConsensusFeature::end()
  {
    return handles_.end();
  }

  ConsensusFeature::const_reverse_iterator ConsensusFeature::rbegin() const
  {
    return handles_.rbegin();
  }

  ConsensusFeature::reverse_iterator ConsensusFeature::rbegin()
  {
    return handles_.rbegin();
  }

  ConsensusFeature::const_reverse_iterator ConsensusFeature::rend() const
  {
    return handles_.rend();
  }

  ConsensusFeature::reverse_iterator ConsensusFeature::rend()
  {
    return handles_.rend();
  }

  void ConsensusFeature::clear()
  {
    handles_.clear();
  }

  bool ConsensusFeature::empty() const
  {
    return handles_.empty();
  }


  std::ostream & operator<<(std::ostream & os, const ConsensusFeature & cons)
  {
    os << "---------- CONSENSUS ELEMENT BEGIN -----------------\n";
    os << "Position: " << cons.getPosition() << std::endl;
    os << "Intensity " << precisionWrapper(cons.getIntensity()) << std::endl;
    os << "Quality " << precisionWrapper(cons.getQuality()) << std::endl;
    os << "Grouped features: " << std::endl;

    for (ConsensusFeature::HandleSetType::const_iterator it = cons.begin(); it != cons.end(); ++it)
    {
      os << " - Map index: " << it->getMapIndex() << std::endl
        << "   Feature id: " << it->getUniqueId() << std::endl
        << "   RT: " << precisionWrapper(it->getRT()) << std::endl
        << "   m/z: " << precisionWrapper(it->getMZ()) << std::endl
        << "   Intensity: " << precisionWrapper(it->getIntensity()) << std::endl;
    }

    os << "Meta information: " << std::endl;
    std::vector<String> keys;
    cons.getKeys(keys);
    for (std::vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
      os << "   " << (*it) << ": " << cons.getMetaValue(*it) << std::endl;
    }
    os << "---------- CONSENSUS ELEMENT END ----------------- " << std::endl;

    return os;
  }

}
