// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{
  ConsensusFeature::ConsensusFeature() :
    BaseFeature(), handles_(), ratios_()
  {
  }

  ConsensusFeature::ConsensusFeature(const BaseFeature& feature) :
    BaseFeature(feature), handles_(), ratios_()
  {
  }

  ConsensusFeature::ConsensusFeature(UInt64 map_index, const Peak2D& element, UInt64 element_index) :
    BaseFeature(element), handles_(), ratios_()
  {
    insert(map_index, element, element_index);
  }

  ConsensusFeature::ConsensusFeature(UInt64 map_index, const BaseFeature& element) :
    BaseFeature(element, map_index), handles_(), ratios_()
  {
    insert(FeatureHandle(map_index, element));
  }

  ConsensusFeature::~ConsensusFeature() = default;

  void ConsensusFeature::insert(const ConsensusFeature& cf)
  {
    handles_.insert(cf.handles_.begin(), cf.handles_.end());
    peptides_.insert(peptides_.end(), cf.getPeptideIdentifications().begin(), cf.getPeptideIdentifications().end());
  }

  void ConsensusFeature::insert(ConsensusFeature&& cf)
  {
    handles_.insert(make_move_iterator(cf.handles_.begin()), make_move_iterator(cf.handles_.end()));
    peptides_.insert(peptides_.end(), make_move_iterator(cf.getPeptideIdentifications().begin()), make_move_iterator(cf.getPeptideIdentifications().end()));
  }

  void ConsensusFeature::insert(const FeatureHandle& handle)
  {
    if (!(handles_.insert(handle).second))
    {
      String key = String("map") + handle.getMapIndex() + "/feature" + handle.getUniqueId();
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The set already contained an element with this key.", key);
    }
  }

  void ConsensusFeature::insert(FeatureHandle&& handle)
  {
    if (!(handles_.insert(std::move(handle)).second))
    {
      String key = String("map") + handle.getMapIndex() + "/feature" + handle.getUniqueId();
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The set already contained an element with this key.", key);
    }
  }

  void ConsensusFeature::insert(const HandleSetType& handle_set)
  {
    for (auto& handle : handle_set)
    {
      insert(handle);
    }
  }

  void ConsensusFeature::insert(HandleSetType&& handle_set)
  {
    for (auto& handle : handle_set)
    {
      insert(std::move(handle));
    }
  }

  void ConsensusFeature::insert(UInt64 map_index, const Peak2D& element, UInt64 element_index)
  {
    insert(FeatureHandle(map_index, element, element_index));
  }

  void ConsensusFeature::insert(UInt64 map_index, const BaseFeature& element)
  {
    insert(FeatureHandle(map_index, element));
    // annotate map index to peptide identification
    std::vector<PeptideIdentification> ids(element.getPeptideIdentifications());
    for (PeptideIdentification& it : ids)
    {
      it.setMetaValue("map_index", map_index);
    }
    peptides_.insert(peptides_.end(), ids.begin(), ids.end());
  }

  void ConsensusFeature::setFeatures(HandleSetType h)
  {
    handles_ = std::move(h);
  }

  const ConsensusFeature::HandleSetType& ConsensusFeature::getFeatures() const
  {
    return handles_;
  }

  std::vector<FeatureHandle> ConsensusFeature::getFeatureList() const
  {
    std::vector<FeatureHandle> tmp;
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      tmp.push_back(*it);
    }
    return tmp;
  }

  DRange<2> ConsensusFeature::getPositionRange() const
  {
    DPosition<2> min = DPosition<2>::maxPositive();
    DPosition<2> max = DPosition<2>::minPositive();
    for (ConsensusFeature::HandleSetType::const_iterator it = handles_.begin(); it != handles_.end(); ++it)
    {
      if (it->getRT() < min[0])
      {
        min[0] = it->getRT();
      }
      if (it->getRT() > max[0])
      {
        max[0] = it->getRT();
      }
      if (it->getMZ() < min[1])
      {
        min[1] = it->getMZ();
      }
      if (it->getMZ() > max[1])
      {
        max[1] = it->getMZ();
      }
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
      {
        min[0] = it->getIntensity();
      }
      if (it->getIntensity() > max[0])
      {
        max[0] = it->getIntensity();
      }
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
      {
        mz = it->getMZ();
      }
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

  void ConsensusFeature::computeDechargeConsensus(const FeatureMap& fm, bool intensity_weighted_averaging)
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
      {
        OPENMS_LOG_WARN << "ConsensusFeature::computeDechargeConsensus() WARNING: Feature's charge is 0! This will lead to M=0!\n";
      }
      double adduct_mass;
      Size index = fm.uniqueIdToIndex(it->getUniqueId());
      if (index > fm.size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, fm.size());
      }
      if (fm[index].metaValueExists("dc_charge_adduct_mass"))
      {
        adduct_mass = (double) fm[index].getMetaValue("dc_charge_adduct_mass");
      }
      else
      {
        adduct_mass = q * proton_mass;
      }

      if (intensity_weighted_averaging)
      {
        weighting_factor = it->getIntensity() / intensity;
      }
      rt += it->getRT() * weighting_factor;
      m += (it->getMZ() * abs(q) - adduct_mass) * weighting_factor;
    }

    // compute the average position and intensity
    setRT(rt);
    setMZ(m);
    setIntensity(intensity);
    setCharge(0);
    return;
  }

  void ConsensusFeature::addRatio(const ConsensusFeature::Ratio& r)
  {
    ratios_.push_back(r);
  }

  void ConsensusFeature::setRatios(std::vector<ConsensusFeature::Ratio>& rs)
  {
    ratios_ = rs;
  }

  std::vector<ConsensusFeature::Ratio>& ConsensusFeature::getRatios()
  {
    return ratios_;
  }

  const std::vector<ConsensusFeature::Ratio>& ConsensusFeature::getRatios() const
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

  std::ostream& operator<<(std::ostream& os, const ConsensusFeature& cons)
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
    for (const String& it : keys)
    {
      os << "   " << (it) << ": " << cons.getMetaValue(it) << std::endl;
    }
    os << "---------- CONSENSUS ELEMENT END ----------------- " << std::endl;

    return os;
  }

}
