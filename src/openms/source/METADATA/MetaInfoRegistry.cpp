// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Hendrik Weisser $
// -------------------------------------------------------------------------

#include <mutex>
#include <sstream>

#include <OpenMS/METADATA/MetaInfoRegistry.h>

using namespace std;

namespace OpenMS
{

  MetaInfoRegistry::MetaInfoRegistry() :
    next_index_(1024), 
    name_to_index_(), 
    index_to_name_(), 
    index_to_description_(), 
    index_to_unit_()
  {
    name_to_index_["isotopic_range"] = 1;
    index_to_name_[1] = "isotopic_range";
    index_to_description_[1] = "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak";
    index_to_unit_[1] = "";

    name_to_index_["cluster_id"] = 2;
    index_to_name_[2] = "cluster_id";
    index_to_description_[2] = "consecutive numbering of isotope clusters in a spectrum";
    index_to_unit_[2] = "";

    name_to_index_["label"] = 3;
    index_to_name_[3] = "label";
    index_to_description_[3] = "label e.g. shown in visualization";
    index_to_unit_[3] = "";

    name_to_index_["icon"] = 4;
    index_to_name_[4] = "icon";
    index_to_description_[4] = "icon shown in visualization";
    index_to_unit_[4] = "";

    name_to_index_["color"] = 5;
    index_to_name_[5] = "color";
    index_to_description_[5] = "color used for visualization e.g. #FF00FF for purple";
    index_to_unit_[5] = "";

    name_to_index_["RT"] = 6;
    index_to_name_[6] = "RT";
    index_to_description_[6] = "the retention time of an identification";
    index_to_unit_[6] = "";

    name_to_index_["MZ"] = 7;
    index_to_name_[7] = "MZ";
    index_to_description_[7] = "the MZ of an identification";
    index_to_unit_[7] = "";

    name_to_index_["predicted_RT"] = 8;
    index_to_name_[8] = "predicted_RT";
    index_to_description_[8] = "the predicted retention time of a peptide hit";
    index_to_unit_[8] = "";

    name_to_index_["predicted_RT_p_value"] = 9;
    index_to_name_[9] = "predicted_RT_p_value";
    index_to_description_[9] = "the predicted RT p-value of a peptide hit";
    index_to_unit_[9] = "";

    name_to_index_["spectrum_reference"] = 10;
    index_to_name_[10] = "spectrum_reference";
    index_to_description_[10] = "Reference to a spectrum or feature number";
    index_to_unit_[10] = "";

    name_to_index_["ID"] = 11;
    index_to_name_[11] = "ID";
    index_to_description_[11] = "Some type of identifier";
    index_to_unit_[11] = "";

    name_to_index_["low_quality"] = 12;
    index_to_name_[12] = "low_quality";
    index_to_description_[12] = "Flag which indicates that some entity has a low quality (e.g. a feature pair)";
    index_to_unit_[12] = "";

    name_to_index_["charge"] = 13;
    index_to_name_[13] = "charge";
    index_to_description_[13] = "Charge of a feature or peak";
    index_to_unit_[13] = "";
  }

  MetaInfoRegistry::MetaInfoRegistry(const MetaInfoRegistry& rhs)
  {
    *this = rhs;
  }

  MetaInfoRegistry::~MetaInfoRegistry() = default;

  MetaInfoRegistry& MetaInfoRegistry::operator=(const MetaInfoRegistry& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }
    std::unique_lock<std::shared_mutex> lock(mutex_, std::defer_lock);
    std::shared_lock<std::shared_mutex> rhs_lock(rhs.mutex_, std::defer_lock);
    std::lock(lock, rhs_lock);
    next_index_ = rhs.next_index_;
    name_to_index_ = rhs.name_to_index_;
    index_to_name_ = rhs.index_to_name_;
    index_to_description_ = rhs.index_to_description_;
    index_to_unit_ = rhs.index_to_unit_;
    return *this;
  }

  UInt MetaInfoRegistry::registerName(const std::string& name, const std::string& description, const std::string& unit)
  {
    // Nearly all calls are for names that are registered already: look those up under the shared lock.
    {
      std::shared_lock<std::shared_mutex> lock(mutex_);
      MapString2IndexType::const_iterator it = name_to_index_.find(name);
      if (it != name_to_index_.end())
      {
        return it->second;
      }
    }
    std::unique_lock<std::shared_mutex> lock(mutex_);
    // another thread may have registered the name in between
    MapString2IndexType::const_iterator it = name_to_index_.find(name);
    if (it != name_to_index_.end())
    {
      return it->second;
    }
    name_to_index_[name] = next_index_;
    index_to_name_[next_index_] = name;
    index_to_description_[next_index_] = description;
    index_to_unit_[next_index_] = unit;
    return next_index_++;
  }

  void MetaInfoRegistry::setDescription(UInt index, const std::string& description)
  {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    MapIndex2StringType::iterator pos = index_to_description_.find(index);
    if (pos == index_to_description_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!",StringUtils::toStr(index));
    }
    pos->second = description;
  }

  void MetaInfoRegistry::setDescription(const std::string& name, const std::string& description)
  {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    MapString2IndexType::iterator pos = name_to_index_.find(name);
    if (pos == name_to_index_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered name!", name);
    }
    index_to_description_[pos->second] = description;
  }

  void MetaInfoRegistry::setUnit(UInt index, const std::string& unit)
  {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    MapIndex2StringType::iterator pos = index_to_unit_.find(index);
    if (pos == index_to_unit_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!",StringUtils::toStr(index));
    }
    pos->second = unit;
  }

  void MetaInfoRegistry::setUnit(const std::string& name, const std::string& unit)
  {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    MapString2IndexType::iterator pos = name_to_index_.find(name);
    if (pos == name_to_index_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered name!", name);
    }
    index_to_unit_[pos->second] = unit;
  }

  UInt MetaInfoRegistry::getIndex(const std::string& name) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapString2IndexType::const_iterator it = name_to_index_.find(name);
    return it != name_to_index_.end() ? it->second : UInt(-1);
  }

  std::string MetaInfoRegistry::getDescription(UInt index) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapIndex2StringType::const_iterator it = index_to_description_.find(index);
    if (it == index_to_description_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!",StringUtils::toStr(index));
    }
    return it->second;
  }

  std::string MetaInfoRegistry::getDescription(const std::string& name) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapString2IndexType::const_iterator it = name_to_index_.find(name);
    if (it == name_to_index_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered Name!", name);
    }
    return index_to_description_.find(it->second)->second;
  }

  std::string MetaInfoRegistry::getUnit(UInt index) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapIndex2StringType::const_iterator it = index_to_unit_.find(index);
    if (it == index_to_unit_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!",StringUtils::toStr(index));
    }
    return it->second;
  }

  std::string MetaInfoRegistry::getUnit(const std::string& name) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapString2IndexType::const_iterator it = name_to_index_.find(name);
    if (it == name_to_index_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered Name!", name);
    }
    return index_to_unit_.find(it->second)->second;
  }

  std::string MetaInfoRegistry::getName(UInt index) const
  {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    MapIndex2StringType::const_iterator it = index_to_name_.find(index);
    if (it == index_to_name_.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!",StringUtils::toStr(index));
    }
    return it->second;
  }

} //namespace

