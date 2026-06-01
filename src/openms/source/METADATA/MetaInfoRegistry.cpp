// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Hendrik Weisser $
// -------------------------------------------------------------------------

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
#pragma omp critical (MetaInfoRegistry)
    {
      next_index_ = rhs.next_index_;
      name_to_index_ = rhs.name_to_index_;
      index_to_name_ = rhs.index_to_name_;
      index_to_description_ = rhs.index_to_description_;
      index_to_unit_ = rhs.index_to_unit_;
    }
    return *this;
  }

  UInt MetaInfoRegistry::registerName(const String& name, const String& description, const String& unit)
  {
    UInt rv;
#pragma omp critical (MetaInfoRegistry)
    {
      MapString2IndexType::iterator it = name_to_index_.find(name);
      if (it == name_to_index_.end())
      {
        name_to_index_[name] = next_index_;
        index_to_name_[next_index_] = name;
        index_to_description_[next_index_] = description;
        index_to_unit_[next_index_] = unit;
        rv = next_index_++;
      }
      else
      {
        rv = it->second;
      }
    }
    return rv;
  }

  void MetaInfoRegistry::setDescription(UInt index, const String& description)
  {
    MapIndex2StringType::iterator pos;
#pragma omp critical (MetaInfoRegistry)
    {
      pos = index_to_description_.find(index);
      if (pos != index_to_description_.end())
      {
        pos->second = description;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!", String(index));
      }
    }
  }

  void MetaInfoRegistry::setDescription(const String& name, const String& description)
  {
    MapString2IndexType::iterator pos;
#pragma omp critical (MetaInfoRegistry)
    {
      pos = name_to_index_.find(name);
      if (pos != name_to_index_.end())
      {
        index_to_description_[pos->second] = description;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered name!", name);
      }
    }
  }

  void MetaInfoRegistry::setUnit(UInt index, const String& unit)
  {
    MapIndex2StringType::iterator pos;
#pragma omp critical (MetaInfoRegistry)
    {
      pos = index_to_unit_.find(index);
      if (pos != index_to_unit_.end())
      {
        pos->second = unit;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!", String(index));
      }
    }
  }

  void MetaInfoRegistry::setUnit(const String& name, const String& unit)
  {
    MapString2IndexType::iterator pos;
#pragma omp critical (MetaInfoRegistry)
    {
      pos = name_to_index_.find(name);
      if (pos != name_to_index_.end())
      {
        index_to_unit_[pos->second] = unit;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered name!", name);
      }
    }
  }

  UInt MetaInfoRegistry::getIndex(const String& name) const
  {
    UInt rv = UInt(-1);
#pragma omp critical (MetaInfoRegistry)
    {
      MapString2IndexType::const_iterator it = name_to_index_.find(name);
      if (it != name_to_index_.end())
      {
        rv = it->second;
      }
    }
    return rv;
  }

  String MetaInfoRegistry::getDescription(UInt index) const
  {
    String result;
#pragma omp critical (MetaInfoRegistry)
    {
      MapIndex2StringType::const_iterator it = index_to_description_.find(index);
      if (it == index_to_description_.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!", String(index));
      }
      result = it->second;
    }
    return result;
  }

  String MetaInfoRegistry::getDescription(const String& name) const
  {
    String rv;
    UInt index = getIndex(name); // this has to be outside the OpenMP "critical" block!
    if (index == UInt(-1)) // not found
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered Name!", name);
    }
    else
    {
#pragma omp critical (MetaInfoRegistry)
      {
        rv = (index_to_description_.find(index))->second;
      }
    }
    return rv;
  }

  String MetaInfoRegistry::getUnit(UInt index) const
  {
    String result;
#pragma omp critical (MetaInfoRegistry)
    {
      MapIndex2StringType::const_iterator it = index_to_unit_.find(index);
      if (it == index_to_unit_.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!", String(index));
      }
      result = it->second;
    }
    return result;
  }

  String MetaInfoRegistry::getUnit(const String& name) const
  {
    String rv;
    UInt index = getIndex(name); // this has to be outside the OpenMP "critical" block!
    if (index == UInt(-1)) // not found
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered Name!", name);
    }
    else
    {
#pragma omp critical (MetaInfoRegistry)
      {
        rv = (index_to_unit_.find(index))->second;
      }
    }
    return rv;
  }

  String MetaInfoRegistry::getName(UInt index) const
  {
    String rv;
#pragma omp critical (MetaInfoRegistry)
    {
      MapIndex2StringType::const_iterator it = index_to_name_.find(index);
      if (it != index_to_name_.end())
      {
        rv = it->second;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unregistered index!", String(index));
      }
    }
    return rv;
  }

} //namespace

