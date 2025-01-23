// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

using namespace std;

namespace OpenMS
{

  bool DataFilters::DataFilter::operator==(const DataFilter & rhs) const
  {
    return field == rhs.field
           && op == rhs.op
           && value == rhs.value
           && value_string == rhs.value_string
           && meta_name == rhs.meta_name
           && value_is_numerical == rhs.value_is_numerical;
  }

  ///Inequality operator
  bool DataFilters::DataFilter::operator!=(const DataFilter & rhs) const
  {
    return !operator==(rhs);
  }


  String DataFilters::DataFilter::toString() const
  {
    String out;
    //field
    if (field == INTENSITY)
    {
      out = "Intensity ";
    }
    else if (field == QUALITY)
    {
      out = "Quality ";
    }
    else if (field == CHARGE)
    {
      out = "Charge ";
    }
    else if (field == SIZE)
    {
      out = "Size ";
    }
    else if (field == META_DATA)
    {
      out = "Meta::" + meta_name + " ";
    }
    //operation
    if (op == GREATER_EQUAL)
    {
      out += ">= ";
    }
    else if (op == EQUAL)
    {
      out += "= ";
    }
    else if (op == LESS_EQUAL)
    {
      out += "<= ";
    }
    else if (op == EXISTS)
    {
      out += "exists";
    }
    //value
    if (field == META_DATA)
    {
      if (op != EXISTS)
      {
        if (value_is_numerical)
        {
          out = out + value;
        }
        else
        {
          out = out + "\"" + value_string + "\"";
        }
      }
      return out;
    }
    out = out + value;
    return out;
  }

  void DataFilters::DataFilter::fromString(const String & filter)
  {
    bool meta = false;
    String tmp = filter;
    tmp.trim();
    StringList parts;
    tmp.split(' ', parts);
    SignedSize size = parts.size();
    if (size < 2)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid filter format", tmp);
    }
    //field
    tmp = parts[0];
    tmp.toLower();
    if (tmp == "intensity")
    {
      field = INTENSITY;
    }
    else if (tmp == "charge")
    {
      field = CHARGE;
    }
    else if (tmp == "size")
    {
      field = SIZE;
    }
    else if (tmp == "quality")
    {
      field = QUALITY;
    }
    else if (tmp.hasPrefix(String("meta::")))
    {
      meta = true;
      field = META_DATA;
      meta_name = tmp.suffix(tmp.size() - 6);
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid field name", tmp);
    }
    //operation
    tmp = parts[1];
    if (tmp == ">=")
    {
      op = GREATER_EQUAL;
    }
    else if (tmp == "=")
    {
      op = EQUAL;
    }
    else if (tmp == "<=")
    {
      op = LESS_EQUAL;
    }
    else if (tmp == "exists" && meta)
    {
      op = EXISTS;
      return;
    }
    else
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid operator", tmp);
    //value
    if (size > 3)     // string values may contain spaces, implode to a single string
    {
      tmp.concatenate(parts.begin() + 2, parts.end(), " ");
    }
    else if (size == 3)
    {
      tmp = parts[2];
    }
    else     // size < 3 && operation is binary (only "exists" is unary) --> invalid
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid filter format", tmp);
    }
    try
    {
      value = tmp.toDouble();
      value_is_numerical = true;
    }
    catch (Exception::ConversionError&)
    {
      value_is_numerical = false;
      if (!(tmp.hasPrefix("\"") && tmp.hasSuffix("\"")))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid value", tmp);
      }
      else
      {
        tmp = tmp.substr(1, tmp.size() - 2);
      }
      if (!meta)       // non meta values must be numerical
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "invalid value", tmp);
      }
      else
      {
        value_string = tmp;
      }
    }
  }

  void DataFilters::add(const DataFilter& filter)
  {
    //activate if not empty
    is_active_ = true;

    filters_.push_back(filter);
    if (filter.field == DataFilters::META_DATA)
    {
      meta_indices_.push_back(MetaInfo::registry().getIndex(filter.meta_name));
    }
    else
    {
      meta_indices_.push_back(0);
    }
  }

  void DataFilters::remove(Size index)
  {
    if (index >= filters_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, filters_.size());
    }
    filters_.erase(filters_.begin() + index);
    meta_indices_.erase(meta_indices_.begin() + index);

    //disable if empty
    if (size() == 0)
    {
      is_active_ = false;
    }
  }

  void DataFilters::replace(Size index, const DataFilter & filter)
  {
    if (index >= filters_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, filters_.size());
    }
    filters_[index] = filter;
    if (filter.field == DataFilters::META_DATA)
    {
      meta_indices_[index] = MetaInfo::registry().getIndex(filter.meta_name);
    }
    else
    {
      meta_indices_[index] = 0;
    }
  }

  void DataFilters::clear()
  {
    filters_.clear();
    meta_indices_.clear();
    is_active_ = false;
  }

  Size DataFilters::size() const
  {
    return filters_.size();
  }

  const DataFilters::DataFilter & DataFilters::operator[](Size index) const
  {
    if (index >= filters_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, filters_.size());
    }
    return filters_[index];
  }

  bool DataFilters::passes(const Feature & feature) const
  {
    if (!is_active_)
    {
      return true;
    }
    for (Size i = 0; i < filters_.size(); i++)
    {
      const DataFilters::DataFilter & filter = filters_[i];
      if (filter.field == INTENSITY)
      {
        if (filter.op == GREATER_EQUAL && feature.getIntensity() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && feature.getIntensity() > filter.value)
        {
          return false;
        }
        else if (filter.op == EQUAL && feature.getIntensity() != filter.value)
        {
          return false;
        }
      }
      else if (filter.field == QUALITY)
      {
        if (filter.op == GREATER_EQUAL && feature.getOverallQuality() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && feature.getOverallQuality() > filter.value)
        {
          return false;
        }
        else if (filter.op == EQUAL && feature.getOverallQuality() != filter.value)
        {
          return false;
        }
      }
      else if (filter.field == CHARGE)
      {
        if (filter.op == EQUAL && feature.getCharge() != filter.value)
        {
          return false;
        }
        else if (filter.op == GREATER_EQUAL && feature.getCharge() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && feature.getCharge() > filter.value)
        {
          return false;
        }
      }
      else if (filter.field == SIZE)
      {
        if (filter.op == EQUAL && feature.getSubordinates().size() != filter.value)
        {
          return false;
        }
        else if (filter.op == GREATER_EQUAL && feature.getSubordinates().size() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && feature.getSubordinates().size() > filter.value)
        {
          return false;
        }
      }
      else if (filter.field == META_DATA)
      {
        const MetaInfoInterface & mii = static_cast<MetaInfoInterface>(feature);
        if (!metaPasses_(mii, filter, meta_indices_[i]))
        {
          return false;
        }
      }
    }
    return true;
  }

  bool DataFilters::passes(const ConsensusFeature & consensus_feature) const
  {
    if (!is_active_)
    {
      return true;
    }
    for (Size i = 0; i < filters_.size(); i++)
    {
      const DataFilters::DataFilter & filter = filters_[i];
      if (filter.field == INTENSITY)
      {
        if (filter.op == GREATER_EQUAL && consensus_feature.getIntensity() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && consensus_feature.getIntensity() > filter.value)
        {
          return false;
        }
        else if (filter.op == EQUAL && consensus_feature.getIntensity() != filter.value)
        {
          return false;
        }
      }
      else if (filter.field == QUALITY)
      {
        if (filter.op == GREATER_EQUAL && consensus_feature.getQuality() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && consensus_feature.getQuality() > filter.value)
        {
          return false;
        }
        else if (filter.op == EQUAL && consensus_feature.getQuality() != filter.value)
        {
          return false;
        }
      }
      else if (filter.field == CHARGE)
      {
        if (filter.op == EQUAL && consensus_feature.getCharge() != filter.value)
        {
          return false;
        }
        else if (filter.op == GREATER_EQUAL && consensus_feature.getCharge() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && consensus_feature.getCharge() > filter.value)
        {
          return false;
        }
      }
      else if (filter.field == SIZE)
      {
        if (filter.op == EQUAL && consensus_feature.size() != filter.value)
        {
          return false;
        }
        else if (filter.op == GREATER_EQUAL && consensus_feature.size() < filter.value)
        {
          return false;
        }
        else if (filter.op == LESS_EQUAL && consensus_feature.size() > filter.value)
        {
          return false;
        }
      }
      else if (filter.field == META_DATA)
      {
        const MetaInfoInterface & mii = static_cast<MetaInfoInterface>(consensus_feature);
        if (!metaPasses_(mii, filter, meta_indices_[i]))
        {
          return false;
        }
      }
    }
    return true;
  }

  bool DataFilters::metaPasses_(const MetaInfoInterface& meta_interface, const DataFilters::DataFilter& filter, Size index) const
  {
    if (!meta_interface.metaValueExists((UInt)index))
    {
      return false;
    }
    else if (filter.op != EXISTS)
    {
      const DataValue& data_value = meta_interface.getMetaValue((UInt)index);
      if (!filter.value_is_numerical)
      {
        if (data_value.valueType() != DataValue::STRING_VALUE)
        {
          return false;
        }
        else
        {
          // for string values, equality is the only valid operation (besides "exists", see above)
          if (filter.op != EQUAL)
          {
            return false;
          }
          else if (filter.value_string != data_value.toString())
          {
            return false;
          }
        }
      }
      else             // value_is_numerical
      {
        if (data_value.valueType() == DataValue::STRING_VALUE || data_value.valueType() == DataValue::EMPTY_VALUE)
        {
          return false;
        }
        else
        {
          if (filter.op == EQUAL && (double)data_value != filter.value)
          {
            return false;
          }
          else if (filter.op == LESS_EQUAL && (double)data_value > filter.value)
          {
            return false;
          }
          else if (filter.op == GREATER_EQUAL && (double)data_value < filter.value)
          {
            return false;
          }
        }
      }
    }
    return true;
  }


  void DataFilters::setActive(bool is_active)
  {
    is_active_ = is_active;
  }

} //Namespace
