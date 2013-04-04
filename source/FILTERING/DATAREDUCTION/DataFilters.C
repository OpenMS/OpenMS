// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  String DataFilters::DataFilter::toString() const
  {
    String out;
    //field
    if (field == INTENSITY)
      out = "Intensity ";
    else if (field == QUALITY)
      out = "Quality ";
    else if (field == CHARGE)
      out = "Charge ";
    else if (field == SIZE)
      out = "Size ";
    else if (field == META_DATA)
      out = "Meta::" + meta_name + " ";
    //operation
    if (op == GREATER_EQUAL)
      out += ">= ";
    else if (op == EQUAL)
      out += "= ";
    else if (op == LESS_EQUAL)
      out += "<= ";
    else if (op == EXISTS)
      out += "exists";
    //value
    if (field == META_DATA)
    {
      if (op != EXISTS)
      {
        if (value_is_numerical)
          out = out + value;
        else
          out = out + "\"" + value_string + "\"";
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
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid filter format.", tmp);
    //field
    tmp = parts[0];
    tmp.toLower();
    if (tmp == "intensity")
      field = INTENSITY;
    else if (tmp == "charge")
      field = CHARGE;
    else if (tmp == "size")
      field = SIZE;
    else if (tmp == "quality")
      field = QUALITY;
    else if (tmp.hasPrefix(String("meta::")))
    {
      meta = true;
      field = META_DATA;
      meta_name = tmp.suffix(tmp.size() - 6);
    }
    else
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid field name.", tmp);
    //operation
    tmp = parts[1];
    if (tmp == ">=")
      op = GREATER_EQUAL;
    else if (tmp == "=")
      op = EQUAL;
    else if (tmp == "<=")
      op = LESS_EQUAL;
    else if (tmp == "exists" && meta)
    {
      op = EXISTS;
      return;
    }
    else
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid operator.", tmp);
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
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid filter format.", tmp);
    }
    try
    {
      value = tmp.toDouble();
      value_is_numerical = true;
    }
    catch (Exception::ConversionError)
    {
      value_is_numerical = false;
      if (!(tmp.hasPrefix("\"") && tmp.hasSuffix("\"")))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid value.", tmp);
      }
      else
      {
        tmp = tmp.substr(1, tmp.size() - 2);
      }
      if (!meta)       // non meta values must be numerical
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid value.", tmp);
      }
      else
      {
        value_string = tmp;
      }
    }
  }

  void DataFilters::add(const DataFilter & filter)
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, filters_.size());
    filters_.erase(filters_.begin() + index);
    meta_indices_.erase(meta_indices_.begin() + index);

    //disable if empty
    if (size() == 0)
      is_active_ = false;
  }

  void DataFilters::replace(Size index, const DataFilter & filter)
  {
    if (index >= filters_.size())
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, filters_.size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, filters_.size());
    return filters_[index];
  }

  bool DataFilters::passes(const Feature & feature) const
  {
    if (!is_active_)
      return true;

    for (Size i = 0; i < filters_.size(); i++)
    {
      const DataFilters::DataFilter & filter = filters_[i];
      if (filter.field == INTENSITY)
      {
        if (filter.op == GREATER_EQUAL && feature.getIntensity() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && feature.getIntensity() > filter.value)
          return false;
        else if (filter.op == EQUAL && feature.getIntensity() != filter.value)
          return false;
      }
      else if (filter.field == QUALITY)
      {
        if (filter.op == GREATER_EQUAL && feature.getOverallQuality() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && feature.getOverallQuality() > filter.value)
          return false;
        else if (filter.op == EQUAL && feature.getOverallQuality() != filter.value)
          return false;
      }
      else if (filter.field == CHARGE)
      {
        if (filter.op == EQUAL && feature.getCharge() != filter.value)
          return false;
        else if (filter.op == GREATER_EQUAL && feature.getCharge() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && feature.getCharge() > filter.value)
          return false;
      }
      else if (filter.field == SIZE)
      {
        if (filter.op == EQUAL && feature.getSubordinates().size() != filter.value)
          return false;
        else if (filter.op == GREATER_EQUAL && feature.getSubordinates().size() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && feature.getSubordinates().size() > filter.value)
          return false;
      }
      else if (filter.field == META_DATA)
      {
        const MetaInfoInterface & mii = static_cast<MetaInfoInterface>(feature);
        if (!metaPasses_(mii, filter, meta_indices_[i]))
          return false;
      }
    }
    return true;
  }

  bool DataFilters::passes(const ConsensusFeature & consensus_feature) const
  {
    if (!is_active_)
      return true;

    for (Size i = 0; i < filters_.size(); i++)
    {
      const DataFilters::DataFilter & filter = filters_[i];
      if (filter.field == INTENSITY)
      {
        if (filter.op == GREATER_EQUAL && consensus_feature.getIntensity() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && consensus_feature.getIntensity() > filter.value)
          return false;
        else if (filter.op == EQUAL && consensus_feature.getIntensity() != filter.value)
          return false;
      }
      else if (filter.field == QUALITY)
      {
        if (filter.op == GREATER_EQUAL && consensus_feature.getQuality() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && consensus_feature.getQuality() > filter.value)
          return false;
        else if (filter.op == EQUAL && consensus_feature.getQuality() != filter.value)
          return false;
      }
      else if (filter.field == CHARGE)
      {
        if (filter.op == EQUAL && consensus_feature.getCharge() != filter.value)
          return false;
        else if (filter.op == GREATER_EQUAL && consensus_feature.getCharge() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && consensus_feature.getCharge() > filter.value)
          return false;
      }
      else if (filter.field == SIZE)
      {
        if (filter.op == EQUAL && consensus_feature.size() != filter.value)
          return false;
        else if (filter.op == GREATER_EQUAL && consensus_feature.size() < filter.value)
          return false;
        else if (filter.op == LESS_EQUAL && consensus_feature.size() > filter.value)
          return false;
      }
      else if (filter.field == META_DATA)
      {
        const MetaInfoInterface & mii = static_cast<MetaInfoInterface>(consensus_feature);
        if (!metaPasses_(mii, filter, meta_indices_[i]))
          return false;
      }
    }
    return true;
  }

  void DataFilters::setActive(bool is_active)
  {
    is_active_ = is_active;
  }

} //Namespace
