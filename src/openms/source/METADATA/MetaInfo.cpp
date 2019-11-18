// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{

  MetaInfoRegistry MetaInfo::registry_ = MetaInfoRegistry();

  MetaInfo::~MetaInfo()
  {
  }

  bool MetaInfo::operator==(const MetaInfo& rhs) const
  {
    return index_to_value_ == rhs.index_to_value_;
  }

  bool MetaInfo::operator!=(const MetaInfo& rhs) const
  {
    return !(operator==(rhs));
  }

  const DataValue& MetaInfo::getValue(const String& name, const DataValue& default_value) const
  {
    MapType::const_iterator it = index_to_value_.find(registry_.getIndex(name));
    if (it != index_to_value_.end())
    {
      return it->second;
    }
    return default_value;
  }

  const DataValue& MetaInfo::getValue(UInt index, const DataValue& default_value) const
  {
    MapType::const_iterator it = index_to_value_.find(index);
    if (it != index_to_value_.end())
    {
      return it->second;
    }
    return default_value;
  }

  void MetaInfo::setValue(const String& name, const DataValue& value)
  {
    UInt index = registry_.registerName(name); // no-op if name is already registered
    setValue(index, value);
  }

  void MetaInfo::setValue(UInt index, const DataValue& value)
  {
    // @TODO: check if that index is registered in MetaInfoRegistry?
    auto it = index_to_value_.find(index);
    if (it != index_to_value_.end())
    {
      it->second = value;
    }
    else
    {
      // Note; we need to create a copy of data value here and can't use the const &
      // The underlying flat_map invalidates references to it if inserting
      // an element leads to relocation (e.g, in constructs like: m.insert(1, m[2]));)
      DataValue tmp = value;
      index_to_value_.insert(std::make_pair(index, tmp));
    }
  }

  MetaInfoRegistry& MetaInfo::registry()
  {
    return registry_;
  }

  bool MetaInfo::exists(const String& name) const
  {
    UInt index = registry_.getIndex(name);
    if (index != UInt(-1))
    {
      return (index_to_value_.find(index) != index_to_value_.end());
    }
    return false;
  }

  bool MetaInfo::exists(UInt index) const
  {
    return (index_to_value_.find(index) != index_to_value_.end());
  }

  void MetaInfo::removeValue(const String& name)
  {
    MapType::iterator it = index_to_value_.find(registry_.getIndex(name));
    if (it != index_to_value_.end())
    {
      index_to_value_.erase(it);
    }
  }

  void MetaInfo::removeValue(UInt index)
  {
    MapType::iterator it = index_to_value_.find(index);
    if (it != index_to_value_.end())
    {
      index_to_value_.erase(it);
    }
  }

  void MetaInfo::getKeys(vector<String>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i = 0;
    for (MapType::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
    {
      keys[i++] = registry_.getName(it->first);
    }
  }

  void MetaInfo::getKeys(vector<UInt>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i = 0;
    for (MapType::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
    {
      keys[i++] = it->first;
    }
  }

  bool MetaInfo::empty() const
  {
    return index_to_value_.empty();
  }

  void MetaInfo::clear()
  {
    index_to_value_.clear();
  }

} //namespace

