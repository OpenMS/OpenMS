// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <algorithm>
#include <map>
#include <vector>

namespace OpenMS
{
  namespace Detail
  {
    template<typename T>
    struct MetaKeyGetter
    {
      static void getKeys(const T& object, std::vector<String>& keys)
      {
        object.getKeys(keys);
      };
    };
  }

  /**
    @brief Utilities operating on containers inheriting from MetaInfoInterface

    @ingroup Metadata
  */
  class /*OPENMS_DLLAPI -- disabled since it's template code only */ MetaInfoInterfaceUtils
  {
public:
    /// hide c'tors to avoid instantiation of utils class
    MetaInfoInterfaceUtils() = delete;
    MetaInfoInterfaceUtils(const MetaInfoInterfaceUtils&) = delete;
    MetaInfoInterfaceUtils& operator=(MetaInfoInterfaceUtils&) = delete;
    // no Move semantics for utils class


    ///@name Methods to find key sets
    //@{
    /**
      @brief Find keys in a collection of MetaInfoInterface objects which reach a certain frequency threshold.

      Searches the given iterator range for the keys of each element's MetaInfoInterface keys and returns those keys, which
      reach a certain frequency threshold. Common use cases 
      are @p min_frequency = 0 (i.e. take any key which occurs)
      and @p min_frequency = 100 (i.e. take only keys which are common to all elements in the iterator range).

      @tparam T_In Input container (e.g. std::vector or alike), containing objects which implement the MetaInfoInterface (i.e. support 'getKeys()')
      @tparam T_Out Output container of type T<String> (e.g. std::set<String>)
      @param it_start Iterator pointing to the initial position to search. (note: this does not need to correspond to the beginning of the container)
      @param it_end Iterator pointing to the end final position to search.
      @param min_frequency Minimum required frequency (in percent). Must be between 0-100. Other values are corrected to the closest value allowed.
      @param getter Helper class, which has a getKeys() member, which can extract the keys for a given MetaInfoInterface-object; see MetaKeyGetter
      @return Returns a vector/list/set of keys passing the frequency criterion.
    */
    template<typename T_In, typename T_Out>
    static T_Out findCommonMetaKeys(const typename T_In::const_iterator& it_start, const typename T_In::const_iterator& it_end, float min_frequency, typename Detail::MetaKeyGetter<typename T_In::value_type> getter = Detail::MetaKeyGetter<typename T_In::value_type>())
    {
      // make sure min_frequency is within [0,100]
      min_frequency = std::min(100.0f, std::max(0.0f, min_frequency));

      std::map<String, UInt> counter;
      typedef std::vector<String> KeysType;
      KeysType keys;
      for (typename T_In::const_iterator it = it_start; it != it_end; ++it)
      {
        getter.getKeys(*it, keys);
        for (KeysType::const_iterator itk = keys.begin(); itk != keys.end(); ++itk)
        {
          ++counter[*itk];
        }
      }
      // pick the keys which occur often enough
      const UInt required_counts = UInt(min_frequency / 100.0 * std::distance(it_start, it_end));
      T_Out common_keys;
      for (const auto& [key, count] : counter)
      {
        if (count >= required_counts) 
        {
          common_keys.insert(common_keys.end(), key);
        }
      }
      return common_keys;
    }
  
  }; // class

} // namespace OPENMS

