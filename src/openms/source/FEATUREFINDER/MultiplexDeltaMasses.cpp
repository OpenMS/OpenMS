// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>

#include <sstream>
#include <utility>

using namespace std;

namespace OpenMS
{

  MultiplexDeltaMasses::DeltaMass::DeltaMass(double dm, LabelSet ls) :
    delta_mass(dm), label_set(std::move(ls))
  {
  }
  
  MultiplexDeltaMasses::DeltaMass::DeltaMass(double dm, const String& l) :
    delta_mass(dm), label_set()
  {
    label_set.insert(l);
  }

  MultiplexDeltaMasses::MultiplexDeltaMasses() = default;

  MultiplexDeltaMasses::MultiplexDeltaMasses(const vector<MultiplexDeltaMasses::DeltaMass>& dm) :
    delta_masses_(dm)
  {
  }

  std::vector<MultiplexDeltaMasses::DeltaMass>& MultiplexDeltaMasses::getDeltaMasses()
  {
    return delta_masses_;
  }

  const std::vector<MultiplexDeltaMasses::DeltaMass>& MultiplexDeltaMasses::getDeltaMasses() const
  {
    return delta_masses_;
  }
  
  String MultiplexDeltaMasses::labelSetToString(const MultiplexDeltaMasses::LabelSet& ls)
  {
    std::stringstream ss;
    
    for (MultiplexDeltaMasses::LabelSet::const_iterator it = ls.begin(); it != ls.end(); ++it)
    {
      if (it != ls.begin())
      {
        ss << " ";
      }
      ss << (*it);
    }

    return String(ss.str());
  }

  bool operator<(const MultiplexDeltaMasses &dm1, const MultiplexDeltaMasses &dm2)
  {
    if (dm1.getDeltaMasses().size() != dm2.getDeltaMasses().size())
    {
      // Search first for complete multiplets, then knock-out cases.
      return (dm1.getDeltaMasses().size() > dm2.getDeltaMasses().size());
    }
    else
    {
      for (unsigned i = 0; i < dm1.getDeltaMasses().size(); ++i)
      {
        double ms1 = dm1.getDeltaMasses()[i].delta_mass - dm1.getDeltaMasses()[0].delta_mass;
        double ms2 = dm2.getDeltaMasses()[i].delta_mass - dm2.getDeltaMasses()[0].delta_mass;
        
        if (ms1 != ms2)
        {
          // Search first for cases without miscleavages.
          return (ms1 < ms2);
        }
      }
    }

    return (false);
  }
  
}
