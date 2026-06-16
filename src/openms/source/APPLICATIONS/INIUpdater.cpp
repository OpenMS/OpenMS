// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/INIUpdater.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <algorithm>

namespace OpenMS
{
  using namespace Internal;

  typedef ToolDescriptionInternal TDE;

  INIUpdater::INIUpdater()
  {
    getNameMapping(); // init map_
  }

  StringList INIUpdater::getToolNamesFromINI(const Param & ini) const
  {
    StringList tool_names;

    for (Param::ParamIterator it = ini.begin(); it != ini.end(); ++it)
    {
      std::string name = it.getName();
      if (std::count(name.begin(), name.end(), ':') == 1 && StringUtils::hasSuffix(name, ":version"))
      {
        tool_names.push_back(StringUtils::prefix(name, ':'));
      }
    }
    return tool_names;
  }

  const ToolMapping & INIUpdater::getNameMapping()
  {
    if (map_.empty())
    {
      map_[TDE("FeatureFinder", ListUtils::create<std::string>("centroided"))] = TDE("FeatureFinderCentroided", ListUtils::create<std::string>(""));
      map_[TDE("FeatureFinder", ListUtils::create<std::string>("isotope_wavelet"))] = TDE("FeatureFinderIsotopeWavelet", ListUtils::create<std::string>(""));
      map_[TDE("FeatureFinder", ListUtils::create<std::string>("mrm"))] = TDE("FeatureFinderMRM", ListUtils::create<std::string>(""));

      map_[TDE("FeatureLinker", ListUtils::create<std::string>("labeled"))] = TDE("FeatureLinkerLabeled", ListUtils::create<std::string>(""));
      map_[TDE("FeatureLinker", ListUtils::create<std::string>("unlabeled"))] = TDE("FeatureLinkerUnlabeled", ListUtils::create<std::string>(""));
      map_[TDE("FeatureLinker", ListUtils::create<std::string>("unlabeled_qt"))] = TDE("FeatureLinkerUnlabeledQT", ListUtils::create<std::string>(""));

      map_[TDE("NoiseFilter", ListUtils::create<std::string>("gaussian"))] = TDE("NoiseFilterGaussian", ListUtils::create<std::string>(""));
      map_[TDE("NoiseFilter", ListUtils::create<std::string>("sgolay"))] = TDE("NoiseFilterSGolay", ListUtils::create<std::string>(""));

      map_[TDE("MapAligner", ListUtils::create<std::string>("apply_given_trafo"))] = TDE("MapRTTransformer", ListUtils::create<std::string>(""));
      map_[TDE("MapAligner", ListUtils::create<std::string>("identification"))] = TDE("MapAlignerIdentification", ListUtils::create<std::string>(""));
      map_[TDE("MapAligner", ListUtils::create<std::string>("pose_clustering"))] = TDE("MapAlignerPoseClustering", ListUtils::create<std::string>(""));

      // SpectraFilter ...

      map_[TDE("PeakPicker", ListUtils::create<std::string>("wavelet"))] = TDE("PeakPickerWavelet", ListUtils::create<std::string>(""));
      map_[TDE("PeakPicker", ListUtils::create<std::string>("high_res"))] = TDE("PeakPickerHiRes", ListUtils::create<std::string>(""));
    }

    return map_;
  }

  bool INIUpdater::getNewToolName(const std::string & old_name, const std::string & tools_type, std::string & new_name)
  {
    new_name = "";
    // try with type (as some new tools for one type might have the exact same name as old ones with several types)
    TDE old_withtype(old_name, ListUtils::create<std::string>(tools_type));
    if (map_.contains(old_withtype))
    {
      new_name = map_[old_withtype].name;
      return true;
    }

    // try without type
    TDE old_notype(old_name, StringList());
    if (map_.contains(old_notype))
    {
      new_name = map_[old_notype].name;
      return true;
    }

    // default to ToolHandler
    const auto& topp = ToolHandler::getTOPPToolList();
    if (topp.contains(old_name))
    {
      new_name = old_name;
      return true;
    }

    return false;
  }

  ToolMapping INIUpdater::map_;



} // namespace
