// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/INIUpdater.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <QtCore/QString>

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
      String name = it.getName();
      if (name.toQString().count(':') == 1 && name.hasSuffix(":version"))
      {
        tool_names.push_back(name.prefix(':'));
      }
    }
    return tool_names;
  }

  const ToolMapping & INIUpdater::getNameMapping()
  {
    if (map_.empty())
    {
      map_[TDE("FeatureFinder", ListUtils::create<String>("centroided"))] = TDE("FeatureFinderCentroided", ListUtils::create<String>(""));
      map_[TDE("FeatureFinder", ListUtils::create<String>("isotope_wavelet"))] = TDE("FeatureFinderIsotopeWavelet", ListUtils::create<String>(""));
      map_[TDE("FeatureFinder", ListUtils::create<String>("mrm"))] = TDE("FeatureFinderMRM", ListUtils::create<String>(""));

      map_[TDE("FeatureLinker", ListUtils::create<String>("labeled"))] = TDE("FeatureLinkerLabeled", ListUtils::create<String>(""));
      map_[TDE("FeatureLinker", ListUtils::create<String>("unlabeled"))] = TDE("FeatureLinkerUnlabeled", ListUtils::create<String>(""));
      map_[TDE("FeatureLinker", ListUtils::create<String>("unlabeled_qt"))] = TDE("FeatureLinkerUnlabeledQT", ListUtils::create<String>(""));

      map_[TDE("NoiseFilter", ListUtils::create<String>("gaussian"))] = TDE("NoiseFilterGaussian", ListUtils::create<String>(""));
      map_[TDE("NoiseFilter", ListUtils::create<String>("sgolay"))] = TDE("NoiseFilterSGolay", ListUtils::create<String>(""));

      map_[TDE("MapAligner", ListUtils::create<String>("apply_given_trafo"))] = TDE("MapRTTransformer", ListUtils::create<String>(""));
      map_[TDE("MapAligner", ListUtils::create<String>("identification"))] = TDE("MapAlignerIdentification", ListUtils::create<String>(""));
      map_[TDE("MapAligner", ListUtils::create<String>("pose_clustering"))] = TDE("MapAlignerPoseClustering", ListUtils::create<String>(""));

      // SpectraFilter ...

      map_[TDE("PeakPicker", ListUtils::create<String>("wavelet"))] = TDE("PeakPickerWavelet", ListUtils::create<String>(""));
      map_[TDE("PeakPicker", ListUtils::create<String>("high_res"))] = TDE("PeakPickerHiRes", ListUtils::create<String>(""));
    }

    return map_;
  }

  bool INIUpdater::getNewToolName(const String & old_name, const String & tools_type, String & new_name)
  {
    new_name = "";
    // try with type (as some new tools for one type might have the exact same name as old ones with several types)
    TDE old_withtype(old_name, ListUtils::create<String>(tools_type));
    if (map_.find(old_withtype) != map_.end())
    {
      new_name = map_[old_withtype].name;
      return true;
    }

    // try without type
    TDE old_notype(old_name, StringList());
    if (map_.find(old_notype) != map_.end())
    {
      new_name = map_[old_notype].name;
      return true;
    }

    // default to ToolHandler
    const auto& topp = ToolHandler::getTOPPToolList(true);
    if (topp.find(old_name) != topp.end())
    {
      new_name = old_name;
      return true;
    }

    return false;
  }

  ToolMapping INIUpdater::map_;



} // namespace
