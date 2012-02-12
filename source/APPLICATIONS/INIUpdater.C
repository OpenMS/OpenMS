// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/INIUpdater.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>

namespace OpenMS
{
  using namespace Internal;

  typedef ToolDescriptionInternal TDE;

  INIUpdater::INIUpdater()
  {
    getNameMapping(); // init map_
  }

  StringList INIUpdater::getToolNamesFromINI(const Param& ini) const
  {
    StringList tool_names;
    std::set <String> tool_names_set;
    for (Param::ParamIterator it = ini.begin(); it!=ini.end(); ++it)
    {
      if (it.getName().toQString().count(':')==1 && it.getName().hasSuffix(":version"))
      {
        tool_names.push_back(it.getName().prefix(':'));
      }
    }
    return tool_names;
  }

  const ToolMapping& INIUpdater::getNameMapping()
  {
    if (map_.empty())
    {
      map_[TDE("FeatureFinder",StringList::create("centroided"))] = TDE("FeatureFinderCentroided",StringList::create(""));
      map_[TDE("FeatureFinder",StringList::create("isotope_wavelet"))] = TDE("FeatureFinderIsotopeWavelet",StringList::create(""));
      map_[TDE("FeatureFinder",StringList::create("mrm"))] = TDE("FeatureFinderMRM",StringList::create(""));

      map_[TDE("FeatureLinker",StringList::create("labeled"))] = TDE("FeatureLinkerLabeled",StringList::create(""));
      map_[TDE("FeatureLinker",StringList::create("unlabeled"))] = TDE("FeatureLinkerUnlabeled",StringList::create(""));
      map_[TDE("FeatureLinker",StringList::create("unlabeled_qt"))] = TDE("FeatureLinkerUnlabeledQT",StringList::create(""));

      map_[TDE("NoiseFilter",StringList::create("gaussian"))] = TDE("NoiseFilterGaussian",StringList::create(""));
      map_[TDE("NoiseFilter",StringList::create("sgolay"))] = TDE("NoiseFilterSGolay",StringList::create(""));

      map_[TDE("MapAligner",StringList::create("apply_given_trafo"))] = TDE("MapRTTransformer",StringList::create(""));
      map_[TDE("MapAligner",StringList::create("identification"))] = TDE("MapAlignerIdentification",StringList::create(""));
      map_[TDE("MapAligner",StringList::create("pose_clustering"))] = TDE("MapAlignerPoseClustering",StringList::create(""));
      map_[TDE("MapAligner",StringList::create("spectrum_alignment"))] = TDE("MapAlignerSpectrum",StringList::create(""));

      map_[TDE("CompNovo",StringList::create("CompNovo"))] = TDE("CompNovo",StringList::create(""));
      map_[TDE("CompNovo",StringList::create("CompNovoCID"))] = TDE("CompNovoCID",StringList::create(""));

      // SpectraFilter ...
      // PILISModel ...

      map_[TDE("PeakPicker",StringList::create("wavelet"))] = TDE("PeakPickerWavelet",StringList::create(""));
      map_[TDE("PeakPicker",StringList::create("high_res"))] = TDE("PeakPickerHiRes",StringList::create(""));

      // ITRAQAnalyzer && MSSimulator: no need to list here, as the type simply was made optional (no additional tools)

    }

    return map_;
  }

  bool INIUpdater::getNewToolName(const String& old_name, const String& tools_type, String& new_name)
  {
    new_name = "";
    // try with type (as some new tools for one type might have the exact same name as old ones with several types
    //                e.g., CompNovo
    TDE old_withtype(old_name, StringList::create(tools_type));
    if (map_.has(old_withtype))
    {
      new_name = map_[old_withtype].name;
      return true;
    }

    // try without type
    TDE old_notype(old_name, StringList());
    if (map_.has(old_notype))
    {
      new_name = map_[old_notype].name;
      return true;
    }

    // default to ToolHandler
		if (ToolHandler::getTOPPToolList(true).has(old_name) || ToolHandler::getUtilList().has(old_name))
    {
      new_name = old_name;
      return true;
    }

    return false;
  }

  ToolMapping INIUpdater::map_;

  

} // namespace
