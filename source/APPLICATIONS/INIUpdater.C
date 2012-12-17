// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

  StringList INIUpdater::getToolNamesFromINI(const Param & ini) const
  {
    StringList tool_names;

    for (Param::ParamIterator it = ini.begin(); it != ini.end(); ++it)
    {
      if (it.getName().toQString().count(':') == 1 && it.getName().hasSuffix(":version"))
      {
        tool_names.push_back(it.getName().prefix(':'));
      }
    }
    return tool_names;
  }

  const ToolMapping & INIUpdater::getNameMapping()
  {
    if (map_.empty())
    {
      map_[TDE("FeatureFinder", StringList::create("centroided"))] = TDE("FeatureFinderCentroided", StringList::create(""));
      map_[TDE("FeatureFinder", StringList::create("isotope_wavelet"))] = TDE("FeatureFinderIsotopeWavelet", StringList::create(""));
      map_[TDE("FeatureFinder", StringList::create("mrm"))] = TDE("FeatureFinderMRM", StringList::create(""));

      map_[TDE("FeatureLinker", StringList::create("labeled"))] = TDE("FeatureLinkerLabeled", StringList::create(""));
      map_[TDE("FeatureLinker", StringList::create("unlabeled"))] = TDE("FeatureLinkerUnlabeled", StringList::create(""));
      map_[TDE("FeatureLinker", StringList::create("unlabeled_qt"))] = TDE("FeatureLinkerUnlabeledQT", StringList::create(""));

      map_[TDE("NoiseFilter", StringList::create("gaussian"))] = TDE("NoiseFilterGaussian", StringList::create(""));
      map_[TDE("NoiseFilter", StringList::create("sgolay"))] = TDE("NoiseFilterSGolay", StringList::create(""));

      map_[TDE("MapAligner", StringList::create("apply_given_trafo"))] = TDE("MapRTTransformer", StringList::create(""));
      map_[TDE("MapAligner", StringList::create("identification"))] = TDE("MapAlignerIdentification", StringList::create(""));
      map_[TDE("MapAligner", StringList::create("pose_clustering"))] = TDE("MapAlignerPoseClustering", StringList::create(""));
      map_[TDE("MapAligner", StringList::create("spectrum_alignment"))] = TDE("MapAlignerSpectrum", StringList::create(""));

      map_[TDE("CompNovo", StringList::create("CompNovo"))] = TDE("CompNovo", StringList::create(""));
      map_[TDE("CompNovo", StringList::create("CompNovoCID"))] = TDE("CompNovoCID", StringList::create(""));

      // SpectraFilter ...
      // PILISModel ...

      map_[TDE("PeakPicker", StringList::create("wavelet"))] = TDE("PeakPickerWavelet", StringList::create(""));
      map_[TDE("PeakPicker", StringList::create("high_res"))] = TDE("PeakPickerHiRes", StringList::create(""));

      // ITRAQAnalyzer && MSSimulator: no need to list here, as the type simply was made optional (no additional tools)

    }

    return map_;
  }

  bool INIUpdater::getNewToolName(const String & old_name, const String & tools_type, String & new_name)
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
