// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <QStringList>
#include <QDir>


#include <vector>

namespace OpenMS
{
	ToolListType ToolHandler::getTOPPToolList(const bool includeGenericWrapper)
	{
		ToolListType tools_map;

    tools_map["AdditiveSeries"] = Internal::ToolDescription("AdditiveSeries", "Quantitation");
		tools_map["BaselineFilter"] = Internal::ToolDescription("BaselineFilter", "Signal processing and preprocessing");
    tools_map["CompNovo"] = Internal::ToolDescription("CompNovo", "Identification");
    tools_map["CompNovoCID"] = Internal::ToolDescription("CompNovoCID", "Identification");
    tools_map["ConsensusID"] = Internal::ToolDescription("ConsensusID", "ID Processing");
		tools_map["ConsensusMapNormalizer"] = Internal::ToolDescription("ConsensusMapNormalizer", "Map Alignment");
		tools_map["DBExporter"] = Internal::ToolDescription("DBExporter", "File Handling");
		tools_map["DBImporter"] = Internal::ToolDescription("DBImporter", "File Handling");
		tools_map["DTAExtractor"] = Internal::ToolDescription("DTAExtractor", "File Handling");
		tools_map["Decharger"] = Internal::ToolDescription("Decharger", "Quantitation");
		tools_map["ExecutePipeline"] = Internal::ToolDescription("ExecutePipeline", "Misc");
    tools_map["EICExtractor"] = Internal::ToolDescription("EICExtractor", "Quantitation");
    tools_map["FalseDiscoveryRate"] = Internal::ToolDescription("FalseDiscoveryRate", "ID Processing");
		tools_map["FeatureFinderCentroided"] = Internal::ToolDescription("FeatureFinderCentroided", "Quantitation");
		tools_map["FeatureFinderIsotopeWavelet"] = Internal::ToolDescription("FeatureFinderIsotopeWavelet", "Quantitation");
		tools_map["FeatureFinderMRM"] = Internal::ToolDescription("FeatureFinderMRM", "Quantitation");
		tools_map["FeatureLinkerLabeled"] = Internal::ToolDescription("FeatureLinkerLabeled", "Map Alignment");
		tools_map["FeatureLinkerUnlabeled"] = Internal::ToolDescription("FeatureLinkerUnlabeled", "Map Alignment");
		tools_map["FeatureLinkerUnlabeledQT"] = Internal::ToolDescription("FeatureLinkerUnlabeledQT", "Map Alignment");
		tools_map["FileConverter"] = Internal::ToolDescription("FileConverter", "File Handling");
		tools_map["FileFilter"] = Internal::ToolDescription("FileFilter", "File Handling");
		tools_map["FileInfo"] = Internal::ToolDescription("FileInfo", "File Handling");
		tools_map["FileMerger"] = Internal::ToolDescription("FileMerger", "File Handling");
    tools_map["HighResPrecursorMassCorrector"] = Internal::ToolDescription("HighResPrecursorMassCorrector", "Signal processing and preprocessing");
		tools_map["IDConflictResolver"] = Internal::ToolDescription("IDConflictResolver", "ID Processing");
    tools_map["IDDecoyProbability"] = Internal::ToolDescription("IDDecoyProbability", "ID Processing");
    tools_map["IDFileConverter"] = Internal::ToolDescription("IDFileConverter", "ID Processing");
    tools_map["IDFilter"] = Internal::ToolDescription("IDFilter", "ID Processing");
    tools_map["IDMapper"] = Internal::ToolDescription("IDMapper", "ID Processing");
		tools_map["IDMerger"] = Internal::ToolDescription("IDMerger", "File Handling");
    tools_map["IDPosteriorErrorProbability"] = Internal::ToolDescription("IDPosteriorErrorProbability", "ID Processing");
    tools_map["IDRTCalibration"] = Internal::ToolDescription("IDRTCalibration", "ID Processing");
		tools_map["ITRAQAnalyzer"] = Internal::ToolDescription("ITRAQAnalyzer", "Quantitation", StringList::create("4plex,8plex"));
		tools_map["InclusionExclusionListCreator"] = Internal::ToolDescription("InclusionExclusionListCreator", "Targeted Experiments");
    tools_map["InspectAdapter"] = Internal::ToolDescription("InspectAdapter", "Identification");
		tools_map["InternalCalibration"] = Internal::ToolDescription("InternalCalibration", "Signal processing and preprocessing");
		tools_map["MapAligner"] = Internal::ToolDescription("MapAligner", "Map Alignment", Factory<MapAlignmentAlgorithm>::registeredProducts());
		tools_map["MapAlignerApplyTransformation"] = Internal::ToolDescription("MapAlignerApplyTransformation", "Map Alignment");
		tools_map["MapAlignerIdentification"] = Internal::ToolDescription("MapAlignerIdentification", "Map Alignment");
		tools_map["MapAlignerPoseClustering"] = Internal::ToolDescription("MapAlignerPoseClustering", "Map Alignment");
		tools_map["MapAlignerSpectrum"] = Internal::ToolDescription("MapAlignerSpectrum", "Map Alignment");
		tools_map["MapNormalizer"] = Internal::ToolDescription("MapNormalizer", "Signal processing and preprocessing");
    tools_map["MapStatistics"] = Internal::ToolDescription("MapStatistics", "File Handling");
    tools_map["MascotAdapter"] = Internal::ToolDescription("MascotAdapter", "Identification");
    tools_map["MascotAdapterOnline"] = Internal::ToolDescription("MascotAdapterOnline", "Identification");
    tools_map["MassTraceExtractor"] = Internal::ToolDescription("MassTraceExtractor", "Detects mass traces in LC-MS data");
		tools_map["NoiseFilterGaussian"] = Internal::ToolDescription("NoiseFilterGaussian", "Signal processing and preprocessing");
		tools_map["NoiseFilterSGolay"] = Internal::ToolDescription("NoiseFilterSGolay", "Signal processing and preprocessing");
    tools_map["OMSSAAdapter"] = Internal::ToolDescription("OMSSAAdapter", "Identification");
    tools_map["PILISIdentification"] = Internal::ToolDescription("PILISIdentification", "Identification");
    tools_map["PILISModel"] = Internal::ToolDescription("PILISModel", "ID Processing");
		tools_map["PTModel"] = Internal::ToolDescription("PTModel", "Peptide property prediction");
		tools_map["PTPredict"] = Internal::ToolDescription("PTPredict", "Peptide property prediction");
		tools_map["PeakPickerHiRes"] = Internal::ToolDescription("PeakPickerHiRes", "Signal processing and preprocessing");
		tools_map["PeakPickerWavelet"] = Internal::ToolDescription("PeakPickerWavelet", "Signal processing and preprocessing");
    tools_map["PepNovoAdapter"] = Internal::ToolDescription("PepNovoAdapter", "Identification");
    tools_map["PeptideIndexer"] = Internal::ToolDescription("PeptideIndexer", "ID Processing");
		tools_map["PrecursorIonSelector"] = Internal::ToolDescription("PrecursorIonSelector", "Targeted Experiments");
    tools_map["PrecursorMassCorrector"] = Internal::ToolDescription("PrecursorMassCorrector", "Signal processing and preprocessing");
    tools_map["ProteinInference"] = Internal::ToolDescription("ProteinInference", "Identification");
		tools_map["ProteinQuantifier"] = Internal::ToolDescription("ProteinQuantifier", "Quantitation");
		tools_map["RTModel"] = Internal::ToolDescription("RTModel", "Peptide property prediction");
		tools_map["RTPredict"] = Internal::ToolDescription("RTPredict", "Peptide property prediction");
		tools_map["Resampler"] = Internal::ToolDescription("Resampler", "Signal processing and preprocessing");
    tools_map["SILACAnalyzer"] = Internal::ToolDescription("SILACAnalyzer", "Quantitation");
		tools_map["SeedListGenerator"] = Internal::ToolDescription("SeedListGenerator", "Quantitation");
    //tools_map["SequestAdapter"] = Internal::ToolDescription("SequestAdapter", "Identification");
    tools_map["SpecLibSearcher"] = Internal::ToolDescription("SpecLibSearcher", "Identification");
    tools_map["SpectraFilter"] = Internal::ToolDescription("SpectraFilter", "Identification", Factory<PreprocessingFunctor>::registeredProducts());
		tools_map["SpectraMerger"] = Internal::ToolDescription("SpectraMerger", "File Handling");
    tools_map["TOFCalibration"] = Internal::ToolDescription("TOFCalibration", "Signal processing and preprocessing");
		tools_map["TextExporter"] = Internal::ToolDescription("TextExporter", "File Handling");
    tools_map["XTandemAdapter"] = Internal::ToolDescription("XTandemAdapter", "Identification");

    // INTERNAL tools
    // this operation is expensive, as we need to parse configuration files (*.ttd)
    std::vector<Internal::ToolDescription> internal_tools = getInternalTools_();
    for (std::vector<Internal::ToolDescription>::const_iterator it = internal_tools.begin(); it != internal_tools.end(); ++it)
    {
      if (tools_map.find(it->name) == tools_map.end())
      {
        tools_map[it->name] = *it;
      } else
      {
        throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Duplicate tool name error: Trying to add internal tool '" + it->name, it->name);
      }
    }

    // EXTERNAL tools
    // this operation is expensive, as we need to parse configuration files (*.ttd)
    if (includeGenericWrapper)
    {
      tools_map["GenericWrapper"] = getExternalTools_();
    }

		return tools_map;
	}

	ToolListType ToolHandler::getUtilList()
	{
		ToolListType util_map;

		util_map["IDMassAccuracy"] = Internal::ToolDescription("IDMassAccuracy", "");
		util_map["DecoyDatabase"] = Internal::ToolDescription("DecoyDatabase", "");
		util_map["MapAlignmentEvaluation"] = Internal::ToolDescription("MapAlignmentEvaluation", "");
		util_map["CaapConvert"] = Internal::ToolDescription("CaapConvert", "");
		util_map["CVInspector"] = Internal::ToolDescription("CVInspector", "");
		util_map["DecoyDatabase"] = Internal::ToolDescription("DecoyDatabase", "");
		util_map["Digestor"] = Internal::ToolDescription("Digestor", "");
		util_map["DigestorMotif"] = Internal::ToolDescription("DigestorMotif", "");
		util_map["FFEval"] = Internal::ToolDescription("FFEval", "");
		util_map["FuzzyDiff"] = Internal::ToolDescription("FuzzyDiff", "");
		util_map["HistView"] = Internal::ToolDescription("HistView", "");
		util_map["IDExtractor"] = Internal::ToolDescription("IDExtractor", "");
		util_map["LabeledEval"] = Internal::ToolDescription("LabeledEval", "");
		util_map["SemanticValidator"] = Internal::ToolDescription("SemanticValidator", "");
		util_map["SequenceCoverageCalculator"] = Internal::ToolDescription("SequenceCoverageCalculator", "");
		util_map["XMLValidator"] = Internal::ToolDescription("XMLValidator", "");
		util_map["IdXMLEvaluation"] = Internal::ToolDescription("IdXMLEvaluation", "");
    util_map["MSSimulator"] = Internal::ToolDescription("MSSimulator", "");
		util_map["ERPairFinder"] = Internal::ToolDescription("ERPairFinder", "");
		util_map["SpecLibCreator"] = Internal::ToolDescription("SpecLibCreator", "");
		util_map["SvmTheoreticalSpectrumGeneratorTrainer"] = Internal::ToolDescription("SvmTheoreticalSpectrumGeneratorTrainer", "");
		util_map["MRMPairFinder"] = Internal::ToolDescription("MRMPairFinder", "");
		util_map["DeMeanderize"] = Internal::ToolDescription("DeMeanderize", "");
		util_map["UniqueIdAssigner"] = Internal::ToolDescription("UniqueIdAssigner", "");
		util_map["ImageCreator"] = Internal::ToolDescription("ImageCreator", "");
		util_map["IDSplitter"] = Internal::ToolDescription("IDSplitter", "");
		util_map["MassCalculator"] = Internal::ToolDescription("MassCalculator", "");
		util_map["TransformationEvaluation"] = Internal::ToolDescription("TransformationEvaluation", "");

		return util_map;
	}

  StringList ToolHandler::getTypes(const String& toolname)
  {
    // for internal tools, query TOPP and UTILS for a match
    Internal::ToolDescription ret;
    if (getUtilList().has(toolname))
    {
      return getUtilList()[toolname].types;
    }
    else
    {
      ToolListType tools;
      if (toolname == "GenericWrapper") tools = getTOPPToolList(true);
      else tools = getTOPPToolList();

      if (tools.has(toolname))
      {
        return tools[toolname].types;
      }
    }
    throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Requested toolname '" + toolname + " does not exist!", toolname);
    return StringList();
  }


  std::vector<Internal::ToolDescription> ToolHandler::getInternalTools_()
  {
    if (!tools_internal_loaded_)
    {
      loadInternalToolConfig_();
      tools_internal_loaded_ = true;
    }
    return tools_internal_;
  }  

  String ToolHandler::getExternalToolsPath()
  {
    return File::getOpenMSDataPath() + "/TOOLS/EXTERNAL";
  }

  String ToolHandler::getInternalToolsPath()
  {
    return File::getOpenMSDataPath() + "/TOOLS/INTERNAL";
  }

  Internal::ToolDescription ToolHandler::getExternalTools_()
  {
    if (!tools_external_loaded_) 
    {
      loadExternalToolConfig_();
      tools_external_loaded_ = true;
    }

    return tools_external_;
  }

  void ToolHandler::loadExternalToolConfig_()
  {
    QStringList files = getExternalToolConfigFiles_();
    for (int i=0;i<files.size();++i)
    {
      ToolDescriptionFile tdf;
      std::vector<Internal::ToolDescription> tools;
      tdf.load(String(files[i]), tools);
      // add every tool from file to list
      for (Size i_t=0; i_t<tools.size(); ++i_t)
      {
        if (i==0 && i_t==0) tools_external_ = tools[i_t]; // init
        else tools_external_.append(tools[i_t]); // append
      }
    }
    tools_external_.name = "GenericWrapper";
    tools_external_.category = "EXTERNAL";
  }

  void ToolHandler::loadInternalToolConfig_()
  {
    QStringList files = getInternalToolConfigFiles_();
    for (int i = 0; i < files.size(); ++i)
    {
      ToolDescriptionFile tdf;
      std::vector<Internal::ToolDescription> tools;
      tdf.load(String(files[i]), tools);
      // add every tool from file to list
      for (Size i_t = 0; i_t < tools.size(); ++i_t)
      {
          tools_internal_.push_back(tools[i_t]);
          tools_external_.category = "INTERNAL";
      }
    }
  }

  QStringList ToolHandler::getExternalToolConfigFiles_()
  {
    
    QStringList paths;
    // *.ttd default path
    paths << getExternalToolsPath().toQString(); 
    // OS-specific path
#ifdef OPENMS_WINDOWSPLATFORM
    paths << (getExternalToolsPath()+"/WINDOWS").toQString();
#else
    paths << (getExternalToolsPath()+"/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_PATH")!=0)
    {
      paths << String(getenv("OPENMS_TTD_PATH")).toQString();
    }
    
    QStringList all_files;
    for (int p=0;p<paths.size();++p)
    {
      QDir dir(paths[p], "*.ttd");
      QStringList files = dir.entryList();
      for (int i=0;i<files.size();++i)
      {
        files[i] = dir.absolutePath()+QDir::separator()+files[i];
      }
      all_files << files;
    }
    //StringList list = StringList::create(getExternalToolsPath() + "/" + "msconvert.ttd");
    return all_files;
  }

  QStringList ToolHandler::getInternalToolConfigFiles_()
  {
    QStringList paths;
    // *.ttd default path
    paths << getInternalToolsPath().toQString();
    // OS-specific path
#ifdef OPENMS_WINDOWSPLATFORM
    paths << (getInternalToolsPath()+"/WINDOWS").toQString();
#else
    paths << (getInternalToolsPath()+"/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_INTERNAL_PATH")!=0)
    {
      paths << String(getenv("OPENMS_TTD_INTERNAL_PATH")).toQString();
    }

    QStringList all_files;
    for (int p = 0; p < paths.size(); ++p)
    {
      QDir dir(paths[p], "*.ttd");
      QStringList files = dir.entryList();
      for (int i = 0; i < files.size(); ++i)
      {
        files[i] = dir.absolutePath()+QDir::separator()+files[i];
      }
      all_files << files;
    }
    return all_files;
  }

  String ToolHandler::getCategory(const String& toolname)
  {
    ToolListType tools = getTOPPToolList(true);
    ToolListType utils = getUtilList();
    String s;
    if (tools.has(toolname))
    {
      s = tools[toolname].category;
    }
    else if (utils.has(toolname))
    {
      s = utils[toolname].category;
    }

    return s;
  }

  // static
  Internal::ToolDescription ToolHandler::tools_external_ = Internal::ToolDescription();  
  std::vector<Internal::ToolDescription> ToolHandler::tools_internal_;
  bool ToolHandler::tools_external_loaded_ = false;
  bool ToolHandler::tools_internal_loaded_ = false;

} // namespace
