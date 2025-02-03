// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QStringList>
#include <QtCore/QDir>

namespace OpenMS
{
  ToolListType ToolHandler::getTOPPToolList(const bool includeGenericWrapper)
  {
    ToolListType tools_map;
    // Note: don't use special characters like slashes in category names (leads to subcategories in KNIME) 
    const auto cat_calibration = "Mass Correction and Calibration";
    const auto cat_centroiding = "Spectrum processing: Centroiding";
    const auto cat_crosslinking = "Cross-Linking";
    const auto cat_dev = "[for Developers]";
    const auto cat_file_converter = "File Converter";
    const auto cat_file_filter_extract_merge = "File Filtering, Extraction and Merging";
    const auto cat_ID_MTX = "Metabolite Identification";
    const auto cat_ID_proc = "Identification Processing";
    const auto cat_ID_search = "Identification of Proteins and Peptides (SearchEngines)";
    const auto cat_linking = "Feature Linking";
    const auto cat_map_align = "Map Alignment";
    const auto cat_misc = "Misc";
    const auto cat_QC = "Quality Control";
    const auto cat_quant = "Quantitation";
    const auto cat_rna = "RNA";
    const auto cat_signal_proc_misc = "Spectrum processing: Misc";
    const auto cat_signal_proc_smooth_normalize = "Spectrum Processing: Peak Smoothing and Normalization";
    const auto cat_targeted = "Targeted Experiments and OpenSWATH";
    const auto cat_topdown = "Top-Down";

    // STOP and read!
    // 1) add your tool in alphabetical order!
    // 2) if you add/change categories, also mirror your changes in doc/doxygen/public/TOPP.doxygen

    tools_map["AccurateMassSearch"] = Internal::ToolDescription("AccurateMassSearch", cat_ID_MTX);
    tools_map["AssayGeneratorMetabo"] = Internal::ToolDescription("AssayGeneratorMetabo", cat_targeted);
    tools_map["AssayGeneratorMetaboSirius"] = Internal::ToolDescription("AssayGeneratorMetaboSirius", cat_targeted);
    tools_map["BaselineFilter"] = Internal::ToolDescription("BaselineFilter", cat_signal_proc_smooth_normalize);
    tools_map["ClusterMassTraces"] = Internal::ToolDescription("ClusterMassTraces", cat_misc);
    tools_map["ClusterMassTracesByPrecursor"] = Internal::ToolDescription("ClusterMassTracesByPrecursor", cat_targeted);
    tools_map["CometAdapter"] = Internal::ToolDescription("CometAdapter", cat_ID_search);
    tools_map["ConsensusID"] = Internal::ToolDescription("ConsensusID", cat_ID_proc);
    tools_map["ConsensusMapNormalizer"] = Internal::ToolDescription("ConsensusMapNormalizer", cat_quant);
    tools_map["CVInspector"] = Internal::ToolDescription("CVInspector", cat_dev);
    tools_map["DatabaseFilter"] = Internal::ToolDescription("DatabaseFilter", cat_file_filter_extract_merge);
    tools_map["DatabaseSuitability"] = Internal::ToolDescription("DatabaseSuitability", cat_QC);
    tools_map["Decharger"] = Internal::ToolDescription("Decharger", cat_quant);
    tools_map["DecoyDatabase"] = Internal::ToolDescription("DecoyDatabase", cat_file_filter_extract_merge);
    tools_map["DeMeanderize"] = Internal::ToolDescription("DeMeanderize", cat_misc);
    tools_map["Digestor"] = Internal::ToolDescription("Digestor", cat_ID_proc);
    tools_map["DigestorMotif"] = Internal::ToolDescription("DigestorMotif", cat_ID_proc);
    tools_map["DTAExtractor"] = Internal::ToolDescription("DTAExtractor", cat_file_filter_extract_merge);
    tools_map["EICExtractor"] = Internal::ToolDescription("EICExtractor", cat_quant);
    tools_map["Epifany"] = Internal::ToolDescription("Epifany", cat_ID_proc);
    tools_map["ExecutePipeline"] = Internal::ToolDescription("ExecutePipeline", cat_misc);
    tools_map["ExternalCalibration"] = Internal::ToolDescription("ExternalCalibration", cat_calibration);
    tools_map["FalseDiscoveryRate"] = Internal::ToolDescription("FalseDiscoveryRate", cat_ID_proc);
    tools_map["FeatureFinderCentroided"] = Internal::ToolDescription("FeatureFinderCentroided", cat_quant);
    tools_map["FeatureFinderIdentification"] = Internal::ToolDescription("FeatureFinderIdentification", cat_quant);
    tools_map["FeatureFinderMetabo"] = Internal::ToolDescription("FeatureFinderMetabo", cat_quant);
    tools_map["FeatureFinderMetaboIdent"] = Internal::ToolDescription("FeatureFinderMetaboIdent", cat_quant);
    tools_map["FeatureFinderMultiplex"] = Internal::ToolDescription("FeatureFinderMultiplex", cat_quant);
    tools_map["FeatureLinkerLabeled"] = Internal::ToolDescription("FeatureLinkerLabeled", cat_linking);
    tools_map["FeatureLinkerUnlabeled"] = Internal::ToolDescription("FeatureLinkerUnlabeled", cat_linking);
    tools_map["FeatureLinkerUnlabeledKD"] = Internal::ToolDescription("FeatureLinkerUnlabeledKD", cat_linking);
    tools_map["FeatureLinkerUnlabeledQT"] = Internal::ToolDescription("FeatureLinkerUnlabeledQT", cat_linking);
    tools_map["FileConverter"] = Internal::ToolDescription("FileConverter", cat_file_converter);
    tools_map["FileFilter"] = Internal::ToolDescription("FileFilter", cat_file_filter_extract_merge);
    tools_map["FileInfo"] = Internal::ToolDescription("FileInfo", cat_file_filter_extract_merge);
    tools_map["FileMerger"] = Internal::ToolDescription("FileMerger", cat_file_filter_extract_merge);
    tools_map["FLASHDeconv"] = Internal::ToolDescription("FLASHDeconv", cat_topdown);
    tools_map["FuzzyDiff"] = Internal::ToolDescription("FuzzyDiff", cat_dev);
    // tools_map["GenericWrapper"] = ... (place any extra handling here)
    tools_map["GNPSExport"] = Internal::ToolDescription("GNPSExport", cat_file_converter);
    tools_map["HighResPrecursorMassCorrector"] = Internal::ToolDescription("HighResPrecursorMassCorrector", cat_calibration);
    tools_map["IDConflictResolver"] = Internal::ToolDescription("IDConflictResolver", cat_ID_proc);
    tools_map["IDDecoyProbability"] = Internal::ToolDescription("IDDecoyProbability", cat_ID_proc);
    tools_map["IDExtractor"] = Internal::ToolDescription("IDExtractor", cat_ID_proc);
    tools_map["IDFileConverter"] = Internal::ToolDescription("IDFileConverter", cat_file_converter);
    tools_map["IDFilter"] = Internal::ToolDescription("IDFilter", cat_file_filter_extract_merge);
    tools_map["IDMapper"] = Internal::ToolDescription("IDMapper", cat_ID_proc);
    tools_map["IDMassAccuracy"] = Internal::ToolDescription("IDMassAccuracy", cat_ID_proc);
    tools_map["IDMerger"] = Internal::ToolDescription("IDMerger", cat_file_filter_extract_merge);
    tools_map["IDPosteriorErrorProbability"] = Internal::ToolDescription("IDPosteriorErrorProbability", cat_ID_proc);
    tools_map["IDRipper"] = Internal::ToolDescription("IDRipper", cat_file_filter_extract_merge);
    tools_map["IDRTCalibration"] = Internal::ToolDescription("IDRTCalibration", cat_calibration);
    tools_map["IDScoreSwitcher"] = Internal::ToolDescription("IDScoreSwitcher", cat_ID_proc);
    tools_map["IDSplitter"] = Internal::ToolDescription("IDSplitter", cat_file_filter_extract_merge);
    tools_map["ImageCreator"] = Internal::ToolDescription("ImageCreator", cat_misc);
    tools_map["INIUpdater"] = Internal::ToolDescription("INIUpdater", cat_misc);
    tools_map["InternalCalibration"] = Internal::ToolDescription("InternalCalibration", cat_calibration);
    tools_map["IonMobilityBinning"] = Internal::ToolDescription("IonMobilityBinning", cat_file_filter_extract_merge);
    tools_map["IsobaricAnalyzer"] = Internal::ToolDescription("IsobaricAnalyzer", cat_quant);
    tools_map["JSONExporter"] = Internal::ToolDescription("JSONExporter", cat_dev);
    tools_map["LuciphorAdapter"] = Internal::ToolDescription("LuciphorAdapter", cat_ID_search);
    tools_map["MapAlignerIdentification"] = Internal::ToolDescription("MapAlignerIdentification", cat_map_align);
    tools_map["MapAlignerPoseClustering"] = Internal::ToolDescription("MapAlignerPoseClustering", cat_map_align);
    tools_map["MapAlignerTreeGuided"] = Internal::ToolDescription("MapAlignerTreeGuided", cat_map_align);
    tools_map["MapNormalizer"] = Internal::ToolDescription("MapNormalizer", cat_signal_proc_smooth_normalize);
    tools_map["MapRTTransformer"] = Internal::ToolDescription("MapRTTransformer", cat_map_align);
    tools_map["MapStatistics"] = Internal::ToolDescription("MapStatistics", cat_file_filter_extract_merge);
    tools_map["MaRaClusterAdapter"] = Internal::ToolDescription("MaRaClusterAdapter", cat_signal_proc_misc);
    tools_map["MascotAdapter"] = Internal::ToolDescription("MascotAdapter", cat_ID_search);
    tools_map["MascotAdapterOnline"] = Internal::ToolDescription("MascotAdapterOnline", cat_ID_search);
    tools_map["MassCalculator"] = Internal::ToolDescription("MassCalculator", cat_misc);
    tools_map["MassTraceExtractor"] = Internal::ToolDescription("MassTraceExtractor", cat_quant);
    tools_map["MetaboliteAdductDecharger"] = Internal::ToolDescription("MetaboliteAdductDecharger", cat_quant);
    tools_map["MetaboliteSpectralMatcher"] = Internal::ToolDescription("MetaboliteSpectralMatcher", cat_ID_MTX);
    tools_map["MetaProSIP"] = Internal::ToolDescription("MetaProSIP", cat_quant);
    tools_map["MRMMapper"] = Internal::ToolDescription("MRMMapper", cat_targeted);
    tools_map["MRMPairFinder"] = Internal::ToolDescription("MRMPairFinder", cat_targeted);
    tools_map["MRMTransitionGroupPicker"] = Internal::ToolDescription("MRMTransitionGroupPicker", cat_targeted);
    tools_map["MSFraggerAdapter"] = Internal::ToolDescription("MSFraggerAdapter", cat_ID_search);
    tools_map["MSGFPlusAdapter"] = Internal::ToolDescription("MSGFPlusAdapter", cat_ID_search);
    tools_map["MSstatsConverter"] = Internal::ToolDescription("MSstatsConverter", cat_file_converter);
    tools_map["MultiplexResolver"] = Internal::ToolDescription("MultiplexResolver", cat_quant);
    tools_map["MzMLSplitter"] = Internal::ToolDescription("MzMLSplitter", cat_file_filter_extract_merge);
    tools_map["MzTabExporter"] = Internal::ToolDescription("MzTabExporter", cat_file_converter);
    tools_map["NoiseFilterGaussian"] = Internal::ToolDescription("NoiseFilterGaussian", cat_signal_proc_smooth_normalize);
    tools_map["NoiseFilterSGolay"] = Internal::ToolDescription("NoiseFilterSGolay", cat_signal_proc_smooth_normalize);
    tools_map["NovorAdapter"] = Internal::ToolDescription("NovorAdapter", cat_ID_search);
    tools_map["NucleicAcidSearchEngine"] = Internal::ToolDescription("NucleicAcidSearchEngine", cat_rna);
    tools_map["OpenMSDatabasesInfo"] = Internal::ToolDescription("OpenMSDatabasesInfo", cat_dev);
    tools_map["OpenMSInfo"] = Internal::ToolDescription("OpenMSInfo", cat_misc);
    tools_map["OpenPepXL"] = Internal::ToolDescription("OpenPepXL", cat_crosslinking);
    tools_map["OpenPepXLLF"] = Internal::ToolDescription("OpenPepXLLF", cat_crosslinking);
    tools_map["OpenSwathAnalyzer"] = Internal::ToolDescription("OpenSwathAnalyzer", cat_targeted);
    tools_map["OpenSwathAssayGenerator"] = Internal::ToolDescription("OpenSwathAssayGenerator", cat_targeted);
    tools_map["OpenSwathChromatogramExtractor"] = Internal::ToolDescription("OpenSwathChromatogramExtractor", cat_targeted);
    tools_map["OpenSwathConfidenceScoring"] = Internal::ToolDescription("OpenSwathConfidenceScoring", cat_targeted);
    tools_map["OpenSwathDecoyGenerator"] = Internal::ToolDescription("OpenSwathDecoyGenerator", cat_targeted);
    tools_map["OpenSwathDIAPreScoring"] = Internal::ToolDescription("OpenSwathDIAPreScoring", cat_targeted);
    tools_map["OpenSwathFeatureXMLToTSV"] = Internal::ToolDescription("OpenSwathFeatureXMLToTSV", cat_targeted);
    tools_map["OpenSwathFileSplitter"] = Internal::ToolDescription("OpenSwathFileSplitter", cat_targeted);
    tools_map["OpenSwathMzMLFileCacher"] = Internal::ToolDescription("OpenSwathMzMLFileCacher", cat_targeted);
    tools_map["OpenSwathRewriteToFeatureXML"] = Internal::ToolDescription("OpenSwathRewriteToFeatureXML", cat_targeted);
    tools_map["OpenSwathRTNormalizer"] = Internal::ToolDescription("OpenSwathRTNormalizer", cat_targeted);
    tools_map["OpenSwathWorkflow"] = Internal::ToolDescription("OpenSwathWorkflow", cat_targeted);
    tools_map["PeakPickerHiRes"] = Internal::ToolDescription("PeakPickerHiRes", cat_centroiding);
    tools_map["PeakPickerIterative"] = Internal::ToolDescription("PeakPickerIterative", cat_centroiding);
    tools_map["PeptideIndexer"] = Internal::ToolDescription("PeptideIndexer", cat_ID_proc);
    tools_map["PercolatorAdapter"] = Internal::ToolDescription("PercolatorAdapter", cat_ID_proc);
    tools_map["PhosphoScoring"] = Internal::ToolDescription("PhosphoScoring", cat_ID_proc);
    tools_map["ProteinInference"] = Internal::ToolDescription("ProteinInference", cat_ID_proc);
    tools_map["ProteinQuantifier"] = Internal::ToolDescription("ProteinQuantifier", cat_quant);
    tools_map["ProteomicsLFQ"] = Internal::ToolDescription("ProteomicsLFQ", cat_quant);
    tools_map["PSMFeatureExtractor"] = Internal::ToolDescription("PSMFeatureExtractor", cat_ID_proc);
    tools_map["QCCalculator"] = Internal::ToolDescription("QCCalculator", cat_QC);
    tools_map["QCEmbedder"] = Internal::ToolDescription("QCEmbedder", cat_QC);
    tools_map["QCExporter"] = Internal::ToolDescription("QCExporter", cat_QC);
    tools_map["QCExtractor"] = Internal::ToolDescription("QCExtractor", cat_QC);
    tools_map["QCImporter"] = Internal::ToolDescription("QCImporter", cat_QC);
    tools_map["QCMerger"] = Internal::ToolDescription("QCMerger", cat_QC);
    tools_map["QCShrinker"] = Internal::ToolDescription("QCExporter", cat_QC);
    tools_map["QualityControl"] = Internal::ToolDescription("QualityControl", cat_QC);
    tools_map["Resampler"] = Internal::ToolDescription("Resampler", cat_signal_proc_misc);
    tools_map["RNADigestor"] = Internal::ToolDescription("RNADigestor", cat_rna);
    tools_map["RNAMassCalculator"] = Internal::ToolDescription("RNAMassCalculator", cat_rna);
    tools_map["RNPxlXICFilter"] = Internal::ToolDescription("RNPxlXICFilter", cat_crosslinking);
    tools_map["SageAdapter"] = Internal::ToolDescription("SageAdapter", cat_ID_search);
    tools_map["SeedListGenerator"] = Internal::ToolDescription("SeedListGenerator", cat_quant);
    tools_map["SemanticValidator"] = Internal::ToolDescription("SemanticValidator", cat_dev);
    tools_map["SequenceCoverageCalculator"] = Internal::ToolDescription("SequenceCoverageCalculator", cat_ID_proc);
    tools_map["SimpleSearchEngine"] = Internal::ToolDescription("SimpleSearchEngine", cat_ID_search);
    tools_map["SiriusExport"] = Internal::ToolDescription("SiriusExport", cat_ID_MTX);
    tools_map["SpecLibCreator"] = Internal::ToolDescription("SpecLibCreator", cat_ID_proc);
    tools_map["SpecLibSearcher"] = Internal::ToolDescription("SpecLibSearcher", cat_ID_search);
    tools_map["SpectraFilterNLargest"] = Internal::ToolDescription("SpectraFilterNLargest", cat_signal_proc_smooth_normalize);
    tools_map["SpectraFilterNormalizer"] = Internal::ToolDescription("SpectraFilterNormalizer", cat_signal_proc_smooth_normalize);
    tools_map["SpectraFilterThresholdMower"] = Internal::ToolDescription("SpectraFilterThresholdMower", cat_signal_proc_smooth_normalize);
    tools_map["SpectraFilterWindowMower"] = Internal::ToolDescription("SpectraFilterWindowMower", cat_signal_proc_smooth_normalize);
    tools_map["SpectraMerger"] = Internal::ToolDescription("SpectraMerger", cat_signal_proc_misc);
    tools_map["SpectraSTSearchAdapter"] = Internal::ToolDescription("SpectraSTSearchAdapter", cat_ID_search);
    tools_map["StaticModification"] = Internal::ToolDescription("StaticModification", cat_ID_proc);
    tools_map["TargetedFileConverter"] = Internal::ToolDescription("TargetedFileConverter", cat_file_converter);
    tools_map["TextExporter"] = Internal::ToolDescription("TextExporter", cat_file_converter);
    tools_map["TICCalculator"] = Internal::ToolDescription("TICCalculator", cat_misc);
    tools_map["TriqlerConverter"] = Internal::ToolDescription("TriqlerConverter", cat_file_converter);
    tools_map["XFDR"] = Internal::ToolDescription("XFDR", cat_crosslinking);
    tools_map["XMLValidator"] = Internal::ToolDescription("XMLValidator", cat_dev);
    tools_map["XTandemAdapter"] = Internal::ToolDescription("XTandemAdapter", cat_ID_search);
    
    // STOP! insert your tool in alphabetical order for easier maintenance (tools requiring the GUI lib should be added below **in addition**)

    // ATTENTION: tools requiring the GUI lib
#ifndef WITH_GUI
    StringList GUI_tools = {
      "ExecutePipeline",
      "ImageCreator",
      "INIUpdater",
      "Resampler",
    };
    std::for_each(GUI_tools.begin(), GUI_tools.end(), [&tools_map](const String& del) {
      tools_map.erase(del);
    });
#endif

    // INTERNAL tools
    // this operation is expensive, as we need to parse configuration files (*.ttd)
    std::vector<Internal::ToolDescription> internal_tools = getInternalTools_();
    for (std::vector<Internal::ToolDescription>::const_iterator it = internal_tools.begin(); it != internal_tools.end(); ++it)
    {
      if (tools_map.find(it->name) == tools_map.end())
      {
        tools_map[it->name] = *it;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Duplicate tool name error: Trying to add internal tool '" + it->name, it->name);
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

  StringList ToolHandler::getTypes(const String& toolname)
  {
    Internal::ToolDescription ret;
    ToolListType tools;
    if (toolname == "GenericWrapper")
    {
      tools = getTOPPToolList(true);
    }
    else
    {
      tools = getTOPPToolList();
    }
    if (tools.find(toolname) != tools.end())
    {
      return tools[toolname].types;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Requested tool '" + toolname + "' does not exist!", toolname);
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
    for (int i = 0; i < files.size(); ++i)
    {
      ToolDescriptionFile tdf;
      std::vector<Internal::ToolDescription> tools;
      tdf.load(String(files[i]), tools);
      // add every tool from file to list
      for (Size i_t = 0; i_t < tools.size(); ++i_t)
      {
        if (i == 0 && i_t == 0)
          {
            tools_external_ = tools[i_t]; // init
          }
        else
          {
            tools_external_.append(tools[i_t]); // append
          }
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
    paths << (getExternalToolsPath() + "/WINDOWS").toQString();
#else
    paths << (getExternalToolsPath() + "/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_PATH") != nullptr)
    {
      paths << String(getenv("OPENMS_TTD_PATH")).toQString();
    }

    QStringList all_files;
    for (int p = 0; p < paths.size(); ++p)
    {
      QDir dir(paths[p], "*.ttd");
      QStringList files = dir.entryList();
      for (int i = 0; i < files.size(); ++i)
      {
        files[i] = dir.absolutePath() + QDir::separator() + files[i];
      }
      all_files << files;
    }
    //StringList list = ListUtils::create<String>(getExternalToolsPath() + "/" + "msconvert.ttd");
    return all_files;
  }

  QStringList ToolHandler::getInternalToolConfigFiles_()
  {
    QStringList paths;
    // *.ttd default path
    paths << getInternalToolsPath().toQString();
    // OS-specific path
#ifdef OPENMS_WINDOWSPLATFORM
    paths << (getInternalToolsPath() + "/WINDOWS").toQString();
#else
    paths << (getInternalToolsPath() + "/LINUX").toQString();
#endif
    // additional environment
    if (getenv("OPENMS_TTD_INTERNAL_PATH") != nullptr)
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
        files[i] = dir.absolutePath() + QDir::separator() + files[i];
      }
      all_files << files;
    }
    return all_files;
  }

  String ToolHandler::getCategory(const String& toolname)
  {
    ToolListType tools = getTOPPToolList(true);
    String s;
    if (tools.find(toolname) != tools.end())
    {
      s = tools[toolname].category;
    }

    return s;
  }

  // static
  Internal::ToolDescription ToolHandler::tools_external_ = Internal::ToolDescription();
  std::vector<Internal::ToolDescription> ToolHandler::tools_internal_;
  bool ToolHandler::tools_external_loaded_ = false;
  bool ToolHandler::tools_internal_loaded_ = false;

} // namespace
