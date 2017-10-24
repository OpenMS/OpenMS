// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <OpenMS/SYSTEM/File.h>

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
    tools_map["CometAdapter"] = Internal::ToolDescription("CometAdapter", "Identification");
    tools_map["Decharger"] = Internal::ToolDescription("Decharger", "Quantitation");
    tools_map["DTAExtractor"] = Internal::ToolDescription("DTAExtractor", "File Handling");
    tools_map["EICExtractor"] = Internal::ToolDescription("EICExtractor", "Quantitation");
    tools_map["ExternalCalibration"] = Internal::ToolDescription("ExternalCalibration", "Signal processing and preprocessing");
    tools_map["FalseDiscoveryRate"] = Internal::ToolDescription("FalseDiscoveryRate", "ID Processing");
    tools_map["FeatureFinderCentroided"] = Internal::ToolDescription("FeatureFinderCentroided", "Quantitation");
    tools_map["FeatureFinderIdentification"] = Internal::ToolDescription("FeatureFinderIdentification", "Quantitation");
    tools_map["FeatureFinderIsotopeWavelet"] = Internal::ToolDescription("FeatureFinderIsotopeWavelet", "Quantitation");
    tools_map["FeatureFinderMetabo"] = Internal::ToolDescription("FeatureFinderMetabo", "Quantitation");
    tools_map["FeatureFinderMultiplex"] = Internal::ToolDescription("FeatureFinderMultiplex", "Quantitation");
    tools_map["FeatureFinderMRM"] = Internal::ToolDescription("FeatureFinderMRM", "Quantitation");
    tools_map["FeatureFinderSuperHirn"] = Internal::ToolDescription("FeatureFinderSuperHirn", "Quantitation");
    tools_map["FeatureLinkerLabeled"] = Internal::ToolDescription("FeatureLinkerLabeled", "Map Alignment");
    tools_map["FeatureLinkerUnlabeled"] = Internal::ToolDescription("FeatureLinkerUnlabeled", "Map Alignment");
    tools_map["FeatureLinkerUnlabeledKD"] = Internal::ToolDescription("FeatureLinkerUnlabeledKD", "Map Alignment");
    tools_map["FeatureLinkerUnlabeledQT"] = Internal::ToolDescription("FeatureLinkerUnlabeledQT", "Map Alignment");
    tools_map["FidoAdapter"] = Internal::ToolDescription("FidoAdapter", "ID Processing");
    tools_map["FileConverter"] = Internal::ToolDescription("FileConverter", "File Handling");
    tools_map["FileFilter"] = Internal::ToolDescription("FileFilter", "File Handling");
    tools_map["FileInfo"] = Internal::ToolDescription("FileInfo", "File Handling");
    tools_map["FileMerger"] = Internal::ToolDescription("FileMerger", "File Handling");
    tools_map["HighResPrecursorMassCorrector"] = Internal::ToolDescription("HighResPrecursorMassCorrector", "Signal processing and preprocessing");
    tools_map["IDConflictResolver"] = Internal::ToolDescription("IDConflictResolver", "ID Processing");
    tools_map["IDFileConverter"] = Internal::ToolDescription("IDFileConverter", "ID Processing");
    tools_map["IDFilter"] = Internal::ToolDescription("IDFilter", "ID Processing");
    tools_map["IDMapper"] = Internal::ToolDescription("IDMapper", "ID Processing");
    tools_map["IDMerger"] = Internal::ToolDescription("IDMerger", "File Handling");
    tools_map["IDPosteriorErrorProbability"] = Internal::ToolDescription("IDPosteriorErrorProbability", "ID Processing");
    tools_map["IDRipper"] = Internal::ToolDescription("IDRipper", "File Handling");
    tools_map["IDRTCalibration"] = Internal::ToolDescription("IDRTCalibration", "ID Processing");
    tools_map["IsobaricAnalyzer"] = Internal::ToolDescription("IsobaricAnalyzer", "Quantitation");
    tools_map["InclusionExclusionListCreator"] = Internal::ToolDescription("InclusionExclusionListCreator", "Targeted Experiments");
    tools_map["InspectAdapter"] = Internal::ToolDescription("InspectAdapter", "Identification");
    tools_map["InternalCalibration"] = Internal::ToolDescription("InternalCalibration", "Signal processing and preprocessing");
    tools_map["LuciphorAdapter"] = Internal::ToolDescription("LuciphorAdapter", "ID Processing");
    tools_map["MapAlignerIdentification"] = Internal::ToolDescription("MapAlignerIdentification", "Map Alignment");
    tools_map["MapAlignerPoseClustering"] = Internal::ToolDescription("MapAlignerPoseClustering", "Map Alignment");
    tools_map["MapAlignerSpectrum"] = Internal::ToolDescription("MapAlignerSpectrum", "Map Alignment");
    tools_map["MapNormalizer"] = Internal::ToolDescription("MapNormalizer", "Signal processing and preprocessing");
    tools_map["MapRTTransformer"] = Internal::ToolDescription("MapRTTransformer", "Map Alignment");
    tools_map["MapStatistics"] = Internal::ToolDescription("MapStatistics", "File Handling");
    tools_map["MascotAdapter"] = Internal::ToolDescription("MascotAdapter", "Identification");
    tools_map["MascotAdapterOnline"] = Internal::ToolDescription("MascotAdapterOnline", "Identification");
    tools_map["MassTraceExtractor"] = Internal::ToolDescription("MassTraceExtractor", "Signal processing and preprocessing");
    tools_map["MRMMapper"] = Internal::ToolDescription("MRMMapper", "Targeted Experiments");
    tools_map["MSGFPlusAdapter"] = Internal::ToolDescription("MSGFPlusAdapter", "Identification");
    tools_map["MyriMatchAdapter"] = Internal::ToolDescription("MyriMatchAdapter", "Identification");
    tools_map["MzTabExporter"] = Internal::ToolDescription("MzTabExporter", "File Handling");
    tools_map["NoiseFilterGaussian"] = Internal::ToolDescription("NoiseFilterGaussian", "Signal processing and preprocessing");
    tools_map["NoiseFilterSGolay"] = Internal::ToolDescription("NoiseFilterSGolay", "Signal processing and preprocessing");
    tools_map["OMSSAAdapter"] = Internal::ToolDescription("OMSSAAdapter", "Identification");
    tools_map["OpenSwathAnalyzer"] = Internal::ToolDescription("OpenSwathAnalyzer", "Targeted Experiments");
    tools_map["OpenSwathAssayGenerator"] = Internal::ToolDescription("OpenSwathAssayGenerator", "Targeted Experiments");
    tools_map["OpenSwathChromatogramExtractor"] = Internal::ToolDescription("OpenSwathChromatogramExtractor", "Targeted Experiments");
    tools_map["OpenSwathConfidenceScoring"] = Internal::ToolDescription("OpenSwathConfidenceScoring", "Targeted Experiments");
    tools_map["OpenSwathDecoyGenerator"] = Internal::ToolDescription("OpenSwathDecoyGenerator", "Targeted Experiments");
    tools_map["OpenSwathFeatureXMLToTSV"] = Internal::ToolDescription("OpenSwathFeatureXMLToTSV", "Targeted Experiments");
    tools_map["OpenSwathRTNormalizer"] = Internal::ToolDescription("OpenSwathRTNormalizer", "Targeted Experiments");
    tools_map["PeakPickerHiRes"] = Internal::ToolDescription("PeakPickerHiRes", "Signal processing and preprocessing");
    tools_map["PeakPickerWavelet"] = Internal::ToolDescription("PeakPickerWavelet", "Signal processing and preprocessing");
    tools_map["PepNovoAdapter"] = Internal::ToolDescription("PepNovoAdapter", "Identification");
    tools_map["PeptideIndexer"] = Internal::ToolDescription("PeptideIndexer", "ID Processing");
    tools_map["PercolatorAdapter"] = Internal::ToolDescription("PercolatorAdapter", "ID Processing");
    tools_map["PhosphoScoring"] = Internal::ToolDescription("PhosphoScoring", "ID Processing");
    tools_map["PrecursorIonSelector"] = Internal::ToolDescription("PrecursorIonSelector", "Targeted Experiments");
    tools_map["PrecursorMassCorrector"] = Internal::ToolDescription("PrecursorMassCorrector", "Signal processing and preprocessing");
    tools_map["ProteinInference"] = Internal::ToolDescription("ProteinInference", "Identification");
    tools_map["ProteinQuantifier"] = Internal::ToolDescription("ProteinQuantifier", "Quantitation");
    tools_map["ProteinResolver"] = Internal::ToolDescription("ProteinResolver", "Quantitation");
    tools_map["PTModel"] = Internal::ToolDescription("PTModel", "Peptide property prediction");
    tools_map["PTPredict"] = Internal::ToolDescription("PTPredict", "Peptide property prediction");
    tools_map["RTModel"] = Internal::ToolDescription("RTModel", "Peptide property prediction");
    tools_map["RTPredict"] = Internal::ToolDescription("RTPredict", "Peptide property prediction");
    tools_map["SeedListGenerator"] = Internal::ToolDescription("SeedListGenerator", "Quantitation");
    //tools_map["SequestAdapter"] = Internal::ToolDescription("SequestAdapter", "Identification");
    tools_map["SpecLibSearcher"] = Internal::ToolDescription("SpecLibSearcher", "Identification");
    tools_map["SpectraFilterBernNorm"] = Internal::ToolDescription("SpectraFilterBernNorm", "Identification");
    tools_map["SpectraFilterMarkerMower"] = Internal::ToolDescription("SpectraFilterMarkerMower", "Identification");
    tools_map["SpectraFilterNLargest"] = Internal::ToolDescription("SpectraFilterNLargest", "Identification");
    tools_map["SpectraFilterNormalizer"] = Internal::ToolDescription("SpectraFilterNormalizer", "Identification");
    tools_map["SpectraFilterParentPeakMower"] = Internal::ToolDescription("SpectraFilterParentPeakMower", "Identification");
    tools_map["SpectraFilterScaler"] = Internal::ToolDescription("SpectraFilterScaler", "Identification");
    tools_map["SpectraFilterSqrtMower"] = Internal::ToolDescription("SpectraFilterSqrtMower", "Identification");
    tools_map["SpectraFilterThresholdMower"] = Internal::ToolDescription("SpectraFilterThresholdMower", "Identification");
    tools_map["SpectraFilterWindowMower"] = Internal::ToolDescription("SpectraFilterWindowMower", "Identification");
    tools_map["SpectraMerger"] = Internal::ToolDescription("SpectraMerger", "Signal processing and preprocessing");
    tools_map["TextExporter"] = Internal::ToolDescription("TextExporter", "File Handling");
    tools_map["TOFCalibration"] = Internal::ToolDescription("TOFCalibration", "Signal processing and preprocessing");
    tools_map["XTandemAdapter"] = Internal::ToolDescription("XTandemAdapter", "Identification");
    // STOP! insert your tool in alphabetical order for easier maintenance (only tools requiring the GUI lib should be added below)

    // ATTENTION: tools requiring the GUI lib
#ifdef WITH_GUI
    tools_map["ExecutePipeline"] = Internal::ToolDescription("ExecutePipeline", "Misc");
    tools_map["Resampler"] = Internal::ToolDescription("Resampler", "Signal processing and preprocessing");
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

  ToolListType ToolHandler::getUtilList()
  {
    ToolListType util_map;
    const String util_category = "Utilities";

    util_map["AccurateMassSearch"] = Internal::ToolDescription("AccurateMassSearch", util_category);
    util_map["CVInspector"] = Internal::ToolDescription("CVInspector", util_category);
    util_map["DecoyDatabase"] = Internal::ToolDescription("DecoyDatabase", util_category);
    util_map["DatabaseFilter"]= Internal::ToolDescription("DatabaseFilter", util_category);
    util_map["DeMeanderize"] = Internal::ToolDescription("DeMeanderize", util_category);
    util_map["Digestor"] = Internal::ToolDescription("Digestor", util_category);
    util_map["DigestorMotif"] = Internal::ToolDescription("DigestorMotif", util_category);
    util_map["ERPairFinder"] = Internal::ToolDescription("ERPairFinder", util_category);
    util_map["FFEval"] = Internal::ToolDescription("FFEval", util_category);
    util_map["FuzzyDiff"] = Internal::ToolDescription("FuzzyDiff", util_category);
    util_map["IDDecoyProbability"] = Internal::ToolDescription("IDDecoyProbability", util_category);
    util_map["IDExtractor"] = Internal::ToolDescription("IDExtractor", util_category);
    util_map["IDMassAccuracy"] = Internal::ToolDescription("IDMassAccuracy", util_category);
    util_map["IDScoreSwitcher"] = Internal::ToolDescription("IDScoreSwitcher", util_category);
    util_map["IDSplitter"] = Internal::ToolDescription("IDSplitter", util_category);
    util_map["LabeledEval"] = Internal::ToolDescription("LabeledEval", util_category);
    util_map["LowMemPeakPickerHiRes"] = Internal::ToolDescription("LowMemPeakPickerHiRes", util_category);
    util_map["LowMemPeakPickerHiRes_RandomAccess"] = Internal::ToolDescription("LowMemPeakPickerHiRes_RandomAccess", util_category);
    util_map["MassCalculator"] = Internal::ToolDescription("MassCalculator", util_category);
    util_map["MetaboliteAdductDecharger"] = Internal::ToolDescription("MetaboliteAdductDecharger", util_category);
    util_map["MetaboliteSpectralMatcher"] = Internal::ToolDescription("MetaboliteSpectralMatcher", util_category);
    util_map["MetaProSIP"] = Internal::ToolDescription("MetaProSIP", util_category);
    util_map["MRMTransitionGroupPicker"] = Internal::ToolDescription("MRMTransitionGroupPicker", util_category);
    util_map["MRMPairFinder"] = Internal::ToolDescription("MRMPairFinder", util_category);
    util_map["MSSimulator"] = Internal::ToolDescription("MSSimulator", util_category);
    util_map["MultiplexResolver"] = Internal::ToolDescription("MultiplexResolver", util_category);
    util_map["MzMLSplitter"] = Internal::ToolDescription("MzMLSplitter", util_category);
    util_map["OpenPepXL"] = Internal::ToolDescription("OpenPepXL", util_category);
    util_map["OpenPepXLLF"] = Internal::ToolDescription("OpenPepXLLF", util_category);
    util_map["OpenSwathWorkflow"] = Internal::ToolDescription("OpenSwathWorkflow", util_category);
    util_map["OpenSwathRewriteToFeatureXML"] = Internal::ToolDescription("OpenSwathRewriteToFeatureXML", "Targeted Experiments");
    util_map["OpenSwathFileSplitter"] = Internal::ToolDescription("OpenSwathFileSplitter", "Targeted Experiments");
    util_map["OpenSwathDIAPreScoring"] = Internal::ToolDescription("OpenSwathDIAPreScoring", "Targeted Experiments");
    util_map["OpenSwathMzMLFileCacher"] = Internal::ToolDescription("OpenSwathMzMLFileCacher", "Targeted Experiments");
    util_map["PeakPickerIterative"] = Internal::ToolDescription("PeakPickerIterative", "Signal processing and preprocessing");
    util_map["TargetedFileConverter"] = Internal::ToolDescription("TargetedFileConverter", "Targeted Experiments");
    //util_map["PeakPickerRapid"] = Internal::ToolDescription("PeakPickerRapid", "Signal processing and preprocessing");
    util_map["PSMFeatureExtractor"] = Internal::ToolDescription("PSMFeatureExtractor", util_category);
    util_map["QCCalculator"] = Internal::ToolDescription("QCCalculator", util_category);
    util_map["QCEmbedder"] = Internal::ToolDescription("QCEmbedder", util_category);
    util_map["QCExtractor"] = Internal::ToolDescription("QCExtractor", util_category);
    util_map["QCExporter"] = Internal::ToolDescription("QCExporter", util_category);
    util_map["QCImporter"] = Internal::ToolDescription("QCImporter", util_category);
    util_map["QCMerger"] = Internal::ToolDescription("QCMerger", util_category);
    util_map["QCShrinker"] = Internal::ToolDescription("QCExporter", util_category);
    util_map["RNADigestor"] = Internal::ToolDescription("RNADigestor", util_category);
    util_map["RNPxl"] = Internal::ToolDescription("RNPxl", util_category);
    util_map["RNPxlSearch"] = Internal::ToolDescription("RNPxlSearch", util_category);
    util_map["RNPxlXICFilter"] = Internal::ToolDescription("RNPxlXICFilter", util_category);
    util_map["RTEvaluation"] = Internal::ToolDescription("RTEvaluation", util_category);
    util_map["SemanticValidator"] = Internal::ToolDescription("SemanticValidator", util_category);
    util_map["SequenceCoverageCalculator"] = Internal::ToolDescription("SequenceCoverageCalculator", util_category);
    util_map["SpecLibCreator"] = Internal::ToolDescription("SpecLibCreator", util_category);
    util_map["SpectraSTSearchAdapter"] = Internal::ToolDescription("SpectraSTSearchAdapter", util_category);
    util_map["SimpleSearchEngine"] = Internal::ToolDescription("SimpleSearchEngine", util_category);
    util_map["SiriusAdapter"] = Internal::ToolDescription("SiriusAdapter", util_category);
    util_map["SvmTheoreticalSpectrumGeneratorTrainer"] = Internal::ToolDescription("SvmTheoreticalSpectrumGeneratorTrainer", util_category);
    util_map["TICCalculator"] = Internal::ToolDescription("TICCalculator", util_category);
    util_map["TransformationEvaluation"] = Internal::ToolDescription("TransformationEvaluation", util_category);
    util_map["XFDR"] = Internal::ToolDescription("XFDR", util_category);
    util_map["XMLValidator"] = Internal::ToolDescription("XMLValidator", util_category);
    // STOP! insert your tool in alphabetical order for easier maintenance (only tools requiring the GUI lib should be added below)

    // ATTENTION: tools requiring the GUI lib
#ifdef WITH_GUI
    util_map["IDEvaluator"] = Internal::ToolDescription("IDEvaluator", util_category);
    util_map["ImageCreator"] = Internal::ToolDescription("ImageCreator", util_category);
    util_map["INIUpdater"] = Internal::ToolDescription("INIUpdater", util_category);
#endif

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
      if (toolname == "GenericWrapper")
        tools = getTOPPToolList(true);
      else
        tools = getTOPPToolList();

      if (tools.has(toolname))
      {
        return tools[toolname].types;
      }
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
          tools_external_ = tools[i_t]; // init
        else
          tools_external_.append(tools[i_t]); // append
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
    if (getenv("OPENMS_TTD_PATH") != 0)
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
    if (getenv("OPENMS_TTD_INTERNAL_PATH") != 0)
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
