// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
// TODO add handler support for other accss
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>
#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/FORMAT/IBSpectraFile.h>
// TODO add handler support for other access
#include <OpenMS/FORMAT/MascotGenericFile.h>
// TODO: remove MZML header after we get cached and Transform working
#include <OpenMS/FORMAT/MzMLFile.h>
// TODO: remove MZXML header after we get cached and Transform working
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/ConversionHelper.h>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FileConverter FileConverter

@brief Converts between different MS file formats.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; FileConverter &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any vendor software exporting supported formats (e.g. mzML) </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on the output format</td>
</tr>
</table>
</CENTER>

The main use of this tool is to convert data from external sources to the formats used by OpenMS/TOPP.
Maybe most importantly, data from MS experiments in a number of different formats can be converted to mzML,
the canonical file format used by OpenMS/TOPP for experimental data. (mzML is the PSI approved format and
supports traceability of analysis steps.)

Thermo raw files can be converted to mzML using the ThermoRawFileParser provided in the THIRDPARTY folder.
On windows, a recent .NET framwork needs to be installed. On linux and mac, the mono runtime needs to be
present and accessible via the -NET_executable parameter. The path to the ThermoRawFileParser can be set
via the -ThermoRaw_executable option.

For MaxQuant-flavoured mzXML the use of the advanced option '-force_MaxQuant_compatibility' is recommended.

Many different format conversions are supported, and some may be more useful than others. Depending on the
file formats involved, information can be lost during conversion, e.g. when converting featureXML to mzData.
In such cases a warning is shown.

The input and output file types are determined from	the file extensions or from the first few lines of the
files. If file type determination is not possible, the input or output file type has to be given explicitly.

Conversion with the same output as input format is supported. In some cases, this can be helpful to remove
errors from files (e.g. the index), to update file formats to new versions, or to check whether information is lost upon
reading or writing.

Some information about the supported input types:
@ref OpenMS::MzMLFile "mzML"
@ref OpenMS::MzXMLFile "mzXML"
@ref OpenMS::MzDataFile "mzData"
@ref OpenMS::MascotGenericFile "mgf"
@ref OpenMS::MSPGenericFile "msp"
@ref OpenMS::DTA2DFile "dta2d"
@ref OpenMS::DTAFile "dta"
@ref OpenMS::FeatureXMLFile "featureXML"
@ref OpenMS::ConsensusXMLFile "consensusXML"
@ref OpenMS::MS2File "ms2"
@ref OpenMS::XMassFile "fid/XMASS"
@ref OpenMS::MsInspectFile "tsv"
@ref OpenMS::SpecArrayFile "peplist"
@ref OpenMS::KroenikFile "kroenik"
@ref OpenMS::EDTAFile "edta"
@ref OpenMS::SqMassFile "sqmass"
@ref OpenMS::OMSFile "oms"

@note See @ref TOPP_IDFileConverter for similar functionality for protein/peptide identification file formats.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FileConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FileConverter.html
*/


String extractCachedMetaFilename(const String& in)
{
  // Special handling of cached mzML as input types:
  // we expect two paired input files which we should read into exp
  std::vector<String> split_out;
  in.split(".cachedMzML", split_out);
  if (split_out.size() != 2)
  {
    OPENMS_LOG_ERROR << "Cannot deduce base path from input '" << in
      << "' (note that '.cachedMzML' should only occur once as the final ending)" << std::endl;
    return "";
  }
  String in_meta = split_out[0] + ".mzML";
  return in_meta;
}

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileConverter :
  public TOPPBase
{
public:
  TOPPFileConverter() :
    TOPPBase("FileConverter", "Converts between different MS file formats.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.");
    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false, false); // optional and not advanced (for workflow engines to show this param)
    vector<String> input_formats = {"mzML", "mzXML", "mgf", "msp", "raw", "cachedMzML", "mzData", "dta", "dta2d", "featureXML", "consensusXML", "ms2", "fid",
#ifdef WITH_OPENTIMS
    "d",
#endif
    "tsv", "peplist", "kroenik", "edta", "oms", "sqMass"};
    setValidFormats_("in", input_formats);
    setValidStrings_("in_type", input_formats);

    registerStringOption_("UID_postprocessing", "<method>", "ensure", "unique ID post-processing for output data.\n'none' keeps current IDs even if invalid.\n'ensure' keeps current IDs but reassigns invalid ones.\n'reassign' assigns new unique IDs.", false, true);
    String method("none,ensure,reassign");
    setValidStrings_("UID_postprocessing", ListUtils::create<String>(method));

    vector<String> output_formats = {"mzML", "mzXML", "cachedMzML", "mgf", "msp", "featureXML", "consensusXML", "edta", "mzData", "dta2d", "csv", "sqMass", "xic", "oms"};
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", output_formats);
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense.", false, false); // optional and not advanced (for workflow engines to show this param)
    setValidStrings_("out_type", output_formats);
    registerFlag_("TIC_DTA2D", "Export the TIC instead of the entire experiment in mzML/mzData/mzXML -> DTA2D conversions.", true);
    registerFlag_("MGF_compact", "Use a more compact format when writing MGF (no zero-intensity peaks, limited number of decimal places)", true);
    registerFlag_("force_MaxQuant_compatibility", "[mzXML output only] Make sure that MaxQuant can read the mzXML and set the msManufacturer to 'Thermo Scientific'.", true);
    registerFlag_("force_TPP_compatibility", "[mzML output only] Make sure that TPP parsers can read the mzML and the precursor ion m/z in the file (otherwise it will be set to zero by the TPP).", true);
    registerFlag_("convert_to_chromatograms", "[mzML output only] Assumes that the provided spectra represent data in SRM mode or targeted MS1 mode and converts them to chromatogram data.", true);

    registerStringOption_("write_scan_index", "<toggle>", "true", "Append an index when writing mzML or mzXML files. Some external tools might rely on it.", false, true);
    setValidStrings_("write_scan_index", ListUtils::create<String>("true,false"));
    registerFlag_("lossy_compression", "Use numpress compression to achieve optimally small file size using linear compression for m/z domain and slof for intensity and float data arrays (attention: may cause small loss of precision; only for mzML data).", true);
    registerDoubleOption_("lossy_mass_accuracy", "<error>", -1.0, "Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).", false, true);

    registerFlag_("process_lowmemory", "Whether to process the file on the fly without loading the whole file into memory first (only for conversions of mzXML/mzML to mzML).\nNote: this flag will prevent conversion from spectra to chromatograms.", true);
    
#ifdef WITH_OPENTIMS
    registerTOPPSubsection_("bruker", "Options for reading Bruker TimsTOF .d files (requires WITH_OPENTIMS)");
    registerDoubleOption_("bruker:calibration_tolerance", "<float>", 0.0, "m/z recalibration tolerance (0 = library default)", false, true);
    setMinFloat_("bruker:calibration_tolerance", 0.0);
    registerStringOption_("bruker:calibrate", "<toggle>", "false", "Enable m/z recalibration (may fail on some datasets)", false, true);
    setValidStrings_("bruker:calibrate", {"true", "false"});
    registerStringOption_("bruker:load_ms1", "<toggle>", "true",
      "Load MS1 spectra. Disable for MS2-only workflows (peptide database search) "
      "where MS1 surveys are not needed — substantially cuts memory and time. Affects all export modes.", false, true);
    setValidStrings_("bruker:load_ms1", {"true", "false"});
    registerStringOption_("bruker:export_mode", "<mode>", "auto", "Export mode: 'auto' detects DDA/DIA acquisition type, "
      "'spectrum' forces per-precursor MS2 spectra (DDA-style), 'frame' returns raw 4D frames without signal processing.", false, true);
    setValidStrings_("bruker:export_mode", {"auto", "spectrum", "frame"});
    registerDoubleOption_("bruker:ms1_centroid_mz_ppm", "<float>", 0.0,
      "MS1 frame IM-centroiding m/z tolerance in ppm. Collapses the ion mobility dimension "
      "by aggregating neighboring peaks. Both this and ms1_centroid_im_pct must be > 0 to enable. "
      "Suggested value: 5.0. Algorithm from Sage (Lazear 2023).", false, true);
    setMinFloat_("bruker:ms1_centroid_mz_ppm", 0.0);
    registerDoubleOption_("bruker:ms1_centroid_im_pct", "<float>", 0.0,
      "MS1 frame IM-centroiding ion mobility tolerance in percent. Both this and ms1_centroid_mz_ppm "
      "must be > 0 to enable. Suggested value: 3.0.", false, true);
    setMinFloat_("bruker:ms1_centroid_im_pct", 0.0);
    registerIntOption_("bruker:dia_ms2_n_neighbors", "<int>", 0,
      "DIA MS2 frame aggregation: number of adjacent frames on each side to sum per SWATH window. "
      "0 = disabled (raw export), 1 = 3-frame sum, 2 = 5-frame sum. "
      "Boosts signal by summing intensity across neighboring RT frames, then removes isolated noise.", false, true);
    setMinInt_("bruker:dia_ms2_n_neighbors", 0);
    registerIntOption_("bruker:dia_ms2_min_support", "<int>", 1,
      "DIA MS2 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point "
      "(center cell excluded from count). Applied after frame aggregation. Only effective when dia_ms2_n_neighbors > 0. "
      "Set to 0 to disable denoising (useful for pure centroiding without noise filtering).", false, true);
    setMinInt_("bruker:dia_ms2_min_support", 0);
    registerStringOption_("bruker:dia_ms2_centroid", "<toggle>", "false",
      "Apply 2D Gaussian smoothing + local maxima peak picking to the denoised DIA MS2 grid. "
      "Produces IM_CENTROIDED spectra with sub-bin (m/z, IM) precision. Only effective when dia_ms2_n_neighbors > 0.", false, true);
    setValidStrings_("bruker:dia_ms2_centroid", {"true", "false"});

    registerIntOption_("bruker:ms1_n_neighbors", "<int>", 0,
      "MS1 frame aggregation: number of adjacent MS1 frames on each side to sum. "
      "0 = disabled (raw export), 1 = 3-frame sum, 2 = 5-frame sum. "
      "Applies to both DIA and DDA; ignored in FRAME export mode.", false, true);
    setMinInt_("bruker:ms1_n_neighbors", 0);
    setMaxInt_("bruker:ms1_n_neighbors", 50);

    registerIntOption_("bruker:ms1_min_support", "<int>", 0,
      "MS1 denoising: minimum occupied neighbor cells in a 3x3 (m/z x IM) grid to keep a point. "
      "Applied after aggregation. 0 = disabled, 8 = all 8 neighbors required (strictest). "
      "Only effective when ms1_n_neighbors > 0. Appropriate for dense survey runs; disable for "
      "rare-species discovery.", false, true);
    setMinInt_("bruker:ms1_min_support", 0);
    setMaxInt_("bruker:ms1_min_support", 8);

    registerDoubleOption_("bruker:ms1_max_rt_distance_sec", "<float>", 0.0,
      "Cap the RT distance (seconds) between a neighbor MS1 frame and the center frame during "
      "aggregation. 0.0 = no cap. Recommended for DDA (e.g. 5.0) where MS1 frame cadence is "
      "irregular. The center frame is always included regardless of this cap.", false, true);
    setMinFloat_("bruker:ms1_max_rt_distance_sec", 0.0);

    registerIntOption_("bruker:ms1_centroid_max_peaks", "<int>", 100000,
      "Cap on the number of centroided peaks retained per MS1 spectrum. Top-intensity peaks "
      "are kept; low-intensity tail is dropped if the limit is hit (a warning is logged in that "
      "case). Only effective when MS1 centroiding is enabled via ms1_centroid_mz_ppm/pct. Raise "
      "for aggregated MS1 (ms1_n_neighbors > 0) on dense surveys; lower to trim long-tail noise.",
      false, true);
    setMinInt_("bruker:ms1_centroid_max_peaks", 1);
#endif

    registerTOPPSubsection_("RawToMzML", "Options for converting raw files to mzML (uses ThermoRawFileParser)");
    registerInputFile_("RawToMzML:NET_executable", "<executable>", "", "The .NET framework executable. Only required on linux and mac.", false, true, {"is_executable"});
    registerInputFile_("RawToMzML:ThermoRaw_executable", "<file>", "ThermoRawFileParser.exe", "The ThermoRawFileParser executable.", false, true, {"is_executable"});
    setValidFormats_("RawToMzML:ThermoRaw_executable", {"exe"});
    registerFlag_("RawToMzML:no_peak_picking", "Disables vendor peak picking for raw files.", true);
    registerFlag_("RawToMzML:no_zlib_compression", "Disables zlib compression for raw file conversion. Enables compatibility with some tools that do not support compressed input files, e.g. X!Tandem.", true);
    registerFlag_("RawToMzML:include_noise", "Include noise data in mzML output.", true);
    
  // OpenSwath / chromatogram options: allow passing a transition library to map extracted ion chromatograms to their matching metadata in the transition list
  registerTOPPSubsection_("OpenSwathWorkflow", "Options for loading OpenSWATH transition libraries used for chromatogram metadata");
  registerInputFile_("OpenSwathWorkflow:tr", "<file>", "", "Transition library (PQP, TSV, or TraML) providing precursor/transition metadata. Required when converting sqMass to CHROMPARQUET (.xic); XICs without associated metadata are not meaningful.", false, true);
  setValidFormats_("OpenSwathWorkflow:tr", {"pqp", "tsv", "traml", "osw"});
  registerStringOption_("OpenSwathWorkflow:tr_type", "<type>", "", "Type hint for the transition file (pqp, tsv, traml). If not provided, the type is inferred from the file extension.", false, true);
  registerFlag_("OpenSwathWorkflow:legacy_traml_id", "When loading PQP libraries: use legacy TraML IDs (TRAML_ID) instead of numeric IDs.", true);
  }

  // Note: subsection defaults are not overridden here; TransitionTSVFile parameters
  // are available via the standard parameter mechanism when requested.

#ifdef WITH_OPENTIMS
  BrukerTimsFile::Config getBrukerConfig_()
  {
    BrukerTimsFile::Config c;
    c.calibration_tolerance = getDoubleOption_("bruker:calibration_tolerance");
    c.calibrate = (getStringOption_("bruker:calibrate") == "true");
    c.load_ms1 = (getStringOption_("bruker:load_ms1") == "true");
    String mode = getStringOption_("bruker:export_mode");
    if (mode == "spectrum") c.export_mode = BrukerTimsFile::Config::SPECTRUM;
    else if (mode == "frame") c.export_mode = BrukerTimsFile::Config::FRAME;
    else c.export_mode = BrukerTimsFile::Config::AUTO;
    c.ms1_centroid_mz_ppm = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_mz_ppm"));
    c.ms1_centroid_im_pct = static_cast<float>(getDoubleOption_("bruker:ms1_centroid_im_pct"));
    c.dia_ms2_n_neighbors = getIntOption_("bruker:dia_ms2_n_neighbors");
    c.dia_ms2_min_support = getIntOption_("bruker:dia_ms2_min_support");
    c.dia_ms2_centroid = (getStringOption_("bruker:dia_ms2_centroid") == "true");
    c.ms1_n_neighbors         = getIntOption_("bruker:ms1_n_neighbors");
    c.ms1_min_support         = getIntOption_("bruker:ms1_min_support");
    c.ms1_max_rt_distance_sec = getDoubleOption_("bruker:ms1_max_rt_distance_sec");
    c.ms1_centroid_max_peaks  = getIntOption_("bruker:ms1_centroid_max_peaks");
    return c;
  }
#endif

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file names
    String in = getStringOption_("in");
    bool write_scan_index = getStringOption_("write_scan_index") == "true" ? true : false;
    bool force_MaxQuant_compatibility = getFlag_("force_MaxQuant_compatibility");
    bool force_TPP_compatibility = getFlag_("force_TPP_compatibility");
    bool convert_to_chromatograms = getFlag_("convert_to_chromatograms");
    bool lossy_compression = getFlag_("lossy_compression");
    double mass_acc = getDoubleOption_("lossy_mass_accuracy");

    // prepare data structures for lossy compression (note that we compress any float data arrays the same as intensity arrays)
    MSNumpressCoder::NumpressConfig npconfig_mz, npconfig_int, npconfig_fda;
    npconfig_mz.estimate_fixed_point = true; // critical
    npconfig_int.estimate_fixed_point = true; // critical
    npconfig_fda.estimate_fixed_point = true; // critical
    npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
    npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
    npconfig_fda.numpressErrorTolerance = -1.0; // skip check, faster
    npconfig_mz.setCompression("linear");
    npconfig_int.setCompression("slof");
    npconfig_fda.setCompression("slof");
    npconfig_mz.linear_fp_mass_acc = mass_acc; // set the desired mass accuracy

    // input file type
    FileHandler fh;
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));
    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
      if (in_type == FileTypes::UNKNOWN)
      {
        writeLogError_("Error: Could not determine input file type!");
        return PARSE_ERROR;
      }
    }

    // output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileHandler::getConsistentOutputfileType(out, getStringOption_("out_type"));
    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine output file type! Please adjust the 'out_type' parameter of this tool.");
      return PARSE_ERROR;
    }

    bool TIC_DTA2D = getFlag_("TIC_DTA2D");
    bool process_lowmemory = getFlag_("process_lowmemory");

    writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 1);

    String uid_postprocessing = getStringOption_("UID_postprocessing");
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment exp;
    assert(exp.empty());
    const MSExperiment empty_exp; ///< to determine if 'exp' was modified (loading and storing an MSExp with metadata but empty spectra/chroms should be valid), i.e. checking exp.empty() is not sufficient

    FeatureMap fm;
    ConsensusMap cm;

    writeDebug_(String("Loading input file"), 1);

    if (in_type == FileTypes::CONSENSUSXML)
    {
      FileHandler().loadConsensusFeatures(in, cm, {FileTypes::CONSENSUSXML}, log_type_);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML) &&
          (out_type != FileTypes::OMS)
          )
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLogWarn_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
    }
    else if (in_type == FileTypes::RAW)
    {
      if (out_type != FileTypes::MZML)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Only conversion to mzML supported at this point.");
      }
      bool no_peak_picking = getFlag_("RawToMzML:no_peak_picking");
      bool no_zlib_compression = getFlag_("RawToMzML:no_zlib_compression");
      bool include_noise = getFlag_("RawToMzML:include_noise");
      writeLogInfo_("RawFileReader reading tool. Copyright 2016 by Thermo Fisher Scientific, Inc. All rights reserved");
      String net_executable = getStringOption_("RawToMzML:NET_executable");
      std::vector<String> arguments;
#ifdef OPENMS_WINDOWSPLATFORM
      if (net_executable.empty())
      { // default on Windows: if NO mono executable is set use the "native" .NET one
        net_executable = getStringOption_("RawToMzML:ThermoRaw_executable");
      }
      else
      { // use e.g., mono
        arguments.push_back(getStringOption_("RawToMzML:ThermoRaw_executable"));
      }
#else
      // default on Mac, Linux: use mono
      net_executable = net_executable.empty() ? "mono" : net_executable;
      arguments.push_back(getStringOption_("RawToMzML:ThermoRaw_executable"));
#endif
      arguments.push_back("--input=" + in);
      arguments.push_back("--output=" + out);
      arguments.push_back("-f=2"); // indexedMzML
      arguments.push_back("-e"); // ignore instrument errors
      if (no_peak_picking)
      {
        arguments.push_back("--noPeakPicking");
      }
      if (no_zlib_compression)
      {
        arguments.push_back("--noZlibCompression");
      }
      if (include_noise)
      {
        arguments.push_back("--noiseData");
      }
      return runExternalProcess_(net_executable, arguments);
    }
    else if (in_type == FileTypes::EDTA)
    {
      FileHandler().loadConsensusFeatures(in, cm, {FileTypes::EDTA}, log_type_);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLogWarn_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
    }
    else if (in_type == FileTypes::FEATUREXML ||
             in_type == FileTypes::TSV ||
             in_type == FileTypes::PEPLIST ||
             in_type == FileTypes::KROENIK)
    {
      fh.loadFeatures(in, fm, {in_type});
      fm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML) &&
          (out_type != FileTypes::OMS))
      {
        // You will lose information and waste memory. Enough reasons to issue a warning!
        writeLogWarn_("Warning: Converting features to peaks. You will lose information! Mass traces are added, if present as 'num_of_masstraces' and 'masstrace_intensity' (X>=0) meta values.");
        exp.set2DData<true>(fm);
      }
    }
    else if (in_type == FileTypes::CACHEDMZML)
    {
      // Determine location of meta information (empty mzML)
      String in_meta = extractCachedMetaFilename(in);
      if (in_meta.empty())
      {
        return ILLEGAL_PARAMETERS;
      }
      Internal::CachedMzMLHandler cacher;
      cacher.setLogType(log_type_);
      PeakMap tmp_exp;

      FileHandler().loadExperiment(in_meta, exp, {FileTypes::MZML}, log_type_);
      cacher.readMemdump(tmp_exp, in);

      // Sanity check
      if (exp.size() != tmp_exp.size())
      {
        OPENMS_LOG_ERROR << "Paired input files do not match, cannot convert: " << in_meta << " and " << in << std::endl;
        return ILLEGAL_PARAMETERS;
      }

      // Populate meta data with actual data points
      for (Size i=0; i < tmp_exp.size(); ++i)
      {
        for (Size j = 0; j < tmp_exp[i].size(); j++)
        {
          exp[i].push_back(tmp_exp[i][j]);
        }
      }
      std::vector<MSChromatogram > old_chromatograms = exp.getChromatograms();
      for (Size i=0; i < tmp_exp.getChromatograms().size(); ++i)
      {
        for (Size j = 0; j < tmp_exp.getChromatograms()[i].size(); j++)
        {
          old_chromatograms[i].push_back(tmp_exp.getChromatograms()[i][j]);
        }
      }
      exp.setChromatograms(old_chromatograms);
    }
    else if (process_lowmemory)
    {
      // Special switch for the low memory options:
      // We can transform the complete experiment directly without first
      // loading the complete data into memory. PlainMSDataWritingConsumer will
      // write out mzML to disk as they are read from the input.

      if ((in_type == FileTypes::MZXML || in_type == FileTypes::MZML) && out_type == FileTypes::MZML)
      {
        // Prepare the consumer
        PlainMSDataWritingConsumer consumer(out);
        consumer.getOptions().setWriteIndex(write_scan_index);
        bool skip_full_count = false;
        // numpress compression
        if (lossy_compression)
        {
          consumer.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
          consumer.getOptions().setNumpressConfigurationIntensity(npconfig_int);
          consumer.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);
          consumer.getOptions().setCompression(true);
        }
        consumer.addDataProcessing(getProcessingInfo_(DataProcessing::CONVERSION_MZML));

        // for different input file type
        if (in_type == FileTypes::MZML)
        {
          MzMLFile mzmlfile;
          mzmlfile.setLogType(log_type_);
          mzmlfile.transform(in, &consumer, skip_full_count);
          return EXECUTION_OK;
        }
        else if (in_type == FileTypes::MZXML)
        {
          MzXMLFile mzxmlfile;
          mzxmlfile.setLogType(log_type_);
          mzxmlfile.transform(in, &consumer, skip_full_count);
          return EXECUTION_OK;
        }
      }
      else if (in_type == FileTypes::MZML && out_type == FileTypes::CACHEDMZML)
      {
        // Determine output path for meta information (empty mzML)
        String out_meta = extractCachedMetaFilename(out);
        if (out_meta.empty())
        {
          return ILLEGAL_PARAMETERS;
        }
        Internal::CachedMzMLHandler cacher;
        cacher.setLogType(log_type_);
        PeakMap exp_meta;

        MSDataCachedConsumer consumer(out);
        MzMLFile().transform(in, &consumer, exp_meta);
        cacher.writeMetadata(exp_meta, out_meta);

        return EXECUTION_OK;
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Process_lowmemory option can only be used with mzML / mzXML input and mzML output data types.");
      }
    }
#ifdef WITH_OPENTIMS
    else if (in_type == FileTypes::BRUKER_TDF)
    {
      auto bruker_config = getBrukerConfig_();
      BrukerTimsFile tims_file;
      tims_file.setLogType(log_type_);
      tims_file.load(in, exp, bruker_config);
    }
#endif
    else
    {
      fh.loadExperiment(in, exp, {in_type}, log_type_, true, true);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    writeDebug_(String("Writing output file"), 1);

    if (out_type == FileTypes::MZML)
    {
      if (exp == empty_exp)
      {
        OPENMS_LOG_ERROR << "No input data: no MS1/MS2 data present! Cannot write mzML. Please use another input/output format combination.";
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }

      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZML));
      FileHandler mzmlFile;
      mzmlFile.getOptions().setWriteIndex(write_scan_index);
      mzmlFile.getOptions().setForceTPPCompatability(force_TPP_compatibility);
      // numpress compression
      if (lossy_compression)
      {
        mzmlFile.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
        mzmlFile.getOptions().setNumpressConfigurationIntensity(npconfig_int);
        mzmlFile.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);
        mzmlFile.getOptions().setCompression(true);
      }

      if (convert_to_chromatograms)
      {
        for (auto & s : exp)
        {
          s.getInstrumentSettings().setScanMode(InstrumentSettings::ScanMode::SRM);
        }
      }

      ChromatogramTools().convertSpectraToChromatograms(exp, true, convert_to_chromatograms);
      mzmlFile.storeExperiment(out, exp, {FileTypes::MZML}, log_type_);
    }
    else if (out_type == FileTypes::MZDATA)
    {
      if (exp == empty_exp)
      {
        OPENMS_LOG_ERROR << "No input data: no MS1/MS2 data present! Cannot write mzData. Please use another input/output format combination.";
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }

      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZDATA));
      ChromatogramTools().convertChromatogramsToSpectra<MSExperiment>(exp);
      FileHandler().storeExperiment(out, exp, {FileTypes::MZDATA}, log_type_);
    }
    else if (out_type == FileTypes::MZXML)
    {
      if (exp == empty_exp)
      {
        OPENMS_LOG_ERROR << "No input data: no MS1/MS2 data present! Cannot write mzXML. Please use another input/output format combination.";
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }

      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZXML));
      FileHandler f;
      f.getOptions().setForceMQCompatability(force_MaxQuant_compatibility);
      f.getOptions().setWriteIndex(write_scan_index);
      f.storeExperiment(out, exp, {FileTypes::MZXML}, log_type_);
    }
    else if (out_type == FileTypes::DTA2D)
    {
      if (exp == empty_exp)
      {
        OPENMS_LOG_ERROR << "No input data: no MS1/MS2 data present! Cannot write DTA2D. Please use another input/output format combination.";
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      DTA2DFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperiment>(exp);
      if (TIC_DTA2D)
      {
        // store the total ion chromatogram (TIC)
        f.storeTIC(out, exp);
      }
      else
      {
        // store entire experiment
        f.store(out, exp);
      }


    }
    else if (out_type == FileTypes::MGF)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      MascotGenericFile f;
      f.setLogType(log_type_);
      f.store(out, exp, getFlag_("MGF_compact"));
    }
    else if (out_type == FileTypes::MSP)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      FileHandler().storeExperiment(out, exp, {FileTypes::MSP}, log_type_);
    }
    else if (out_type == FileTypes::FEATUREXML)
    {
      if ((in_type == FileTypes::FEATUREXML) || (in_type == FileTypes::TSV) ||
          (in_type == FileTypes::PEPLIST) || (in_type == FileTypes::KROENIK))
      {
        if (uid_postprocessing == "ensure")
        {
          fm.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
        }
        else if (uid_postprocessing == "reassign")
        {
          fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
        }
      }
      else if (in_type == FileTypes::CONSENSUSXML || in_type == FileTypes::EDTA)
      {
        MapConversion::convert(cm, true, fm);
      }
      else if (in_type == FileTypes::OMS)
      {
        FileHandler().loadFeatures(in, fm, {FileTypes::OMS}, log_type_);
        IdentificationDataConverter::exportFeatureIDs(fm);
      }
      else // not loaded as feature map or consensus map
      {
        // The feature specific information is only defaulted. Enough reasons to issue a warning!
        writeLogWarn_("Warning: Converting peaks to features will lead to incomplete features!");
        fm.clear();
        fm.reserve(exp.getSize());
        Feature feature;
        feature.setQuality(0, 1); // override default
        feature.setQuality(1, 1); // override default
        feature.setOverallQuality(1); // override default
        for (const MSSpectrum& spec : exp)
        {
          feature.setRT(spec.getRT());
          for (const Peak1D& peak : spec)
          {
            feature.setMZ(peak.getMZ());
            feature.setIntensity(peak.getIntensity());
            feature.setUniqueId();
            fm.push_back(feature);
          }
        }
        fm.updateRanges();
      }

      addDataProcessing_(fm, getProcessingInfo_(DataProcessing::
                                                FORMAT_CONVERSION));
      FileHandler().storeFeatures(out, fm, {FileTypes::FEATUREXML}, log_type_);
    }
    else if (out_type == FileTypes::CONSENSUSXML)
    {
      if ((in_type == FileTypes::FEATUREXML) || (in_type == FileTypes::TSV) ||
          (in_type == FileTypes::PEPLIST) || (in_type == FileTypes::KROENIK))
      {
        if (uid_postprocessing == "ensure")
        {
          fm.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
        }
        else if (uid_postprocessing == "reassign")
        {
          fm.applyMemberFunction(&UniqueIdInterface::setUniqueId);
        }
        MapConversion::convert(0, fm, cm);
      }
      // nothing to do for consensus input
      else if (in_type == FileTypes::CONSENSUSXML || in_type == FileTypes::EDTA)
      {
      }
      else // experimental data
      {
        MapConversion::convert(0, exp, cm, exp.size());
      }
      for (auto& pepID : cm.getUnassignedPeptideIdentifications())
      {
        pepID.setMetaValue("map_index", 0);
      }

      addDataProcessing_(cm, getProcessingInfo_(DataProcessing::
                                                FORMAT_CONVERSION));
      FileHandler().storeConsensusFeatures(out, cm, {FileTypes::CONSENSUSXML}, log_type_);
    }
    else if (out_type == FileTypes::EDTA)
    {
      if (!fm.empty() && !cm.empty())
      {
        OPENMS_LOG_ERROR << "Internal error: cannot decide on container (Consensus or Feature)! This is a bug. Please report it!";
        return INTERNAL_ERROR;
      }
      if (fm.empty() && cm.empty())
      {
        OPENMS_LOG_ERROR << "No input data: either Consensus or Feature data present! Cannot write EDTA. Please use another input/output format combination.";
        return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }
      if (!fm.empty())
      {
        FileHandler().storeFeatures(out, fm, {FileTypes::EDTA}, log_type_);
      }
      else if (!cm.empty())
      {
        FileHandler().storeConsensusFeatures(out, cm, {FileTypes::EDTA}, log_type_);
      }
    }
    else if (out_type == FileTypes::CACHEDMZML)
    {
      // Determine output path for meta information (empty mzML)
      String out_meta = extractCachedMetaFilename(out);
      if (out_meta.empty())
      {
        return ILLEGAL_PARAMETERS;
      }
      Internal::CachedMzMLHandler().writeMetadata(exp, out_meta);
      Internal::CachedMzMLHandler().writeMemdump(exp, out);
    }
    else if (out_type == FileTypes::CSV)
    {
      // as ibspectra is currently the only csv/text based format we assume
      // that out_type == FileTypes::CSV means ibspectra, if more formats
      // are added we need a more intelligent strategy to decide which
      // conversion is requested

      // IBSpectra selected as output type
      if (in_type != FileTypes::CONSENSUSXML)
      {
        OPENMS_LOG_ERROR << "Incompatible input data: FileConverter can only convert consensusXML files to ibspectra format.";
        return INCOMPATIBLE_INPUT_DATA;
      }

      IBSpectraFile ibfile;
      ibfile.store(out, cm);
    }
    else if (out_type == FileTypes::CHROMPARQUET)
    {
      // Convert to OpenSWATH Parquet chromatogram file (.xic)
      if (in_type == FileTypes::SQMASS)
      {
        // A transition library is required: XICs without precursor/transition
        // metadata are not meaningful.
        String tr_file = getStringOption_("OpenSwathWorkflow:tr");
        if (tr_file.empty())
        {
          writeLogError_("Error: Converting sqMass to CHROMPARQUET (.xic) requires a transition library "
                         "supplied via -OpenSwathWorkflow:tr. XICs without associated metadata are not meaningful.");
          return ILLEGAL_PARAMETERS;
        }

        // Resolve the file type: honour the explicit type hint first, then
        // fall back to extension/content detection.
        FileHandler fh_local;
        FileTypes::Type tr_type = FileTypes::UNKNOWN;
        const String tr_type_hint = getStringOption_("OpenSwathWorkflow:tr_type");
        if (!tr_type_hint.empty())
        {
          tr_type = FileTypes::nameToType(tr_type_hint);
          if (tr_type == FileTypes::UNKNOWN)
          {
            writeLogError_(String("Error: Unsupported value for -OpenSwathWorkflow:tr_type: '") + tr_type_hint + "'.");
            return ILLEGAL_PARAMETERS;
          }
        }
        else
        {
          tr_type = fh_local.getType(tr_file);
        }

        OpenSwath::LightTargetedExperiment transition_exp;
        try
        {
          if (tr_type == FileTypes::PQP || tr_type == FileTypes::OSW)
          {
            TransitionPQPFile pqp_reader;
            bool legacy = getFlag_("OpenSwathWorkflow:legacy_traml_id");
            pqp_reader.setLogType(log_type_);
            pqp_reader.convertPQPToTargetedExperiment(tr_file.c_str(), transition_exp, legacy);
          }
          else if (tr_file.hasSuffix(".oswpq") || tr_file.hasSuffix(".oswpq.zip"))
          {
            writeLogError_(String("Error: .oswpq libraries are not supported by this build of FileConverter. "
                                  "Cannot convert without transition metadata (file: '") + tr_file + "').");
            return ILLEGAL_PARAMETERS;
          }
          else if (tr_type == FileTypes::TSV)
          {
            TransitionTSVFile tsv_reader;
            Param reader_parameters = getParam_().copy("OpenSwathWorkflow:", true);
            tsv_reader.setLogType(log_type_);
            tsv_reader.setParameters(reader_parameters);
            tsv_reader.convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp);
          }
          else if (tr_type == FileTypes::TRAML)
          {
            TargetedExperiment targeted_exp;
            FileHandler().loadTransitions(tr_file, targeted_exp, {FileTypes::TRAML});
            OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
          }
          else
          {
            writeLogError_(String("Error: Unrecognized transition library type for '") + tr_file +
                           "'. Cannot convert sqMass to CHROMPARQUET without valid transition metadata.");
            return ILLEGAL_PARAMETERS;
          }
        }
        catch (const Exception::BaseException& e)
        {
          writeLogError_(String("Error: Failed to load transition library '") + tr_file + "': " + e.what());
          return ILLEGAL_PARAMETERS;
        }

        MSChromatogramParquetConsumer consumer(out, 0, in, transition_exp);
        // Stream directly from the sqMass file (memory-efficient)
        SqMassFile().transform(in, &consumer, /*skip_full_count=*/false, /*skip_first_pass=*/false);
        consumer.finalize();
      }
      else
      {
        // For non-sqMass inputs, CHROMPARQUET output requires transition metadata.
        // Only sqMass → CHROMPARQUET is currently supported.
        writeLogError_("Error: CHROMPARQUET (.xic) output is only supported when converting from sqMass input. "
                       "Supply an sqMass file via -in.");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
    else if (out_type == FileTypes::SQMASS)
    {
      FileHandler().storeExperiment(out, exp, {FileTypes::SQMASS}, log_type_);
    }
    else if (out_type == FileTypes::OMS)
    {
      if (in_type == FileTypes::FEATUREXML)
      {
        IdentificationDataConverter::importFeatureIDs(fm);
        FileHandler().storeFeatures(out, fm, {FileTypes::OMS}, log_type_);
      }
      else if (in_type == FileTypes::CONSENSUSXML)
      {
        IdentificationDataConverter::importConsensusIDs(cm);
        FileHandler().storeConsensusFeatures(out, cm, {FileTypes::OMS}, log_type_);
      }
      else
      {
        OPENMS_LOG_ERROR << "Incompatible input data: FileConverter can only convert featureXML and consensusXML files to oms format.";
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
    else
    {
      writeLogError_("Error: Unknown output file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // last check if output file was written:
    if (!File::exists(out))
    {
      OPENMS_LOG_ERROR << "Internal error: Conversion did not create an output file! This is a bug. Please report it!";
      return INTERNAL_ERROR;
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFileConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
