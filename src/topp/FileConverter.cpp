// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Marc Sturm, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/IBSpectraFile.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>


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
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FileConverter \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_GenericWrapper (e.g. for calling external converters) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on the output format</td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any vendor software exporting supported formats (e.g. mzML) </td>
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

namespace OpenMS {
enum IMType
{
  IM_NONE, ///< no ion mobility
  IM_STACKED, ///< ion mobility frame is stacked in a single spectrum
  IM_MULTIPLE_SPECTRA ///< ion mobility is recorded as multiple spectra per frame
};

IMType determineIMType(const MSExperiment& exp)
{
  for (Size k = 0; k < exp.size(); k++)
  {
    if (!exp[k].getFloatDataArrays().empty() && 
        exp[k].getFloatDataArrays()[0].getName().find("Ion Mobility") == 0)
    {
      return IMType::IM_STACKED;
    }
    else if (exp[k].getDriftTime() >= 0.0) 
    {
      return IMType::IM_MULTIPLE_SPECTRA;
    }
  }
  return IMType::IM_NONE;
}

/// Process a stack of drift time spectra
void processDriftTimeStack(const std::vector<MSSpectrum>& stack, std::vector<MSSpectrum>& result)
{
  OPENMS_PRECONDITION(!stack.empty(), "Stack cannot be empty")

  // either no drift time or different RT!
  if (stack[0].getDriftTime() >= 0.0)
  {
    // copy meta data without the raw data and without the IM array
    MSSpectrum new_spec = stack[0]; // there is always one spectrum on the stack
    new_spec.clear(false);
    std::vector<OpenMS::DataArrays::FloatDataArray> empty;
    new_spec.setFloatDataArrays(empty);

    OpenMS::DataArrays::FloatDataArray fda;
    String name = "Ion Mobility";
    if (new_spec.getDriftTimeUnit() == MSSpectrum::DriftTimeUnit::MILLISECOND)
    {
      name += " (MS:1002476)";
    }
    else if (new_spec.getDriftTimeUnit() == MSSpectrum::DriftTimeUnit::VSSC)
    {
      name += " (MS:1002815)";
    }
    fda.setName(name);
    for (const auto& s : stack)
    {
      new_spec.insert(new_spec.end(), s.begin(), s.end());
      fda.insert(fda.end(), s.size(), s.getDriftTime());
    }
    new_spec.setFloatDataArrays({fda});
    new_spec.setDriftTime(-1); // drift time is now encoded in the FloatDataArray
    new_spec.setDriftTimeUnit(MSSpectrum::DriftTimeUnit::NONE); // drift time is now encoded in the FloatDataArray
    result.push_back(new_spec);
  }
  else
  {
    // no drift time for these spectra, simply append to the result
    result.insert(result.end(), stack.begin(), stack.end());
  }

}

/**
  @brief Expands a single MSSpectrum (single frame) into individual ion mobility spectrum

  @param tmps The input spectrum (a single spectrum per frame)
  @param result The output spectra with multiple spectra per frame
*/
void expandIMSpectrum(const MSSpectrum& tmps, std::vector<MSSpectrum>& result)
{
  // copy meta data without the raw data and without the IM array
  MSSpectrum settings = tmps;
  settings.clear(false);
  std::vector<OpenMS::DataArrays::FloatDataArray> empty;
  settings.setFloatDataArrays(empty);

  double IM_BINNING = 1e5;

  // Fill temporary spectral map (mobility -> Spectrum) with data from current spectrum
  std::map< int, MSSpectrum > im_map;
  auto im_arr = tmps.getFloatDataArrays()[0];

  // Determine unit name
  String im_name = im_arr.getName();
  auto unit = MSSpectrum::DriftTimeUnit::MILLISECOND;
  if (im_name == "Ion Mobility (MS:1002476)")
  {
    unit = MSSpectrum::DriftTimeUnit::MILLISECOND;
  }
  else if (im_name == "Ion Mobility (MS:1002815)")
  {
    unit = MSSpectrum::DriftTimeUnit::VSSC;
  }

  for (Size k = 0;  k < tmps.size(); k++)
  {
    double im = im_arr[ k ];
    if (im_map.find( int(im*IM_BINNING) ) == im_map.end() )
    {
      // use meta data from combined spectrum, set new name and current drift time
      MSSpectrum news = settings;
      news.setDriftTime(im);
      news.setDriftTimeUnit(unit);
      news.setName(tmps.getName() + "_combined"); // we will not recover original scan ids
      im_map[ int(im*IM_BINNING) ] = news;
    }
    im_map[ int(im*IM_BINNING) ].push_back( tmps[k] );
  }

  // Add spectra to result, note that this is guaranteed to be
  // sorted by ion mobility (through std::map).
  for (const auto& s : im_map)
  {
    result.push_back(s.second);
  }
}

/**
  @brief Collapses multiple IM spectra from the same frame into a single MSSpectrum

  @param exp The input experiment with multiple spectra per frame
  @param result The output spectra collapsed to a single spectrum per frame

  @note: this requires that all spectra from the same frame have the same RT ("scan start time")
*/
void collapseIMSpectrum(const MSExperiment& exp, std::vector<MSSpectrum>& result)
{
  if (exp.empty()) return;

  std::vector<MSSpectrum> stack;
  double curr_rt = exp[0].getRT();
  stack.push_back(exp[0]);
  for (Size k = 1; k < exp.size(); k++)
  {
    // spectra from the same frame (drift time scan) will have the same retention time start time
    if (exp[k].getDriftTime() >= 0.0 && (fabs(exp[k].getRT() - curr_rt) < 1e-5) )
    {
      stack.push_back(exp[k]);
    }
    else
    {
      processDriftTimeStack(stack, result);
      // push next spectrum
      stack.clear();
      stack.push_back(exp[k]);
      curr_rt = exp[k].getRT();
    }
  }

  // stack will always have at least one spectrum
  processDriftTimeStack(stack, result);
}

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
    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false, true); // for TOPPAS
    String formats("mzML,mzXML,mgf,raw,cachedMzML,mzData,dta,dta2d,featureXML,consensusXML,ms2,fid,tsv,peplist,kroenik,edta");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));
    
    registerStringOption_("UID_postprocessing", "<method>", "ensure", "unique ID post-processing for output data.\n'none' keeps current IDs even if invalid.\n'ensure' keeps current IDs but reassigns invalid ones.\n'reassign' assigns new unique IDs.", false, true);
    String method("none,ensure,reassign");
    setValidStrings_("UID_postprocessing", ListUtils::create<String>(method));

    formats = "mzData,mzXML,mzML,cachedMzML,dta2d,mgf,featureXML,consensusXML,edta,csv";
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense.", false, true);
    setValidStrings_("out_type", ListUtils::create<String>(formats));
    registerFlag_("TIC_DTA2D", "Export the TIC instead of the entire experiment in mzML/mzData/mzXML -> DTA2D conversions.", true);
    registerFlag_("MGF_compact", "Use a more compact format when writing MGF (no zero-intensity peaks, limited number of decimal places)", true);
    registerFlag_("force_MaxQuant_compatibility", "[mzXML output only] Make sure that MaxQuant can read the mzXML and set the msManufacturer to 'Thermo Scientific'.", true);
    registerFlag_("convert_to_chromatograms", "[mzML output only] Assumes that the provided spectra represent data in SRM mode or targeted MS1 mode and converts them to chromatogram data.", true);
    registerFlag_("force_TPP_compatibility", "[mzML output only] Make sure that TPP parsers can read the mzML and the precursor ion m/z in the file (otherwise it will be set to zero by the TPP).", true);

    registerStringOption_("change_im_format", "<toogle>", "none", "[mzML output only] How to store ion mobility scans (none: no change in format, multiple: store each frame as multiple scans, one per drift time value, single: store whole frame as single scan with IM values in a FloatDataArray", false, true);
    setValidStrings_("change_im_format", ListUtils::create<String>("none,multiple,single"));

    registerStringOption_("write_scan_index", "<toogle>", "true", "Append an index when writing mzML or mzXML files. Some external tools might rely on it.", false, true);
    setValidStrings_("write_scan_index", ListUtils::create<String>("true,false"));
    registerFlag_("lossy_compression", "Use numpress compression to achieve optimally small file size using linear compression for m/z domain and slof for intensity and float data arrays (attention: may cause small loss of precision; only for mzML data).", true);
    registerDoubleOption_("lossy_mass_accuracy", "<error>", -1.0, "Desired (absolute) m/z accuracy for lossy compression (e.g. use 0.0001 for a mass accuracy of 0.2 ppm at 500 m/z, default uses -1.0 for maximal accuracy).", false, true);

    registerFlag_("process_lowmemory", "Whether to process the file on the fly without loading the whole file into memory first (only for conversions of mzXML/mzML to mzML).\nNote: this flag will prevent conversion from spectra to chromatograms.", true);
    registerInputFile_("NET_executable", "<executable>", "", "The .NET framework executable. Only required on linux and mac.", false, true, ListUtils::create<String>("skipexists"));
    registerInputFile_("ThermoRaw_executable", "<file>", "ThermoRawFileParser.exe", "The ThermoRawFileParser executable.", false, true, ListUtils::create<String>("skipexists"));
    registerFlag_("no_peak_picking", "Disables vendor peak picking for raw files.", true);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file names
    String in = getStringOption_("in");
    bool write_scan_index = getStringOption_("write_scan_index") == "true" ? true : false;
    String store_im = getStringOption_("change_im_format");
    bool force_MaxQuant_compatibility = getFlag_("force_MaxQuant_compatibility");
    bool force_TPP_compatibility = getFlag_("force_TPP_compatibility");
    bool convert_to_chromatograms = getFlag_("convert_to_chromatograms");
    bool lossy_compression = getFlag_("lossy_compression");
    double mass_acc = getDoubleOption_("lossy_mass_accuracy");
    bool no_peak_picking = getFlag_("no_peak_picking");

    //input file type
    FileHandler fh;
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

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

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }


    // output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    bool TIC_DTA2D = getFlag_("TIC_DTA2D");
    bool process_lowmemory = getFlag_("process_lowmemory");

    writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 1);

    String uid_postprocessing = getStringOption_("UID_postprocessing");
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    typedef PeakMap MSExperimentType;
    MSExperimentType exp;

    typedef MSExperimentType::SpectrumType SpectrumType;

    typedef FeatureMap FeatureMapType;

    FeatureMapType fm;
    ConsensusMap cm;

    writeDebug_(String("Loading input file"), 1);

    if (in_type == FileTypes::CONSENSUSXML)
    {
      ConsensusXMLFile().load(in, cm);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting consensus features to peaks. You will lose information!");
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
      writeLog_("RawFileReader reading tool. Copyright 2016 by Thermo Fisher Scientific, Inc. All rights reserved");
      String net_executable = getStringOption_("NET_executable");
      TOPPBase::ExitCodes exit_code;
      QStringList arguments;
#ifdef OPENMS_WINDOWSPLATFORM      
      if (net_executable.empty())
      { // default on Windows: if no mono executable is set use the "native" .NET one
        arguments << String("-i=" + in).toQString()
                  << String("--output_file=" + out).toQString()
                  << String("-f=2").toQString() // indexedMzML
                  << String("-e").toQString(); // ignore instrument errors
        if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
        exit_code = runExternalProcess_(getStringOption_("ThermoRaw_executable").toQString(), arguments);
      }
      else
      { // use e.g., mono
        arguments << getStringOption_("ThermoRaw_executable").toQString()
                  << String("-i=" + in).toQString()
                  << String("--output_file=" + out).toQString()
                  << String("-f=2").toQString()
                  << String("-e").toQString();
        if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
        exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
      }      
#else
      // default on Mac, Linux: use mono
      net_executable = net_executable.empty() ? "mono" : net_executable;
      arguments << getStringOption_("ThermoRaw_executable").toQString()
                << String("-i=" + in).toQString()
                << String("--output_file=" + out).toQString()
                << String("-f=2").toQString()
                << String("-e").toQString();
      if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
      exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
#endif            
      return exit_code;
    }
    else if (in_type == FileTypes::EDTA)
    {
      EDTAFile().load(in, cm);
      cm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You you will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting consensus features to peaks. You will lose information!");
        exp.set2DData(cm);
      }
    }
    else if (in_type == FileTypes::FEATUREXML ||
             in_type == FileTypes::TSV ||
             in_type == FileTypes::PEPLIST ||
             in_type == FileTypes::KROENIK)
    {
      fh.loadFeatures(in, fm, in_type);
      fm.sortByPosition();
      if ((out_type != FileTypes::FEATUREXML) &&
          (out_type != FileTypes::CONSENSUSXML))
      {
        // You will lose information and waste memory. Enough reasons to issue a warning!
        writeLog_("Warning: Converting features to peaks. You will lose information! Mass traces are added, if present as 'num_of_masstraces' and 'masstrace_intensity' (X>=0) meta values.");
        exp.set2DData<true>(fm);
      }
    }
    else if (in_type == FileTypes::CACHEDMZML)
    {
      // Determine location of meta information (empty mzML)
      String in_meta = extractCachedMetaFilename(in);
      if (in_meta.empty()) return ILLEGAL_PARAMETERS;

      MzMLFile f;
      f.setLogType(log_type_);
      Internal::CachedMzMLHandler cacher;
      cacher.setLogType(log_type_);
      PeakMap tmp_exp;

      f.load(in_meta, exp);
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

      if (store_im != "none")
      {
        std::cout << "Converting IM formats is currently not implemented for low-memory processing" << std::endl;
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

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
        if (out_meta.empty()) return ILLEGAL_PARAMETERS;

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
    else
    {
      fh.loadExperiment(in, exp, in_type, log_type_, true, true);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    writeDebug_(String("Writing output file"), 1);

    if (out_type == FileTypes::MZML)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZML));
      MzMLFile f;
      f.setLogType(log_type_);
      f.getOptions().setWriteIndex(write_scan_index);
      f.getOptions().setForceTPPCompatability(force_TPP_compatibility);
      // numpress compression
      if (lossy_compression)
      {
        f.getOptions().setNumpressConfigurationMassTime(npconfig_mz);
        f.getOptions().setNumpressConfigurationIntensity(npconfig_int);
        f.getOptions().setNumpressConfigurationFloatDataArray(npconfig_fda);
        f.getOptions().setCompression(true);
      }

      if (convert_to_chromatograms)
      {
        for (auto & s : exp)
        {
          s.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
        }
      }

      if (store_im != "none")
      {
        IMType itype = determineIMType(exp);

        if (itype == IMType::IM_NONE)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Requested conversion to different ion mobility format, but no ion mobility data is present.");
        }
        else if (store_im == "multiple" && itype == IMType::IM_MULTIPLE_SPECTRA)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Requested conversion to 'multiple' ion mobility format, but data is already in this format.");
        }
        else if (store_im == "single" && itype == IMType::IM_STACKED)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Requested conversion to 'single' ion mobility format, but data is already in this format.");
        }

        if (store_im == "multiple" && itype == IMType::IM_STACKED)
        {
          std::vector<MSSpectrum> result;
          for (Size k = 0; k < exp.size(); k++)
          {
            // For data without ion mobility, simply append the result (only
            // collapse for scans that actually have a float data array).
            if (exp[k].getFloatDataArrays().empty())
            {
              result.push_back(exp[k]);
            }
            else
            {
              expandIMSpectrum(exp[k], result);
            }
          }
          exp.setSpectra(result); // swap data
        }
        else if (store_im == "single" && !exp.empty() && itype == IMType::IM_MULTIPLE_SPECTRA)
        {
          std::vector<MSSpectrum> result;
          collapseIMSpectrum(exp, result);
          exp.setSpectra(result); // swap data
        }
      }
      ChromatogramTools().convertSpectraToChromatograms(exp, true, convert_to_chromatograms);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::MZDATA)
    {
      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZDATA));
      MzDataFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::MZXML)
    {
      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 CONVERSION_MZXML));
      MzXMLFile f;
      f.setLogType(log_type_);
      f.getOptions().setForceMQCompatability(force_MaxQuant_compatibility);
      f.getOptions().setWriteIndex(write_scan_index);
      //ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
      f.store(out, exp);
    }
    else if (out_type == FileTypes::DTA2D)
    {
      //add data processing entry
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::
                                                 FORMAT_CONVERSION));
      DTA2DFile f;
      f.setLogType(log_type_);
      ChromatogramTools().convertChromatogramsToSpectra<MSExperimentType>(exp);
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
      else // not loaded as feature map or consensus map
      {
        // The feature specific information is only defaulted. Enough reasons to issue a warning!
        writeLog_("Warning: Converting peaks to features will lead to incomplete features!");
        fm.clear();
        fm.reserve(exp.getSize());
        typedef FeatureMapType::FeatureType FeatureType;
        FeatureType feature;
        feature.setQuality(0, 1); // override default
        feature.setQuality(1, 1); // override default
        feature.setOverallQuality(1); // override default
        for (MSExperimentType::ConstIterator spec_iter = exp.begin();
             spec_iter != exp.end();
             ++spec_iter
             )
        {
          feature.setRT(spec_iter->getRT());
          for (SpectrumType::ConstIterator peak1_iter = spec_iter->begin();
               peak1_iter != spec_iter->end();
               ++peak1_iter
               )
          {
            feature.setMZ(peak1_iter->getMZ());
            feature.setIntensity(peak1_iter->getIntensity());
            feature.setUniqueId();
            fm.push_back(feature);
          }
        }
        fm.updateRanges();
      }

      addDataProcessing_(fm, getProcessingInfo_(DataProcessing::
                                                FORMAT_CONVERSION));
      FeatureXMLFile().store(out, fm);
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
      ConsensusXMLFile().store(out, cm);
    }
    else if (out_type == FileTypes::EDTA)
    {
      if (fm.size() > 0 && cm.size() > 0)
      {
        OPENMS_LOG_ERROR << "Internal error: cannot decide on container (Consensus or Feature)! This is a bug. Please report it!";
        return INTERNAL_ERROR;
      }
      if (fm.size() > 0) EDTAFile().store(out, fm);
      else if (cm.size() > 0) EDTAFile().store(out, cm);
    }
    else if (out_type == FileTypes::CACHEDMZML)
    {
      // Determine output path for meta information (empty mzML)
      String out_meta = extractCachedMetaFilename(out);
      if (out_meta.empty()) return ILLEGAL_PARAMETERS;

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
    else
    {
      writeLog_("Unknown output file type given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
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
