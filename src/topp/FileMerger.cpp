// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FileMerger FileMerger

@brief Merges several files. Multiple output formats supported, depending on the input format.

<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; FileMerger &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool/instrument producing mergeable files </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating merged files (e.g. @ref TOPP_XTandemAdapter for mzML, @ref TOPP_ProteinQuantifier for consensusXML) </td>
</tr>
</table>
</center>

Special attention should be given to the append_method for consensusXMLs. One column corresponds to one channel/label + raw file. Rows are quantified and linked features.
More details on the use cases can be found at the parameter description.

For non-consensusXML or consensusXML merging with append_rows, the meta information that is valid for the whole experiment (e.g. MS instrument and sample)
is taken from the first file only.

For spectrum-containing formats (no feature/consensusXML), the retention times for the individual scans are taken from either:
<ul>
<li>the input file meta data (e.g. mzML)
<li>from the input file names (name must contain 'rt' directly followed by a number, e.g. 'myscan_rt3892.98_MS2.dta')
<li>as a list (one RT for each file)
<li>or are auto-generated (starting at 1 with 1 second increment).
</ul>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FileMerger.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FileMerger.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFileMerger :
  public TOPPBase
{
public:

  TOPPFileMerger() :
    TOPPBase("FileMerger", "Merges several MS files into one file."),
    rt_gap_(0.0), rt_offset_(0.0)
  {
  }

protected:

  double rt_gap_, rt_offset_; // parameters for RT concatenation

  void registerOptionsAndFlags_() override
  {
    StringList valid_in = ListUtils::create<String>("mzData,mzXML,mzML,dta,dta2d,mgf,featureXML,consensusXML,fid,traML,fasta");
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", valid_in);
    registerStringOption_("in_type", "<type>", "", "Input file type (default: determined from file extension or content)", false);
    setValidStrings_("in_type", valid_in);
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("mzML,featureXML,consensusXML,traML,fasta"));

    registerFlag_("annotate_file_origin", "Store the original filename in each feature using meta value \"file_origin\" (for featureXML and consensusXML only).");
    registerStringOption_("append_method", "<choice>", "append_rows", "(ConsensusXML-only) Append quantitative information about features row-wise or column-wise.\n"
    "- 'append_rows' is usually used when the inputs come from the same MS run (e.g. caused by manual splitting or multiple algorithms on the same file)\n"
    "- 'append_cols' when you want to combine consensusXMLs from e.g. different fractions to be summarized in ProteinQuantifier or jointly exported with MzTabExporter."
    , false);
    setValidStrings_("append_method", ListUtils::create<String>("append_rows,append_cols"));
    
    registerTOPPSubsection_("rt_concat", "Options for concatenating files in the retention time (RT) dimension. The RT ranges of inputs are adjusted so they don't overlap in the merged file (traML input not supported)");
    registerDoubleOption_("rt_concat:gap", "<sec>", 0.0, "The amount of gap (in seconds) to insert between the RT ranges of different input files. RT concatenation is enabled if a value > 0 is set.", false);
    registerOutputFileList_("rt_concat:trafo_out", "<files>", vector<String>(), "Output of retention time transformations that were applied to the input files to produce non-overlapping RT ranges. If used, one output file per input file is required.", false);
    setValidFormats_("rt_concat:trafo_out", ListUtils::create<String>("trafoXML"));

    registerTOPPSubsection_("raw", "Options for raw data input/output (primarily for DTA files)");
    registerFlag_("raw:rt_auto", "Assign retention times automatically (integers starting at 1)");
    registerDoubleList_("raw:rt_custom", "<rts>", DoubleList(), "List of custom retention times that are assigned to the files. The number of given retention times must be equal to the number of input files.", false);
    registerFlag_("raw:rt_filename", "Try to guess the retention time of a file based on the filename. This option is useful for merging DTA files, where filenames should contain the string 'rt' directly followed by a floating point number, e.g. 'my_spectrum_rt2795.15.dta'");
    registerIntOption_("raw:ms_level", "<num>", 0, "If 1 or higher, this number is assigned to spectra as the MS level. This option is useful for DTA files which do not contain MS level information.", false);
  }

  template <class MapType>
  void adjustRetentionTimes_(MapType& map, const String& trafo_out,
                             bool first_file)
  {
    map.updateRanges();
    TransformationDescription trafo;
    if (first_file) // no transformation necessary
    {
      rt_offset_ = map.getMaxRT() + rt_gap_;
      trafo.fitModel("identity");
    }
    else // subsequent file -> apply transformation
    {
      TransformationDescription::DataPoints points(2);
      double rt_min = map.getMinRT(), rt_max = map.getMaxRT();
      points[0] = make_pair(rt_min, rt_offset_);
      rt_offset_ += rt_max - rt_min;
      points[1] = make_pair(rt_max, rt_offset_);
      trafo.setDataPoints(points);
      trafo.fitModel("linear");
      MapAlignmentTransformer::transformRetentionTimes(map, trafo, true);
      rt_offset_ += rt_gap_;
    }
    if (!trafo_out.empty())
    {
      FileHandler().storeTransformations(trafo_out, trafo, {FileTypes::TRANSFORMATIONXML});
    }
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    // file list
    StringList file_list = getStringList_("in");

    // file type
    FileHandler file_handler;
    FileTypes::Type force_type;
    if (!getStringOption_("in_type").empty())
    {
      force_type = FileTypes::nameToType(getStringOption_("in_type"));
    }
    else
    {
      force_type = file_handler.getType(file_list[0]);
    }

    // output file names and types
    String out_file = getStringOption_("out");
   
    // append method
    bool append_rows = false;
    bool append_cols = false;
    String append_method = getStringOption_("append_method");
    append_method == "append_rows" ? append_rows = true : append_cols = true; 
   
    bool annotate_file_origin =  getFlag_("annotate_file_origin");
    rt_gap_ = getDoubleOption_("rt_concat:gap");
    vector<String> trafo_out = getStringList_("rt_concat:trafo_out");
    if (trafo_out.empty())
    {
      // resize now so we don't have to worry about indexing out of bounds:
      trafo_out.resize(file_list.size());
    }
    else if (trafo_out.size() != file_list.size())
    {
      writeLogError_("Error: Number of transformation output files must equal the number of input files (parameters 'rt_concat:trafo_out'/'in')!");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    if (force_type == FileTypes::FEATUREXML)
    {
      FeatureMap out;
      FileHandler fh;
      for (Size i = 0; i < file_list.size(); ++i)
      {
        FeatureMap map;
        fh.loadFeatures(file_list[i], map, {FileTypes::FEATUREXML});

        if (annotate_file_origin)
        {
          for (FeatureMap::iterator it = map.begin(); it != map.end(); ++it)
          {
            it->setMetaValue("file_origin", DataValue(file_list[i]));
          }
        }

        if (rt_gap_ > 0.0) // concatenate in RT
        {
          adjustRetentionTimes_(map, trafo_out[i], i == 0);
        }

        out += map;
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info
      addDataProcessing_(out, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

      fh.storeFeatures(out_file, out, {FileTypes::FEATUREXML});
    }

    else if (force_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap out;
      FileHandler fh;
      // load the metadata from the first file
      fh.loadConsensusFeatures(file_list[0], out, {FileTypes::CONSENSUSXML});
      // but annotate the origins

      if (append_rows) {
          if (annotate_file_origin)
          {
            for (ConsensusFeature& cm : out)
            {
              cm.setMetaValue("file_origin", DataValue(file_list[0]));
            }
          }

          // skip first file for adding
          for (Size i = 1; i < file_list.size(); ++i)
          {
            ConsensusMap map;
            fh.loadConsensusFeatures(file_list[i], map, {FileTypes::CONSENSUSXML});

            if (annotate_file_origin)
            {
              for (ConsensusFeature& cm : map)
              {
                cm.setMetaValue("file_origin", DataValue(file_list[i]));
              }  
            } 

            if (rt_gap_ > 0.0) // concatenate in RT
            {  
              adjustRetentionTimes_(map, trafo_out[i], i == 0);
            }

            out.appendRows(map);
          }
      }
      
      if (append_cols)
      { 
          // skip first file for adding
          for (Size i = 1; i < file_list.size(); ++i)
          {
            ConsensusMap map;
            fh.loadConsensusFeatures(file_list[i], map, {FileTypes::CONSENSUSXML});
            out.appendColumns(map);
          }
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info
      addDataProcessing_(out, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

      fh.storeConsensusFeatures(out_file, out,{FileTypes::CONSENSUSXML});
    }

    else if (force_type == FileTypes::FASTA)
    {
      FASTAFile infile;
      FASTAFile outfile;
      vector <FASTAFile::FASTAEntry> entries;
      vector <FASTAFile::FASTAEntry> temp_entries;
      vector <FASTAFile::FASTAEntry>::iterator loopiter;
      vector <FASTAFile::FASTAEntry>::iterator iter;

      for (Size i = 0; i < file_list.size(); ++i)
      {
        infile.load(file_list[i], temp_entries);
        entries.insert(entries.end(), temp_entries.begin(), temp_entries.end());
      }

      for (loopiter = entries.begin(); loopiter != entries.end(); loopiter = std::next(loopiter))
      {

        iter = find_if(entries.begin(), loopiter, [&loopiter](const FASTAFile::FASTAEntry& entry) { return entry.headerMatches(*loopiter); });

        if (iter != loopiter)
        {
          std::cout << "Warning: Duplicate header, Number: " << std::distance(entries.begin(), loopiter) + 1 << ", ID: " << loopiter->identifier << " is same as Number: " << std::distance(entries.begin(), iter) << ", ID: " << iter->identifier << "\n";
        }

        iter = find_if(entries.begin(), loopiter, [&loopiter](const FASTAFile::FASTAEntry& entry) { return entry.sequenceMatches(*loopiter); });

        if (iter != loopiter && iter != entries.end())
        {
          std::cout << "Warning: Duplicate sequence, Number: " << std::distance(entries.begin(), loopiter) + 1 << ", ID: " << loopiter->identifier << " is same as Number: " << std::distance(entries.begin(), iter) << ", ID: " << iter->identifier << "\n";
        }
      }

      outfile.store(out_file, entries);
    }

    else if (force_type == FileTypes::TRAML)
    {
      TargetedExperiment out;
      FileHandler fh;
      for (Size i = 0; i < file_list.size(); ++i)
      {
        TargetedExperiment map;
        fh.loadTransitions(file_list[i], map, {FileTypes::TRAML});
        out += map;
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info
      Software software;
      software.setName("FileMerger");
      software.setVersion(VersionInfo::getVersion());
      out.addSoftware(software);

      fh.storeTransitions(out_file, out, {FileTypes::TRAML});
    }
    else // raw data input (e.g. mzML)
    {
      // RT
      bool rt_auto_number = getFlag_("raw:rt_auto");
      bool rt_filename = getFlag_("raw:rt_filename");
      bool rt_custom = false;
      DoubleList custom_rts = getDoubleList_("raw:rt_custom");
      if (!custom_rts.empty())
      {
        rt_custom = true;
        if (custom_rts.size() != file_list.size())
        {
          writeLogError_("Error: Custom retention time list (parameter 'raw:rt_custom') must have as many elements as there are input files (parameter 'in')!");
          return ILLEGAL_PARAMETERS;
        }
      }

      // MS level
      Int ms_level = getIntOption_("raw:ms_level");

      PeakMap out;
      UInt rt_auto = 0;
      UInt native_id = 0;
      for (Size i = 0; i < file_list.size(); ++i)
      {
        String filename = file_list[i];

        // load file
        force_type = file_handler.getType(file_list[i]);
        PeakMap in;
        file_handler.loadExperiment(filename, in, {force_type}, log_type_, true, true);

        if (in.empty() && in.getChromatograms().empty())
        {
          writeLogWarn_(String("Warning: Empty file '") + filename + "'!");
          continue;
        }
        out.reserve(out.size() + in.size());

        // warn if custom RT and more than one scan in input file
        if (rt_custom && in.size() > 1)
        {
          writeLogWarn_(String("Warning: More than one scan in file '") + filename + "'! All scans will have the same retention time!");
        }

        // handle special raw data options:
        for (MSSpectrum& spec : in)
        {
          float rt_final = spec.getRT();
          if (rt_auto_number)
          {
            rt_final = ++rt_auto;
          }
          else if (rt_custom)
          {
            rt_final = custom_rts[i];
          }
          else if (rt_filename)
          {
            static const boost::regex re(R"(rt(\d+(\.\d+)?))");
            boost::smatch match;
            bool found = boost::regex_search(filename, match, re);
            if (found)
            {
              rt_final = String(match[1]).toFloat();
            }
            else
            {
              writeLogWarn_("Warning: could not extract retention time from filename '" + filename + "'");
            }
          }

          // none of the rt methods were successful
          if (rt_final < 0)
          {
            writeLogWarn_(String("Warning: No valid retention time for output scan '") + rt_auto + "' from file '" + filename + "'");
          }

          spec.setRT(rt_final);
          spec.setNativeID("spectrum=" + String(native_id));
          if (ms_level > 0)
          {
            spec.setMSLevel(ms_level);
          }
          ++native_id;
        }

        // if we have only one spectrum, we can annotate it directly, for more spectra, we just name the source file leaving the spectra unannotated (to avoid a long and redundant list of sourceFiles)
        if (in.size() == 1)
        {
          in[0].setSourceFile(in.getSourceFiles()[0]);
          in.getSourceFiles().clear(); // delete source file annotated from source file (it's in the spectrum anyways)
        }

        if (rt_gap_ > 0.0) // concatenate in RT
        {
          adjustRetentionTimes_(in, trafo_out[i], i == 0);
        }

        // add spectra to output
        for (const MSSpectrum& spec : in)
        {
          out.addSpectrum(spec);
        }
        // also add the chromatograms
        for (vector<MSChromatogram >::const_iterator
               chrom_it = in.getChromatograms().begin(); chrom_it != 
               in.getChromatograms().end(); ++chrom_it)
        {
          out.addChromatogram(*chrom_it);
        }

        // copy experimental settings from first file
        if (i == 0)
        {
          out.ExperimentalSettings::operator=(in);
        }
        else // otherwise append
        {
          out.getSourceFiles().insert(out.getSourceFiles().end(), in.getSourceFiles().begin(), in.getSourceFiles().end()); // could be emtpty if spectrum was annotated above, but that's ok then
        }
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      // annotate output with data processing info
      addDataProcessing_(out, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

      FileHandler().storeExperiment(out_file, out,{FileTypes::MZML}, log_type_);
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFileMerger tool;
  return tool.main(argc, argv);
}

/// @endcond
