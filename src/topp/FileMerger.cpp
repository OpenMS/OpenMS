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
// $Authors: Marc Sturm, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FileMerger FileMerger

  @brief Merges several files. Multiple output format supported, depending on input format.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileMerger \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool/instrument producing merge able files </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating merged files (e.g. @ref TOPP_XTandemAdapter) </td>
  </tr>
  </table>
  </center>

  The meta information that is valid for the whole experiment (e.g. MS instrument and sample)
  is taken from the first file.

  The retention times for the individual scans are taken from either:
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
    StringList valid_in = ListUtils::create<String>("mzData,mzXML,mzML,dta,dta2d,mgf,featureXML,consensusXML,fid,traML,FASTA");
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", valid_in);
    registerStringOption_("in_type", "<type>", "", "Input file type (default: determined from file extension or content)", false);
    setValidStrings_("in_type", valid_in);
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("mzML,featureXML,consensusXML,traML"));

    registerFlag_("annotate_file_origin", "Store the original filename in each feature using meta value \"file_origin\" (for featureXML and consensusXML only).");
    
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
      rt_offset_ = map.getMax()[0] + rt_gap_;
      trafo.fitModel("identity");
    }
    else // subsequent file -> apply transformation
    {
      TransformationDescription::DataPoints points(2);
      double rt_min = map.getMin()[0], rt_max = map.getMax()[0];
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
      TransformationXMLFile().store(trafo_out, trafo);
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
    if (getStringOption_("in_type").size() > 0)
    {
      force_type = FileTypes::nameToType(getStringOption_("in_type"));
    }
    else
    {
      force_type = file_handler.getType(file_list[0]);
    }

    // output file names and types
    String out_file = getStringOption_("out");

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
      writeLog_("Error: Number of transformation output files must equal the number of input files (parameters 'rt_concat:trafo_out'/'in')!");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    if (force_type == FileTypes::FEATUREXML)
    {
      FeatureMap out;
      FeatureXMLFile fh;
      for (Size i = 0; i < file_list.size(); ++i)
      {
        FeatureMap map;
        fh.load(file_list[i], map);

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

      fh.store(out_file, out);
    }

    else if (force_type == FileTypes::CONSENSUSXML)
    {
      ConsensusMap out;
      ConsensusXMLFile fh;
      // load the metadata from the first file
      fh.load(file_list[0], out);
      // but annotate the origins
      if (annotate_file_origin)
      {
        for (ConsensusMap::iterator it = out.begin(); it != out.end(); ++it)
        {
          it->setMetaValue("file_origin", DataValue(file_list[0]));
        }
      }

      // skip first file for adding
      for (Size i = 1; i < file_list.size(); ++i)
      {
        ConsensusMap map;
        fh.load(file_list[i], map);

        if (annotate_file_origin)
        {
          for (ConsensusMap::iterator it = map.begin(); it != map.end(); ++it)
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

      fh.store(out_file, out);
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
      TraMLFile fh;
      for (Size i = 0; i < file_list.size(); ++i)
      {
        TargetedExperiment map;
        fh.load(file_list[i], map);
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

      fh.store(out_file, out);
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
          writeLog_("Custom retention time list (parameter 'raw:rt_custom') must have as many elements as there are input files (parameter 'in')!");
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
        file_handler.loadExperiment(filename, in, force_type, log_type_);

        if (in.empty() && in.getChromatograms().empty())
        {
          writeLog_(String("Warning: Empty file '") + filename + "'!");
          continue;
        }
        out.reserve(out.size() + in.size());

        // warn if custom RT and more than one scan in input file
        if (rt_custom && in.size() > 1)
        {
          writeLog_(String("Warning: More than one scan in file '") + filename + "'! All scans will have the same retention time!");
        }

        // handle special raw data options:
        for (PeakMap::iterator spec_it = in.begin();
             spec_it != in.end(); ++spec_it)
        {
          float rt_final = spec_it->getRT();
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
            static const boost::regex re("rt(\\d+(\\.\\d+)?)");
            boost::smatch match;
            bool found = boost::regex_search(filename, match, re);
            if (found)
            {
              rt_final = String(match[1]).toFloat();
            }
            else
            {
              writeLog_("Warning: could not extract retention time from filename '" + filename + "'");
            }
          }

          // none of the rt methods were successful
          if (rt_final < 0)
          {
            writeLog_(String("Warning: No valid retention time for output scan '") + rt_auto + "' from file '" + filename + "'");
          }

          spec_it->setRT(rt_final);
          spec_it->setNativeID("spectrum=" + String(native_id));
          if (ms_level > 0)
          {
            spec_it->setMSLevel(ms_level);
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
        for (PeakMap::const_iterator spec_it = in.begin();
             spec_it != in.end(); ++spec_it)
        {
          out.addSpectrum(*spec_it);
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

      MzMLFile f;
      f.setLogType(log_type_);
      f.store(out_file, out);
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
