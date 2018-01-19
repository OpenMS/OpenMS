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
// $Maintainer: Lars Nilse $
// $Authors: Marc Sturm, Clemens Groepl, Lars Nilse $
// --------------------------------------------------------------------------

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/IndexedMzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FileInfo FileInfo
  @brief Shows basic information about the data in an OpenMS readable file.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FileInfo \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none ; console or text file</td>
  </tr>
  </table>
  </CENTER>

  This tool can show basic information about the data in several peak, feature and consensus feature files. It can
  - show information about the data range of a file (m/z, RT, intensity)
  - show a statistical summary for intensities, qualities, feature widths
  - show an overview of the metadata
  - validate several XML formats against their XML schema
  - check for corrupt data in a file (e.g., duplicate spectra)

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FileInfo.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_FileInfo.html

  In order to enrich the resulting data of your analysis pipeline or to quickly compare different outcomes of your pipeline you can invoke the aforementioned information of your input data and (intermediary) results.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
  //helper struct for identification data
  struct IdData
  {
    String identifier;
    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
  };

  /// Write SomeStatistics to a stream.
  template <class T>
  static ostream &operator<<(ostream &os, const Math::SummaryStatistics<T> &rhs)
  {
    return os << "  num. of values: " << rhs.count << "\n"
              << "  mean:           " << rhs.mean << "\n"
              << "  minimum:        " << rhs.min << "\n"
              << "  lower quartile: " << rhs.lowerq << "\n"
              << "  median:         " << rhs.median << "\n"
              << "  upper quartile: " << rhs.upperq << "\n"
              << "  maximum:        " << rhs.max << "\n"
              << "  variance:       " << rhs.variance << "\n";
  }
} // namespace

class TOPPFileInfo : public TOPPBase
{
public:
  TOPPFileInfo() : TOPPBase("FileInfo", "Shows basic information about the file, such as data ranges and file type.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzData,mzXML,mzML,dta,dta2d,mgf,featureXML,consensusXML,idXML,pepXML,fid,mzid,trafoXML,fasta"));
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content", false);
    setValidStrings_("in_type", ListUtils::create<String>("mzData,mzXML,mzML,dta,dta2d,mgf,featureXML,consensusXML,idXML,pepXML,fid,mzid,trafoXML"));
    registerOutputFile_("out", "<file>", "", "Optional output file. If left out, the output is written to the command line.", false);
    setValidFormats_("out", ListUtils::create<String>("txt"));
    registerOutputFile_("out_tsv", "<file>", "", "Second optional output file. Tab separated flat text file.", false, true);
    setValidFormats_("out_tsv", ListUtils::create<String>("csv"));
    registerFlag_("m", "Show meta information about the whole experiment");
    registerFlag_("p", "Shows data processing information");
    registerFlag_("s", "Computes a five-number statistics of intensities, qualities, and widths");
    registerFlag_("d", "Show detailed listing of all spectra and chromatograms (peak files only)");
    registerFlag_("c", "Check for corrupt data in the file (peak files only)");
    registerFlag_("v", "Validate the file only (for mzML, mzData, mzXML, featureXML, idXML, consensusXML, pepXML)");
    registerFlag_("i", "Check whether a given mzML file contains valid indices (conforming to the indexedmzML standard)");
  }

  template <class Map>
  void writeRangesHumanReadable_(const Map& map, ostream &os)
  {
    os << "Ranges:"
       << "\n"
       << "  retention time: " << String::number(map.getMin()[Peak2D::RT], 2) << " .. " << String::number(map.getMax()[Peak2D::RT], 2) << "\n"
       << "  mass-to-charge: " << String::number(map.getMin()[Peak2D::MZ], 2) << " .. " << String::number(map.getMax()[Peak2D::MZ], 2) << "\n"
       << "  intensity:      " << String::number(map.getMinInt(), 2) << " .. " << String::number(map.getMaxInt(), 2) << "\n"
       << "\n";
  }

  template <class Map>
  void writeRangesMachineReadable_(const Map& map, ostream &os)
  {
    os << "retention time (min)"
       << "\t" << String::number(map.getMin()[Peak2D::RT], 2) << "\n"
       << "retention time (max)"
       << "\t" << String::number(map.getMax()[Peak2D::RT], 2) << "\n"
       << "mass-to-charge (min)"
       << "\t" << String::number(map.getMin()[Peak2D::MZ], 2) << "\n"
       << "mass-to-charge (max)"
       << "\t" << String::number(map.getMax()[Peak2D::MZ], 2) << "\n"
       << "intensity (min)"
       << "\t" << String::number(map.getMinInt(), 2) << "\n"
       << "intensity (max)"
       << "\t" << String::number(map.getMaxInt(), 2) << "\n";
  }

  ExitCodes outputTo_(ostream &os, ostream &os_tsv)
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------

    // File names
    String in = getStringOption_("in");

    // File type
    FileHandler fh;
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

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

    os << "\n"
       << "-- General information --"
       << "\n"
       << "\n"
       << "File name: " << in << "\n"
       << "File type: " << FileTypes::typeToName(in_type) << "\n";

    os_tsv << "file name"
           << "\t" << in << "\n"
           << "file type"
           << "\t" << FileTypes::typeToName(in_type) << "\n";

    PeakMap exp;
    FeatureMap feat;
    ConsensusMap cons;
    IdData id_data;

    //-------------------------------------------------------------
    // Validation
    //-------------------------------------------------------------
    if (getFlag_("v"))
    {
      bool valid = true;
      os << "\n"
         << "Validating " << FileTypes::typeToName(in_type) << " file";
      switch (in_type)
      {
      case FileTypes::MZDATA:
        os << " against XML schema version " << MzDataFile().getVersion() << "\n";
        valid = MzDataFile().isValid(in, os);
        break;

      case FileTypes::MZML:
        os << " against XML schema version " << MzMLFile().getVersion() << "\n";
        valid = MzMLFile().isValid(in, os);
        break;

      case FileTypes::FEATUREXML:
        os << " against XML schema version " << FeatureXMLFile().getVersion() << "\n";
        valid = FeatureXMLFile().isValid(in, os);
        break;

      case FileTypes::IDXML:
        os << " against XML schema version " << IdXMLFile().getVersion() << "\n";
        valid = IdXMLFile().isValid(in, os);
        break;

      case FileTypes::MZIDENTML:
        os << " against XML schema version " << MzIdentMLFile().getVersion() << "\n";
        valid = MzIdentMLFile().isValid(in, os);
        break;

      case FileTypes::CONSENSUSXML:
        os << " against XML schema version " << ConsensusXMLFile().getVersion() << "\n";
        valid = ConsensusXMLFile().isValid(in, os);
        break;

      case FileTypes::MZXML:
        os << " against XML schema version " << MzXMLFile().getVersion() << "\n";
        valid = MzXMLFile().isValid(in, os);
        break;

      case FileTypes::PEPXML:
        os << " against XML schema version " << PepXMLFile().getVersion() << "\n";
        valid = PepXMLFile().isValid(in, os);
        break;

      case FileTypes::TRANSFORMATIONXML:
        os << " against XML schema version " << TransformationXMLFile().getVersion() << "\n";
        valid = TransformationXMLFile().isValid(in, os);
        break;

      default:
        os << "\n"
           << "Aborted: Validation of this file type is not supported!"
           << "\n";
        return EXECUTION_OK;
      }

      if (valid)
      {
        os << "Success - the file is valid!"
           << "\n";
      }
      else
      {
        os << "Failed - errors are listed above!"
           << "\n";
      }

      // semantic validation:
      if ((in_type == FileTypes::MZML) || (in_type == FileTypes::MZDATA))
      {
        if (!valid)
        {
          os << "\n"
             << "Semantic validation is not performed due to previous errors!"
             << "\n";
        }
        else
        {
          os << "\n"
             << "Semantically validating " << FileTypes::typeToName(in_type)
             << " file";
          if (in_type == FileTypes::MZDATA)
            os << " (EXPERIMENTAL)";
          os << ":"
             << "\n";

          StringList errors, warnings;
          if (in_type == FileTypes::MZML)
          {
            valid = MzMLFile().isSemanticallyValid(in, errors, warnings);
          }
          else
          {
            valid = MzDataFile().isSemanticallyValid(in, errors, warnings);
          }

          for (Size i = 0; i < warnings.size(); ++i)
          {
            os << "Warning: " << warnings[i] << "\n";
          }
          for (Size i = 0; i < errors.size(); ++i)
          {
            os << "Error: " << errors[i] << "\n";
          }
          if (valid)
          {
            os << "Success - the file is semantically valid!"
               << "\n";
          }
          else
          {
            os << "Failed - errors are listed above!"
               << "\n";
          }
        }
      }

      return EXECUTION_OK;
    }

    //-------------------------------------------------------------
    // Validation of indices
    //-------------------------------------------------------------
    if (getFlag_("i"))
    {
      if (in_type != FileTypes::MZML)
      {
        writeLog_("Error: Can only validate indices for mzML files");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      std::cout << "Checking mzML file for valid indices ... " << std::endl;
      IndexedMzMLFile ifile;
      ifile.openFile(in);
      if (ifile.getParsingSuccess())
      {
        // Validate that we can access each single spectrum and chromatogram
        for (int i = 0; i < (int)ifile.getNrSpectra(); i++)
        {
          OpenMS::Interfaces::SpectrumPtr p = ifile.getSpectrumById(i);
        }
        for (int i = 0; i < (int)ifile.getNrChromatograms(); i++)
        {
          OpenMS::Interfaces::ChromatogramPtr p = ifile.getChromatogramById(i);
        }

        std::cout << "Found a valid indexed mzML XML File with " << ifile.getNrSpectra() << " spectra and " << ifile.getNrChromatograms() << " chromatograms." << std::endl
                  << std::endl;
      }
      else
      {
        std::cout << "Could not detect a valid index for the mzML file " << in << "\nEither the index is not present or is not correct." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    //-------------------------------------------------------------
    // Content statistics
    //-------------------------------------------------------------
    Map<String, int> meta_names;

    if (in_type == FileTypes::FASTA)
    {
      vector<FASTAFile::FASTAEntry> entries;
      FASTAFile file;

      map<char, int> aacids;
      size_t number_of_aacids = 0;

      SysInfo::MemUsage mu;
      //loading input
      file.load(in, entries);
      std::cout << "\n\n" << mu.delta("loading FASTA") << std::endl;

      os << "Number of sequences: " << entries.size() << "\n"
         << "\n";

      for (auto loopiter = entries.begin(); loopiter != entries.end(); ++loopiter)
      {
        auto iter = find_if(entries.begin(), loopiter, [&loopiter](const FASTAFile::FASTAEntry& entry) { return entry.headerMatches(*loopiter); });
        if (iter != loopiter)
        {
          os << "Warning: Duplicate header, Number: " << std::distance(entries.begin(), loopiter) << ", ID: " << loopiter->identifier << " is same as Number: " << std::distance(entries.begin(), iter) << ", ID: " << iter->identifier << "\n";
        }

        iter = find_if(entries.begin(), loopiter, [&loopiter](const FASTAFile::FASTAEntry& entry) { return entry.sequenceMatches(*loopiter); });
        if (iter != loopiter && iter != entries.end())
        {
          os << "Warning: Duplicate sequence, Number: " << std::distance(entries.begin(), loopiter) << ", ID: " << loopiter->identifier << " is same as Number: " << std::distance(entries.begin(), iter) << ", ID: " << iter->identifier << "\n";
        }

        for (Size i = 0; i < loopiter->sequence.size(); i++)
        {
          ++aacids[loopiter->sequence[i]];
        }
        number_of_aacids += loopiter->sequence.size();
      }

      os << "Total amino acids: " << number_of_aacids << "\n\n";
      os << "Amino acid counts: "
         << "\n";

      for (auto it = aacids.begin(); it != aacids.end(); ++it)
      {
        os << it->first << '\t' << it->second << "\n";
      }
    }

    else if (in_type == FileTypes::FEATUREXML) //features
    {
      FeatureXMLFile ff;
      ff.getOptions().setLoadConvexHull(false);   // CH's currently not needed here
      ff.getOptions().setLoadSubordinates(false); // SO's currently not needed here

      SysInfo::MemUsage mu;
      // reading input
      ff.load(in, feat);
      std::cout << "\n\n" << mu.delta("loading featureXML") << std::endl;

      feat.updateRanges();

      os << "Number of features: " << feat.size() << "\n"
         << "\n";
      writeRangesHumanReadable_(feat, os);
      writeRangesMachineReadable_(feat, os_tsv);

      // Charge distribution and TIC
      Map<Int, UInt> charges;
      Map<size_t, UInt> numberofids;
      double tic = 0.0;
      for (Size i = 0; i < feat.size(); ++i)
      {
        ++charges[feat[i].getCharge()];
        tic += feat[i].getIntensity();
        const vector<PeptideIdentification> &peptide_ids = feat[i].getPeptideIdentifications();
        ++numberofids[peptide_ids.size()];
      }

      os << "Total ion current in features: " << tic << "\n";
      os << "\n"
         << "Charge distribution:"
         << "\n";
      for (auto it = charges.begin(); it != charges.end(); ++it)
      {
        os << "  charge " << it->first << ": " << it->second << "\n";
      }

      os << "\n"
         << "Distribution of peptide identifications (IDs) per feature:\n";
      for (auto it = numberofids.begin(); it != numberofids.end(); ++it)
      {
        os << "  " << it->first << " IDs: " << it->second << "\n";
      }

      os << "\n"
         << "Unassigned peptide identifications: " << feat.getUnassignedPeptideIdentifications().size() << "\n";
    }
    else if (in_type == FileTypes::CONSENSUSXML) //consensus features
    {

      SysInfo::MemUsage mu;
      // reading input
      ConsensusXMLFile().load(in, cons);
      std::cout << "\n\n" << mu.delta("loading consensusXML") << std::endl;

      cons.updateRanges();

      map<Size, UInt> num_consfeat_of_size;
      for (ConsensusMap::const_iterator cmit = cons.begin(); cmit != cons.end(); ++cmit)
      {
        ++num_consfeat_of_size[cmit->size()];
      }
      if (num_consfeat_of_size.empty())
      {
        os << "\n"
           << "Number of consensus features: 0"
           << "\n";
        os << "No consensus features found, map is empty!"
           << "\n\n";
      }
      else
      {
        Size field_width = num_consfeat_of_size.rbegin()->first / 10 + 1;
        os << "\n"
           << "Number of consensus features:"
           << "\n";
        for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin(); i != num_consfeat_of_size.rend(); ++i)
        {
          os << "  of size " << setw(field_width) << i->first << ": " << i->second << "\n";
        }
        os << "  total:    " << string(field_width, ' ') << cons.size() << "\n"
           << "\n";

        writeRangesHumanReadable_(cons, os);
        writeRangesMachineReadable_(cons, os_tsv);
      }

      // file descriptions
      const ConsensusMap::FileDescriptions &descs = cons.getFileDescriptions();
      if (!descs.empty())
      {
        os << "File descriptions:"
           << "\n";
        for (ConsensusMap::FileDescriptions::const_iterator it = descs.begin(); it != descs.end(); ++it)
        {
          os << "  " << it->second.filename << ":"
             << "\n"
             << "    identifier: " << it->first << "\n"
             << "    label:      " << it->second.label << "\n"
             << "    size:       " << it->second.size << "\n";
        }
        os << "\n";
      }
    }
    else if (in_type == FileTypes::IDXML || in_type == FileTypes::MZIDENTML) //identifications
    {
      UInt spectrum_count(0);
      Size peptide_hit_count(0);
      UInt runs_count(0);
      Size protein_hit_count(0);
      set<String> peptides;
      set<String> proteins;
      Size modified_peptide_count(0);
      Map<String, int> mod_counts;
      vector<uint16_t> peptide_length;

      // reading input
      SysInfo::MemUsage mu;
      if (in_type == FileTypes::MZIDENTML)
      {
        MzIdentMLFile().load(in, id_data.proteins, id_data.peptides);
      }
      else
      {
        IdXMLFile().load(in, id_data.proteins, id_data.peptides, id_data.identifier);
      }
      std::cout << "\n\n" << mu.delta("loading idXML") << std::endl;

      // export metadata to second output stream
      os_tsv << "database"
             << "\t" << id_data.proteins[0].getSearchParameters().db << "\n"
             << "database version"
             << "\t" << id_data.proteins[0].getSearchParameters().db_version << "\n"
             << "taxonomy"
             << "\t" << id_data.proteins[0].getSearchParameters().taxonomy << "\n";

      // calculations
      for (Size i = 0; i < id_data.peptides.size(); ++i)
      {
        if (!id_data.peptides[i].empty())
        {
          ++spectrum_count;
          peptide_hit_count += id_data.peptides[i].getHits().size();
          const vector<PeptideHit> &temp_hits = id_data.peptides[i].getHits();
          // collect stats about modifications from TOP HIT!
          if (temp_hits[0].getSequence().isModified())
          {
            ++modified_peptide_count;
            AASequence aa = temp_hits[0].getSequence();
            for (Size ia = 0; ia < aa.size(); ++ia)
            {
              if (aa[ia].isModified())
                ++mod_counts[aa[ia].getModificationName()];
            }
          }
          for (Size j = 0; j < temp_hits.size(); ++j)
          {
            peptides.insert(temp_hits[j].getSequence().toString());
            peptide_length.push_back((uint16_t)temp_hits[j].getSequence().size());
          }
        }
      }
      for (Size i = 0; i < id_data.proteins.size(); ++i)
      {
        ++runs_count;
        protein_hit_count += id_data.proteins[i].getHits().size();
        const vector<ProteinHit> &temp_hits = id_data.proteins[i].getHits();
        for (Size j = 0; j < temp_hits.size(); ++j)
        {
          proteins.insert(temp_hits[j].getAccession());
        }
      }
      if (peptide_length.empty())
      { // avoid invalid-range exception when computing mean()
        peptide_length.push_back(0); 
      }

      os << "Number of:"
         << "\n";
      os << "  runs:                       " << runs_count << "\n";
      os << "  protein hits:               " << protein_hit_count << "\n";
      os << "  non-redundant protein hits: " << proteins.size() << "\n";
      os << "  (only hits that differ in the accession)"
         << "\n";
      os << "\n";
      os << "  matched spectra:    " << spectrum_count << "\n";
      os << "  peptide hits:               " << peptide_hit_count << " (avg. length: " << Math::round(Math::mean(peptide_length.begin(), peptide_length.end())) << ")\n";
      os << "  modified top-hits:          " << modified_peptide_count << "/" << spectrum_count << (spectrum_count > 0 ? String(" (") + Math::round(modified_peptide_count * 1000.0 / spectrum_count) / 10 + "%)" : "") << "\n";
      os << "  non-redundant peptide hits: " << peptides.size() << "\n";
      os << "  (only hits that differ in sequence and/or modifications)"
         << "\n";
      for (Map<String, int>::ConstIterator it = mod_counts.begin(); it != mod_counts.end(); ++it)
      {
        if (it != mod_counts.begin())
          os << ", ";
        else
          os << "  Modifications (top-hits only): ";
        os << it->first << "(" << it->second << ")";
      }

      os_tsv << "peptide hits"
             << "\t" << peptide_hit_count << "\n";
      os_tsv << "non-redundant peptide hits (only hits that differ in sequence and/or modifications): "
             << "\t" << peptides.size() << "\n";
      os_tsv << "protein hits"
             << "\t" << protein_hit_count << "\n";
      os_tsv << "non-redundant protein hits (only hits that differ in the accession)"
             << "\t" << proteins.size() << "\n";
    }
    else if (in_type == FileTypes::PEPXML)
    {
      os << "\nFor pepXML files, only validation against the XML schema is implemented at this point."
         << "\n";
    }
    else if (in_type == FileTypes::TRANSFORMATIONXML)
    {
      TransformationDescription trafo;
      TransformationXMLFile().load(in, trafo);
      os << "\nTransformation model: " << trafo.getModelType() << "\n";
      trafo.printSummary(os);
    }
    else // peaks
    {

      SysInfo::MemUsage mu;
      if (!fh.loadExperiment(in, exp, in_type, log_type_, false, false))
      {
        writeLog_("Unsupported or corrupt input file. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      // update range information and retrieve which MS levels were recorded
      exp.updateRanges();
      vector<UInt> levels = exp.getMSLevels();

      // report memory consumption
      std::cout << "\n\n" << mu.delta("loading MS data") << std::endl;

      os << "\n";

      // check if the meta data indicates that this is peak data
      // and count how many spectra per MS level there are
      map<Size, UInt> level_annotated_picked;
      map<Size, UInt> level_estimated_picked;
      map<Size, UInt> counts;
      for (Size i = 0; i != exp.size(); ++i) 
      { 
        // read stored metadata
        auto peak_type = exp[i].getType();


        const Size level = exp[i].getMSLevel();

        ++counts[level];  // count MS level

        // annotate peak type (profile / centroided) from meta data
        if (level_annotated_picked.count(level) == 0)
        {
          // Some conversion software only annotate "MS:1000525 spectrum representation" leading to an UNKNOWN type
          // Fortunately, some store a data processing item that indicates that the data has been picked
          if (peak_type == SpectrumSettings::UNKNOWN)
          {
            for (auto dp : exp[i].getDataProcessing())
            {
              if (dp->getProcessingActions().count(DataProcessing::PEAK_PICKING) == 1)
              {
                peak_type = SpectrumSettings::CENTROID;
              }
            }
          }
          level_annotated_picked[level] = peak_type;
        }
        
        // estimate peak type once for every level (take a spectrum with enough peaks for stable estimation) 
        if (level_estimated_picked.count(level) == 0 && exp[i].size() > 10)
        {
          level_estimated_picked[level] = PeakTypeEstimator().estimateType(exp[i].begin(), exp[i].end());
        }
      }

      os << "MS levels: ";
      if (!levels.empty())
      {
        os << *(levels.begin());
        for (vector<UInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
        {
          os << ", " << *it;
        }
      }
      os << "\n";

      // basic info
      os << "Number of spectra: " << exp.size() << "\n";
      os << "Number of peaks: " << exp.getSize() << "\n"
         << "\n";
      os_tsv << "number of spectra"
             << "\t" << exp.size() << "\n"
             << "number of peaks"
             << "\t" << exp.getSize() << "\n";

      // output
      if (!counts.empty())
      {
        os << "Number of spectra per MS level:"
           << "\n";
        for (auto it = counts.begin(); it != counts.end(); ++it)
        {
          os << "  level " << it->first << ": " << it->second << "\n";
          os_tsv << "number of MS" << it->first << " spectra"
                 << "\t" << it->second << "\n";
        }
        os << "\n";
      }

      writeRangesHumanReadable_(exp, os);
      writeRangesMachineReadable_(exp, os_tsv);

      // write peak types (centroided / profile mode)
      os << "Peak type metadata (estimated)\n"; 
      for (auto const l : levels)
      {
        os << "  level " << l << ": " 
           << SpectrumSettings::NamesOfSpectrumType[level_annotated_picked[l]] << " ("
           << SpectrumSettings::NamesOfSpectrumType[level_estimated_picked[l]] << ")\n";
      }

      // show meta data array names
      for (PeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
      {
        for (Size i = 0; i < it->getFloatDataArrays().size(); ++i)
        {
          String name = it->getFloatDataArrays()[i].getName();
          if (meta_names.has(name))
          {
            meta_names[name]++;
          }
          else
          {
            meta_names[name] = 1;
          }
        }
        for (Size i = 0; i < it->getIntegerDataArrays().size(); ++i)
        {
          String name = it->getIntegerDataArrays()[i].getName();
          if (meta_names.has(name))
          {
            meta_names[name]++;
          }
          else
          {
            meta_names[name] = 1;
          }
        }
        for (Size i = 0; i < it->getStringDataArrays().size(); ++i)
        {
          String name = it->getStringDataArrays()[i].getName();
          if (meta_names.has(name))
          {
            meta_names[name]++;
          }
          else
          {
            meta_names[name] = 1;
          }
        }
      }
      if (!meta_names.empty())
      {
        // nice formatting:
        Size max_length = 0;
        for (Map<String, int>::ConstIterator it = meta_names.begin(); it != meta_names.end(); ++it)
        {
          if (it->first.size() > max_length)
            max_length = it->first.size();
        }
        os << "Meta data array:"
           << "\n";
        for (Map<String, int>::ConstIterator it = meta_names.begin(); it != meta_names.end(); ++it)
        {
          String padding(max_length - it->first.size(), ' ');
          os << "  " << it->first << ": " << padding << it->second << " spectra"
             << "\n";
        }
        os << "\n";
      }

      // some chromatogram information
      if (!exp.getChromatograms().empty())
      {
        os << "Number of chromatograms: " << exp.getChromatograms().size() << "\n";
        os_tsv << "number of chromatograms"
               << "\t" << exp.getChromatograms().size() << "\n";

        Size num_chrom_peaks(0);
        Map<ChromatogramSettings::ChromatogramType, Size> chrom_types;
        for (vector<MSChromatogram>::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
        {
          num_chrom_peaks += it->size();
          if (chrom_types.has(it->getChromatogramType()))
          {
            chrom_types[it->getChromatogramType()]++;
          }
          else
          {
            chrom_types[it->getChromatogramType()] = 1;
          }
        }
        os << "Number of chromatographic peaks: " << num_chrom_peaks << "\n"
           << "\n";
        os_tsv << "number of chromatographic peaks"
               << "\t" << num_chrom_peaks << "\n";

        os << "Number of chromatograms per type: "
           << "\n";
        for (Map<ChromatogramSettings::ChromatogramType, Size>::const_iterator it = chrom_types.begin(); it != chrom_types.end(); ++it)
        {
          os << String("  ") + ChromatogramSettings::ChromatogramNames[it->first] + ":                         "
             << it->second << "\n";
        }
        if (getFlag_("d") && chrom_types.has(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM))
        {
          os << "\n"
             << " -- Detailed chromatogram listing -- "
             << "\n";
          os << "\nSelected Reaction Monitoring Transitions:"
             << "\n";
          os << "Q1 Q3 RT_begin RT_end name comment"
             << "\n";
          for (vector<MSChromatogram>::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
          {
            if (it->getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
            {
              os << it->getPrecursor().getMZ() << " " << it->getProduct().getMZ() << " " << it->front().getRT() << " " << it->back().getRT() << " " << it->getName() << " " << it->getComment() << "\n";
            }
          }
        }
      }

      // Detailed listing of scans
      if (getFlag_("d") && exp.size() > 0)
      {
        os << "\n"
           << "-- Detailed spectrum listing --"
           << "\n";
        UInt count = 0;
        for (auto const& spectrum : exp)
        {
          ++count;
          os << "\n"
             << "Spectrum " << count << ":"
             << "\n"
             << "  mslevel:    " << spectrum.getMSLevel() << "\n"
             << "  scanMode:   " << InstrumentSettings::NamesOfScanMode[spectrum.getInstrumentSettings().getScanMode()] << "\n"
             << "  peaks:      " << spectrum.size() << "\n"
             << "  RT:         " << spectrum.getRT() << "\n"
             << "  m/z:        ";

          if (!spectrum.empty())
          {
            os << spectrum.begin()->getMZ() << " .. " << spectrum.rbegin()->getMZ() << "\n";
          }

          os << "Precursors:  " << spectrum.getPrecursors().size() <<  "\n";

          auto pc_count = UInt{0};
          for (auto const& pc : spectrum.getPrecursors())
          {
            os << "Precursor[" << pc_count << "]\n"
               << "  charge: " << pc.getCharge() << "\n"
               << "  mz:     " << pc.getMZ() << "\n"
               << "  activation methods: \n";
            for (auto const& am : pc.getActivationMethods())
            {
              os << "    " << Precursor::NamesOfActivationMethodShort[am] << " (" << Precursor::NamesOfActivationMethod[am] << ")\n";
            }

            os << "\n";

            ++pc_count;
          }
        }
      }

      // Check for corrupt data
      if (getFlag_("c"))
      {
        os << "\n"
           << "-- Checking for corrupt data --"
           << "\n"
           << "\n";
        // RTs sorted?
        if (!exp.isSorted(false))
        {
          os << "Error: Spectrum retention times are not sorted in ascending order"
             << "\n";
        }
        vector<double> ms1_rts;
        ms1_rts.reserve(exp.size());
        for (Size s = 0; s < exp.size(); ++s)
        {
          // ms level = 0
          if (exp[s].getMSLevel() == 0)
          {
            os << "Error: MS-level 0 in spectrum (RT: " << exp[s].getRT() << ")"
               << "\n";
          }
          //scan size = 0
          if (exp[s].empty())
          {
            os << "Warning: No peaks in spectrum (RT: " << exp[s].getRT() << ")"
               << "\n";
          }
          //duplicate meta data array names
          Map<String, int> names;
          for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
          {
            String name = exp[s].getFloatDataArrays()[m].getName();
            if (names.has(name))
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << "\n";
            }
            else
            {
              names[name] = 0;
            }
          }
          for (Size m = 0; m < exp[s].getIntegerDataArrays().size(); ++m)
          {
            String name = exp[s].getIntegerDataArrays()[m].getName();
            if (names.has(name))
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << "\n";
            }
            else
            {
              names[name] = 0;
            }
          }
          for (Size m = 0; m < exp[s].getStringDataArrays().size(); ++m)
          {
            String name = exp[s].getStringDataArrays()[m].getName();
            if (names.has(name))
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << "\n";
            }
            else
            {
              names[name] = 0;
            }
          }
          //duplicate scans (part 1)
          if (exp[s].getMSLevel() == 1)
          {
            ms1_rts.push_back(exp[s].getRT());
          }
        }
        //duplicate scans (part 2)
        sort(ms1_rts.begin(), ms1_rts.end());
        for (Size i = 1; i < ms1_rts.size(); ++i)
        {
          if (ms1_rts[i - 1] == ms1_rts[i])
            os << "Error: Duplicate spectrum retention time: " << ms1_rts[i] << "\n";
        }
        //check peaks
        for (Size s = 0; s < exp.size(); ++s)
        {
          //peaks sorted?
          if (!exp[s].isSorted())
          {
            os << "Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: " << exp[s].getRT() << ")"
               << "\n";
          }
          vector<double> mzs;
          mzs.reserve(exp[s].size());
          for (Size p = 0; p < exp[s].size(); ++p)
          {
            //negative intensity
            if (exp[s][p].getIntensity() < 0.0)
            {
              os << "Warning: Negative peak intensity peak (RT: " << exp[s].getRT() << " MZ: " << exp[s][p].getMZ() << " intensity: " << exp[s][p].getIntensity() << ")"
                 << "\n";
            }
            //duplicate m/z (part 1)
            mzs.push_back(exp[s][p].getMZ());
          }
          //duplicate m/z (part 2)
          sort(mzs.begin(), mzs.end());
          for (Size i = 1; i < mzs.size(); ++i)
          {
            if (mzs[i - 1] == mzs[i])
              os << "Error: Duplicate peak m/z " << mzs[i] << " in spectrum (RT: " << exp[s].getRT() << ")"
                 << "\n";
          }
        }
      }
    }

    //-------------------------------------------------------------
    // meta information
    //-------------------------------------------------------------
    if (getFlag_("m") || getStringOption_("out_tsv") != "")
    {

      //basic info
      os << "\n"
         << "-- Meta information --"
         << "\n"
         << "\n";

      if (in_type == FileTypes::FEATUREXML) //features
      {
        os << "Document ID: " << feat.getIdentifier() << "\n"
           << "\n";
      }
      else if (in_type == FileTypes::CONSENSUSXML) //consensus features
      {
        os << "Document ID: " << cons.getIdentifier() << "\n"
           << "\n";
      }
      else if (in_type == FileTypes::IDXML) //identifications
      {
        os << "Document ID: " << id_data.identifier << "\n"
           << "\n";
      }
      else if (in_type == FileTypes::PEPXML)
      {
        // TODO
      }
      else if (in_type == FileTypes::FASTA)
      {
      }
      else //peaks
      {

        os << "Document ID:        " << exp.getIdentifier() << "\n"
           << "Date:               " << exp.getDateTime().get() << "\n";
        os_tsv << "document id"
               << "\t" << exp.getIdentifier() << "\n"
               << "date"
               << "\t" << exp.getDateTime().get() << "\n";

        //basic info
        os << "\n"
           << "Sample:"
           << "\n"
           << "  name:             " << exp.getSample().getName() << "\n"
           << "  organism:         " << exp.getSample().getOrganism() << "\n"
           << "  comment:          " << exp.getSample().getComment() << "\n";
        os_tsv << "sample name"
               << "\t" << exp.getSample().getName() << "\n"
               << "sample organism"
               << "\t" << exp.getSample().getOrganism() << "\n"
               << "sample comment"
               << "\t" << exp.getSample().getComment() << "\n";

        //instrument info
        os << "\n"
           << "Instrument:"
           << "\n"
           << "  name:             " << exp.getInstrument().getName() << "\n"
           << "  model:            " << exp.getInstrument().getModel() << "\n"
           << "  vendor:           " << exp.getInstrument().getVendor() << "\n"
           << "  ion source(s):    ";
        os_tsv << "instrument name"
               << "\t" << exp.getInstrument().getName() << "\n"
               << "instrument model"
               << "\t" << exp.getInstrument().getModel() << "\n"
               << "instrument vendor"
               << "\t" << exp.getInstrument().getVendor() << "\n";
        for (Size i = 0; i < exp.getInstrument().getIonSources().size(); ++i)
        {
          os << IonSource::NamesOfIonizationMethod[exp.getInstrument().getIonSources()[i].getIonizationMethod()];
          if (i != exp.getInstrument().getIonSources().size() - 1)
          {
            os << ", ";
          }
        }
        os << "\n"
           << "  mass analyzer(s): ";
        for (Size i = 0; i < exp.getInstrument().getMassAnalyzers().size(); ++i)
        {
          os << MassAnalyzer::NamesOfAnalyzerType[exp.getInstrument().getMassAnalyzers()[i].getType()];
          if (i != exp.getInstrument().getMassAnalyzers().size() - 1)
          {
            os << ", ";
          }
        }
        os << "\n"
           << "  detector(s):      ";
        for (Size i = 0; i < exp.getInstrument().getIonDetectors().size(); ++i)
        {
          os << IonDetector::NamesOfType[exp.getInstrument().getIonDetectors()[i].getType()];
          if (i != exp.getInstrument().getIonDetectors().size() - 1)
            os << ", ";
        }
        os << "\n"
           << "\n";

        //contact persons
        for (Size i = 0; i < exp.getContacts().size(); ++i)
        {
          os << "Contact person:"
             << "\n"
             << "  first name:     " << exp.getContacts()[i].getFirstName() << "\n"
             << "  last name:      " << exp.getContacts()[i].getLastName() << "\n"
             << "  email:          " << exp.getContacts()[i].getEmail() << "\n"
             << "\n";
        }
      }
    }

    //-------------------------------------------------------------
    // data processing
    //-------------------------------------------------------------
    if (getFlag_("p"))
    {
      //basic info
      os << "\n"
         << "-- Data processing information --"
         << "\n"
         << "\n";

      //get data processing info
      vector<DataProcessing> dp;
      if (in_type == FileTypes::FEATUREXML) //features
      {
        dp = feat.getDataProcessing();
      }
      else if (in_type == FileTypes::CONSENSUSXML) //consensus features
      {
        dp = cons.getDataProcessing();
      }
      else if (in_type == FileTypes::IDXML) //identifications
      {
      }
      else if (in_type == FileTypes::PEPXML)
      {
        // TODO
      }
      else if (in_type == FileTypes::FASTA)
      {
      }
      else //peaks
      {
        if (!exp.empty())
        {
          os << "Note: The data is taken from the first spectrum!"
             << "\n"
             << "\n";
          for (Size i = 0; i < exp[0].getDataProcessing().size(); i++)
          {
            dp.push_back(*exp[0].getDataProcessing()[i].get());
          }
        }
      }

      //print data
      if (dp.empty())
      {
        os << "No information about data processing available!"
           << "\n"
           << "\n";
      }
      else
      {
        for (Size i = 0; i < dp.size(); ++i)
        {
          os << "Processing " << (i + 1) << ":"
             << "\n";
          os << "  software name:    " << dp[i].getSoftware().getName() << "\n";
          os << "  software version: " << dp[i].getSoftware().getVersion() << "\n";
          os << "  completion time:  " << dp[i].getCompletionTime().get() << "\n";
          os << "  actions:          ";
          for (set<DataProcessing::ProcessingAction>::const_iterator it = dp[i].getProcessingActions().begin();
               it != dp[i].getProcessingActions().end(); ++it)
          {
            if (it != dp[i].getProcessingActions().begin())
              os << ", ";
            os << DataProcessing::NamesOfProcessingAction[*it];
          }
          os << "\n"
             << "\n";
        }
      }
    }

    //-------------------------------------------------------------
    // statistics
    //-------------------------------------------------------------
    if (getFlag_("s"))
    {
      os << "\n"
         << "-- Statistics --"
         << "\n"
         << "\n";

      if (in_type == FileTypes::FEATUREXML) //features
      {
        Size size = feat.size();

        vector<double> intensities(size);
        vector<double> overall_qualities(size);
        vector<double> mz_qualities(size);
        vector<double> rt_qualities(size);
        vector<double> peak_widths(size);

        Size idx = 0;
        for (FeatureMap::const_iterator fm_iter = feat.begin();
             fm_iter != feat.end(); ++fm_iter, ++idx)
        {
          intensities[idx] = fm_iter->getIntensity();
          overall_qualities[idx] = fm_iter->getOverallQuality();
          rt_qualities[idx] = fm_iter->getQuality(Feature::RT);
          mz_qualities[idx] = fm_iter->getQuality(Feature::MZ);
          peak_widths[idx] = fm_iter->getWidth();
        }

        os.precision(writtenDigits<>(Feature::IntensityType()));
        os << "Intensities:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(intensities) << "\n";

        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Feature FWHM in RT dimension:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(peak_widths) << "\n";

        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Overall qualities:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(overall_qualities) << "\n";

        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Qualities in retention time dimension:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(rt_qualities) << "\n";

        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Qualities in mass-to-charge dimension:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(mz_qualities) << "\n";
      }
      else if (in_type == FileTypes::CONSENSUSXML) //consensus features
      {
        Size size = cons.size();

        vector<double> intensities;
        intensities.reserve(size);
        vector<double> qualities(size);
        qualities.reserve(size);
        vector<double> widths(size);
        widths.reserve(size);

        vector<double> rt_delta_by_elems;
        vector<double> rt_aad_by_elems;
        vector<double> rt_aad_by_cfs;
        rt_aad_by_cfs.reserve(size);

        vector<double> mz_delta_by_elems;
        vector<double> mz_aad_by_elems;
        vector<double> mz_aad_by_cfs;
        mz_aad_by_cfs.reserve(size);

        vector<double> it_delta_by_elems;
        vector<double> it_aad_by_elems;
        vector<double> it_aad_by_cfs;
        it_aad_by_cfs.reserve(size);

        for (ConsensusMap::const_iterator cm_iter = cons.begin();
             cm_iter != cons.end(); ++cm_iter)
        {
          double rt_aad = 0;
          double mz_aad = 0;
          double it_aad = 0;
          intensities.push_back(cm_iter->getIntensity());
          qualities.push_back(cm_iter->getQuality());
          widths.push_back(cm_iter->getWidth());
          for (ConsensusFeature::HandleSetType::const_iterator hs_iter = cm_iter->begin();
               hs_iter != cm_iter->end(); ++hs_iter)
          {
            double rt_diff = hs_iter->getRT() - cm_iter->getRT();
            rt_delta_by_elems.push_back(rt_diff);
            if (rt_diff < 0)
            {
              rt_diff = -rt_diff;
            }
            rt_aad_by_elems.push_back(rt_diff);
            rt_aad += rt_diff;
            double mz_diff = hs_iter->getMZ() - cm_iter->getMZ();
            mz_delta_by_elems.push_back(mz_diff);
            if (mz_diff < 0)
            {
              mz_diff = -mz_diff;
            }
            mz_aad_by_elems.push_back(mz_diff);
            mz_aad += mz_diff;
            double it_ratio = hs_iter->getIntensity() / (cm_iter->getIntensity() ? cm_iter->getIntensity() : 1.);
            it_delta_by_elems.push_back(it_ratio);
            if (it_ratio < 1.)
            {
              it_ratio = 1. / it_ratio;
            }
            it_aad_by_elems.push_back(it_ratio);
            it_aad += it_ratio;
          }
          if (!cm_iter->empty())
          {
            rt_aad /= cm_iter->size();
            mz_aad /= cm_iter->size();
            it_aad /= cm_iter->size();
          } // otherwise rt_aad etc. are 0 anyway
          rt_aad_by_cfs.push_back(rt_aad);
          mz_aad_by_cfs.push_back(mz_aad);
          it_aad_by_cfs.push_back(it_aad);
        }

        os.precision(writtenDigits(ConsensusFeature::IntensityType()));
        os << "Intensities of consensus features:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(intensities) << "\n";

        os.precision(writtenDigits(ConsensusFeature::QualityType()));
        os << "Qualities of consensus features:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(qualities) << "\n";

        os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
        os << "Retention time differences (\"element - center\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(rt_delta_by_elems) << "\n";
        os << "Absolute retention time differences (\"|element - center|\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(rt_aad_by_elems) << "\n";
        os << "Average absolute differences of retention time within consensus features (\"|element - center|\", weight 1 per consensus features):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(rt_aad_by_cfs) << "\n";

        os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
        os << "Mass-to-charge differences (\"element - center\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(mz_delta_by_elems) << "\n";
        os << "Absolute differences of mass-to-charge (\"|element - center|\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(mz_aad_by_elems) << "\n";
        os << "Average absolute differences of mass-to-charge within consensus features (\"|element - center|\", weight 1 per consensus features):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(mz_aad_by_cfs) << "\n";

        os.precision(writtenDigits(ConsensusFeature::IntensityType()));
        os << "Intensity ratios (\"element / center\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(it_delta_by_elems) << "\n";
        os << "Relative intensity error (\"max{(element / center), (center / element)}\", weight 1 per element):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(it_aad_by_elems) << "\n";
        os << "Average relative intensity error within consensus features (\"max{(element / center), (center / element)}\", weight 1 per consensus features):"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(it_aad_by_cfs) << "\n";
      }
      else if (in_type == FileTypes::IDXML) //identifications
      {
        //TODO
      }
      else if (in_type == FileTypes::FASTA)
      {
      }
      else if (in_type == FileTypes::PEPXML)
      {
        // TODO
      }
      else //peaks
      {
        //copy intensities of  MS-level 1 peaks
        exp.updateRanges(1);
        Size size = exp.getSize();
        vector<double> intensities;
        intensities.reserve(size);
        for (PeakMap::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
        {
          if (spec->getMSLevel() != 1)
          {
            continue;
          }
          for (PeakMap::SpectrumType::const_iterator it = spec->begin(); it != spec->end(); ++it)
          {
            intensities.push_back(it->getIntensity());
          }
        }

        sort(intensities.begin(), intensities.end());
        os.precision(writtenDigits(Peak1D::IntensityType()));
        os << "Intensities:"
           << "\n"
           << Math::SummaryStatistics<vector<double>>(intensities) << "\n";

        //Statistics for meta information
        for (Map<String, int>::ConstIterator it = meta_names.begin(); it != meta_names.end(); ++it)
        {
          String name = it->first;
          vector<double> m_values;
          for (PeakMap::const_iterator spec = exp.begin(); spec != exp.end(); ++spec)
          {
            for (Size meta = 0; meta < spec->getFloatDataArrays().size(); ++meta)
            {
              if (spec->getFloatDataArrays()[meta].getName() != name)
                continue;
              for (Size peak = 0; peak < spec->getFloatDataArrays()[meta].size(); ++peak)
              {
                m_values.push_back(spec->getFloatDataArrays()[meta][peak]);
              }
            }
            for (Size meta = 0; meta < spec->getIntegerDataArrays().size(); ++meta)
            {
              if (spec->getIntegerDataArrays()[meta].getName() != name)
              {
                continue;
              }
              for (Size peak = 0; peak < spec->getIntegerDataArrays()[meta].size(); ++peak)
              {
                m_values.push_back(spec->getIntegerDataArrays()[meta][peak]);
              }
            }
          }
          os << "Meta data: " << name << "\n"
             << Math::SummaryStatistics<vector<double>>(m_values) << "\n";
        }
      }
    }

    os << "\n"
       << "\n";

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char **) override
  {
    String out = getStringOption_("out");
    String out_tsv = getStringOption_("out_tsv");

    TOPPBase::ExitCodes ret;
    if (out != "" && out_tsv != "")
    {
      ofstream os(out.c_str());
      ofstream os_tsv(out_tsv.c_str());
      ret = outputTo_(os, os_tsv);
      os.close();
      os_tsv.close();
    }
    else if (out != "" && out_tsv == "")
    {
      ofstream os(out.c_str());
      // Output stream with null output
      boost::iostreams::filtering_ostream os_tsv;
      os_tsv.push(boost::iostreams::null_sink());
      ret = outputTo_(os, os_tsv);
      os.close();
    }
    else if (out == "" && out_tsv != "")
    {
      ofstream os_tsv(out_tsv.c_str());
      ret = outputTo_(LOG_INFO, os_tsv);
      os_tsv.close();
    }
    else
    {
      // Output stream with null output
      boost::iostreams::filtering_ostream os_tsv;
      os_tsv.push(boost::iostreams::null_sink());
      ret = outputTo_(LOG_INFO, os_tsv);
    }
    return ret;
  }
};


int main(int argc, const char **argv)
{
  TOPPFileInfo tool;
  return tool.main(argc, argv);
}

/// @endcond
