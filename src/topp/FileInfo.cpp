// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Marc Sturm, Clemens Groepl, Lars Nilse, Chris Bielow $
// --------------------------------------------------------------------------

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <OpenMS/config.h>



#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h> // for operator<< on StringList
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/SYSTEM/SysInfo.h>


#include <unordered_map>
#include <iomanip>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FileInfo FileInfo
@brief Shows basic information about the data in an %OpenMS readable file.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; FileInfo &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format) </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none ; console or text file</td>
</tr>
</table>
</CENTER>

This tool can show basic information about the data in different file types, such as raw peak, featureXML and consensusXML files. It can
- show information about the data range of a file (m/z, RT, ion mobility, intensity)
- show a statistical summary for intensities, qualities, feature widths, precursor charges, activation methods
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
  /// iterate over all items in CONTAINER and
  /// use LAMBDA function to extract the charge from each item
  template <class CONTAINER, typename LAMBDA>
  void printChargeDistribution(const CONTAINER& data, LAMBDA lam, ostream& os, ostream& os_tsv, const String& header = "Charge")
  {
    std::map<Int, UInt> charges;
    Int q;
    for (const auto& item : data)
    {
      if (lam(item, q)) ++charges[q];
    }

    os << header << " distribution:"
      << '\n';
    for (const auto& ch : charges)
    {
      os << "  charge " << ch.first << ": " << ch.second << "x\n";
      os_tsv << "general: charge distribution: charge: "
        << ch.first << '\t'
        << ch.second << '\n';
    }
    os << '\n';


  };
  

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
    return os << "  num. of values: " << rhs.count << '\n'
              << "  mean:           " << rhs.mean << '\n'
              << "  minimum:        " << rhs.min << '\n'
              << "  lower quartile: " << rhs.lowerq << '\n'
              << "  median:         " << rhs.median << '\n'
              << "  upper quartile: " << rhs.upperq << '\n'
              << "  maximum:        " << rhs.max << '\n'
              << "  variance:       " << rhs.variance << '\n';
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
    StringList in_types = { "mzData", "mzXML", "mzML", "sqMass", "dta", "dta2d", "mgf", "featureXML", "consensusXML", "idXML", "pepXML", "mzTab", "fid", "mzid", "trafoXML", "fasta", "pqp" };
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", in_types);
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content", false);
    setValidStrings_("in_type", in_types);
    registerOutputFile_("out", "<file>", "", "Optional output file. If left out, the output is written to the command line.", false);
    setValidFormats_("out", {"txt"});
    registerOutputFile_("out_tsv", "<file>", "", "Second optional output file. Tab separated flat text file.", false, true);
    setValidFormats_("out_tsv", {"tsv"});
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
    os << "Ranges:" << '\n'
       << "  retention time: " << String::number(map.getMinRT(), 2) << " .. " << String::number(map.getMaxRT(), 2) << " sec ("
       << String::number((map.getMaxRT() - map.getMinRT()) / 60, 1) << " min)\n"
       << "  mass-to-charge: " << String::number(map.getMinMZ(), 2) << " .. " << String::number(map.getMaxMZ(), 2) << '\n';
    if constexpr (std::is_base_of < RangeMobility, Map>())
    {
      os << "    ion mobility: ";
      if (map.RangeMobility::isEmpty()) os << "<none>\n";
      else os << String::number(map.getMinMobility(), 2) << " .. " << String::number(map.getMaxMobility(), 2) << '\n';
    }
    os << "       intensity: " << String::number(map.getMinIntensity(), 2) << " .. " << String::number(map.getMaxIntensity(), 2) << '\n'
       << '\n';
  }

  template <class Map>
  void writeRangesMachineReadable_(const Map& map, ostream &os)
  {
    os << "general: ranges: retention time: min" << '\t' << String::number(map.getMinRT(), 2) << '\n'
       << "general: ranges: retention time: max" << '\t' << String::number(map.getMaxRT(), 2) << '\n'
       << "general: ranges: mass-to-charge: min" << '\t' << String::number(map.getMinMZ(), 2) << '\n'
       << "general: ranges: mass-to-charge: max" << '\t' << String::number(map.getMaxMZ(), 2) << '\n';
    if constexpr (std::is_base_of < RangeMobility, Map>())
    {
      if (!map.RangeMobility::isEmpty())
      {
        os << "general: ranges: ion-mobility: min" << '\t' << String::number(map.getMinMobility(), 2) << '\n'
           << "general: ranges: ion-mobility: max" << '\t' << String::number(map.getMaxMobility(), 2) << '\n';
      }
    }
    os << "general: ranges: intensity: min"
       << '\t' << String::number(map.getMinIntensity(), 2) << '\n'
       << "general: ranges: intensity: max"
       << '\t' << String::number(map.getMaxIntensity(), 2) << '\n';
  }

  template <class T>
  void writeSummaryStatisticsMachineReadable_(const Math::SummaryStatistics<T> &stats, ostream &os, String title)
  {
    os << "statistics: " << title << ": num. of values" << '\t' << stats.count    << '\n'
       << "statistics: " << title << ": mean"           << '\t' << stats.mean     << '\n'
       << "statistics: " << title << ": minimum"        << '\t' << stats.min      << '\n'
       << "statistics: " << title << ": lower quartile" << '\t' << stats.lowerq   << '\n'
       << "statistics: " << title << ": median"         << '\t' << stats.median   << '\n'
       << "statistics: " << title << ": upper quartile" << '\t' << stats.upperq   << '\n'
       << "statistics: " << title << ": maximum"        << '\t' << stats.max      << '\n'
       << "statistics: " << title << ": variance"       << '\t' << stats.variance << '\n';
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
      in_type = FileHandler::getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    os << '\n'
       << "-- General information --"
       << '\n'
       << '\n'
       << "File name: " << in << '\n'
       << "File type: " << FileTypes::typeToName(in_type) << '\n';

    os_tsv << "general: file name"
           << '\t' << in << '\n'
           << "general: file type"
           << '\t' << FileTypes::typeToName(in_type) << '\n';

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
      os << '\n'
         << "Validating " << FileTypes::typeToName(in_type) << " file";
      switch (in_type)
      {
      case FileTypes::MZDATA:
        os << " against XML schema version " << MzDataFile().getVersion() << '\n';
        valid = MzDataFile().isValid(in, os);
        break;

      case FileTypes::MZML:
        os << " against XML schema version " << MzMLFile().getVersion() << '\n';
        valid = MzMLFile().isValid(in, os);
        break;

      case FileTypes::FEATUREXML:
        os << " against XML schema version " << FeatureXMLFile().getVersion() << '\n';
        valid = FeatureXMLFile().isValid(in, os);
        break;

      case FileTypes::IDXML:
        os << " against XML schema version " << IdXMLFile().getVersion() << '\n';
        valid = IdXMLFile().isValid(in, os);
        break;

      case FileTypes::MZIDENTML:
        os << " against XML schema version " << MzIdentMLFile().getVersion() << '\n';
        valid = MzIdentMLFile().isValid(in, os);
        break;

      case FileTypes::CONSENSUSXML:
        os << " against XML schema version " << ConsensusXMLFile().getVersion() << '\n';
        valid = ConsensusXMLFile().isValid(in, os);
        break;

      case FileTypes::MZXML:
        os << " against XML schema version " << MzXMLFile().getVersion() << '\n';
        valid = MzXMLFile().isValid(in, os);
        break;

      case FileTypes::PEPXML:
        os << " against XML schema version " << PepXMLFile().getVersion() << '\n';
        valid = PepXMLFile().isValid(in, os);
        break;

      case FileTypes::TRANSFORMATIONXML:
        os << " against XML schema version " << TransformationXMLFile().getVersion() << '\n';
        valid = TransformationXMLFile().isValid(in, os);
        break;

      default:
        os << '\n'
           << "Aborted: Validation of this file type is not supported!"
           << '\n';
        return EXECUTION_OK;
      }

      if (valid)
      {
        os << "Success - the file is valid!"
           << '\n';
      }
      else
      {
        os << "Failed - errors are listed above!"
           << '\n';
      }

      // semantic validation:
      if ((in_type == FileTypes::MZML) || (in_type == FileTypes::MZDATA))
      {
        if (!valid)
        {
          os << '\n'
             << "Semantic validation is not performed due to previous errors!"
             << '\n';
        }
        else
        {
          os << '\n'
             << "Semantically validating " << FileTypes::typeToName(in_type)
             << " file";
          if (in_type == FileTypes::MZDATA)
            os << " (EXPERIMENTAL)";
          os << ":"
             << '\n';

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
            os << "Warning: " << warnings[i] << '\n';
          }
          for (Size i = 0; i < errors.size(); ++i)
          {
            os << "Error: " << errors[i] << '\n';
          }
          if (valid)
          {
            os << "Success - the file is semantically valid!"
               << '\n';
          }
          else
          {
            os << "Failed - errors are listed above!"
               << '\n';
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
        writeLogError_("Error: Can only validate indices for mzML files");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      std::cout << "Checking mzML file for valid indices ... " << std::endl;
      Internal::IndexedMzMLHandler ifile;
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
    std::map<String, int> meta_names;

    if (in_type == FileTypes::FASTA)
    {
      vector<FASTAFile::FASTAEntry> entries;
      FASTAFile file;
      file.setLogType(log_type_);

      SysInfo::MemUsage mu;
      // loading input
      file.load(in, entries);
      std::cout << "\n\n" << mu.delta("loading FASTA") << std::endl;

      std::map<char, int> aacids;// required for default construction of non-existing keys
      size_t number_of_aacids = 0;

      Size dup_header(0);
      Size dup_seq(0);

      typedef std::unordered_map<size_t, vector<ptrdiff_t> > SHashmap;
      SHashmap m_headers;
      SHashmap m_seqs;

      std::hash<string> s_hash;
      for (auto loopiter = entries.begin(); loopiter != entries.end(); ++loopiter)
      {
        {
          size_t id_hash = s_hash(loopiter->identifier);
          auto it_id = m_headers.find(id_hash);
          if (it_id != m_headers.end())
          { // hash matches ... test the real indices to make sure
            vector<ptrdiff_t>::const_iterator iter = find_if(it_id->second.begin(), it_id->second.end(), [&loopiter, &entries](const ptrdiff_t& idx) { return entries[idx].headerMatches(*loopiter); });
            if (iter != it_id->second.end())
            {
              os << "Warning: Duplicate header, #" << std::distance(entries.begin(), loopiter) << ", ID: " << loopiter->identifier << " = #" << *iter << ", ID: " << entries[*iter].identifier << '\n';
              ++dup_header;
            }

          }
          // add our own hash
          m_headers[id_hash] = { std::distance(entries.begin(), loopiter) };
        }

        {
          size_t id_seq = s_hash(loopiter->sequence);
          auto it_id = m_seqs.find(id_seq);
          if (it_id != m_seqs.end())
          { // hash matches ... test the real indices to make sure
            vector<ptrdiff_t>::const_iterator iter = find_if(it_id->second.begin(), it_id->second.end(), [&loopiter, &entries](const ptrdiff_t& idx) { return entries[idx].sequenceMatches(*loopiter); });
            if (iter != it_id->second.end())
            {
              os << "Warning: Duplicate sequence, #" << std::distance(entries.begin(), loopiter) << ", ID: " << loopiter->identifier << " == #" << *iter << ", ID: " << entries[*iter].identifier << '\n';
              ++dup_seq;
            }

          }
          // add our own hash
          m_seqs[id_seq] = { std::distance(entries.begin(), loopiter) };
        }

        for (char a : loopiter->sequence)
        {
          ++aacids[a];
        }
        number_of_aacids += loopiter->sequence.size();
      }

      os << '\n';
      os << "Number of sequences   : " << entries.size() << '\n';
      os << "# duplicated headers  : " << dup_header << " (" << (entries.empty() ? 0 :
                                                                 static_cast<Size>(dup_header * 1000 / entries.size()) / 10.0) << "%)\n";
      os << "# duplicated sequences: " << dup_seq << " (" << (entries.empty() ? 0 : static_cast<Size>(dup_seq * 1000 / entries.size()) / 10.0) << "%) [by exact string matching]\n";
      os << "Total amino acids     : " << number_of_aacids << "\n\n";
      os << "Amino acid counts: \n";

      for (auto it = aacids.begin(); it != aacids.end(); ++it)
      {
        os << it->first << '\t' << it->second << '\n';
      }
      size_t amb = aacids['B'] + aacids['Z'] + aacids['X'] + aacids['b'] + aacids['z'] + aacids['x'];
      size_t amb_I = amb + aacids['I'] + aacids['i'];
      os << "Ambiguous amino acids (B/Z/X)  : " << amb   << " (" << (amb > 0 ? (static_cast<Size>(amb * 10000 / number_of_aacids) / 100.0) : 0) << "%)\n";
      os << "                      (B/Z/X/I): " << amb_I << " (" << (amb_I > 0 ? (static_cast<Size>(amb_I * 10000 / number_of_aacids) / 100.0) : 0) << "%)\n\n";
    }
    else if (in_type == FileTypes::FEATUREXML) //features
    {
      FileHandler ff;
      ff.getFeatOptions().setLoadConvexHull(false);   // CH's currently not needed here
      ff.getFeatOptions().setLoadSubordinates(false); // SO's currently not needed here

      SysInfo::MemUsage mu;
      // reading input
      ff.loadFeatures(in, feat, {FileTypes::FEATUREXML});
      std::cout << "\n\n" << mu.delta("loading featureXML") << std::endl;

      feat.updateRanges();

      os << "Number of features: " << feat.size() << '\n'
         << '\n';
      os_tsv << "general: number of features" << '\t'
             << feat.size() << '\n';

      writeRangesHumanReadable_(feat, os);
      writeRangesMachineReadable_(feat, os_tsv);

      // Charge distribution and TIC
      std::map<Int, UInt> charges;
      std::map<size_t, UInt> numberofids;
      double tic = 0.0;
      for (Size i = 0; i < feat.size(); ++i)
      {
        ++charges[feat[i].getCharge()];
        tic += feat[i].getIntensity();
        const vector<PeptideIdentification> &peptide_ids = feat[i].getPeptideIdentifications();
        ++numberofids[peptide_ids.size()];
      }

      os << "Total ion current in features: " << tic << '\n';
      os_tsv << "general: total ion current in features" << '\t'
             << tic << '\n';
      os << '\n';

      printChargeDistribution(feat, [](const Feature& f, Int& q) { q = f.getCharge(); return true;}, os, os_tsv);

      os << "Distribution of peptide identifications (IDs) per feature:\n";
      for (auto it = numberofids.begin(); it != numberofids.end(); ++it)
      {
        os << "  " << it->first << " IDs: " << it->second << '\n';
        os_tsv << "general: distribution of peptide identifications (IDs) per feature: IDs: "
               << it->first << '\t'
               << it->second << '\n';
      }

      os << '\n'
         << "Unassigned peptide identifications: " << feat.getUnassignedPeptideIdentifications().size() << '\n';
      os_tsv << "general: unassigned peptide identifications" << '\t'
             << feat.getUnassignedPeptideIdentifications().size() << '\n';
    }
    else if (in_type == FileTypes::CONSENSUSXML) //consensus features
    {

      SysInfo::MemUsage mu;
      // reading input
      FileHandler().loadConsensusFeatures(in, cons, {FileTypes::CONSENSUSXML});
      std::cout << "\n\n" << mu.delta("loading consensusXML") << std::endl;

      cons.updateRanges();

      map<Size, UInt> num_consfeat_of_size;
      map<Size, UInt> num_consfeat_of_size_with_id;

      map<pair<String, UInt>, vector<int> > seq_charge2map_occurence;
      for (const ConsensusFeature& cm : cons)
      {
        ++num_consfeat_of_size[cm.size()];
        const auto& pids = cm.getPeptideIdentifications();
        if (!pids.empty())
        {
          ++num_consfeat_of_size_with_id[cm.size()];

          // count how often a peptide/charge pair has been observed in the different maps
          const vector<PeptideHit>& phits = pids[0].getHits();
          if (!phits.empty())
          {
            const String s = phits[0].getSequence().toString();
            const int z = phits[0].getCharge();

            if (seq_charge2map_occurence[make_pair(s,z)].empty())
            {
              seq_charge2map_occurence[make_pair(s,z)] = vector<int>(cons.getColumnHeaders().size(), 0);
            }

            // assign id to all dimensions in the consensus feature
            for (auto const & f : cm.getFeatures())
            {
              Size map_index = f.getMapIndex();
              seq_charge2map_occurence[make_pair(s,z)][map_index] += 1;
            }
          }
        }
      }

      // now at the level of peptides (different charges and modifications are counted separately) 
      // to get a number independent of potential alignment/link errors
      // Note: 
      // we determine the size of a consensus feature we would obtain if we would link just be sequence and charge
      // and we sum up all sub features for these consensus feature that has at least one id
      // (as we assume that the ID is transfered to all sub features)
      map<Size, Size> num_aggregated_consfeat_of_size_with_id;
      map<Size, Size> num_aggregated_feat_of_size_with_id;
      for (auto & a : seq_charge2map_occurence)
      {
        const vector<int>& occurrences = a.second;
        UInt n(0); // dimensions with at least one peptide id assigned
        UInt f(0); // number of subfeatures with a least one peptide id assigned
        for (int i : occurrences) 
        { 
          if (i != 0) ++n; 
          f += i;
        }
        num_aggregated_consfeat_of_size_with_id[n] += 1;
        num_aggregated_feat_of_size_with_id[n] += f;	
      }

      if (num_consfeat_of_size.empty())
      {
        os << '\n'
           << "Number of consensus features: 0"
           << '\n';
        os << "No consensus features found, map is empty!"
           << "\n\n";
      }
      else
      {
        Size field_width = num_consfeat_of_size.rbegin()->first / 10 + 1;
        os << '\n'
           << "Number of consensus features:"
           << '\n';
        
        Size number_features{0};
        Size number_cons_features_with_id{0};
        Size number_features_with_id{0};
        for (auto it = num_consfeat_of_size.rbegin(); it != num_consfeat_of_size.rend(); ++it)
        {
          const Size csize = it->first;
          const Size nfeatures = it->first * it->second;
          const Size nc_with_id = num_consfeat_of_size_with_id[it->first];
          number_features += nfeatures;
          number_features_with_id += csize * nc_with_id;
          number_cons_features_with_id += nc_with_id;

          os << "  of size " << setw(field_width) << csize << ": " << it->second 
             << "\t (features: " << nfeatures << " )"
             << "\t with at least one ID: " << nc_with_id
             << "\t (features: " << csize * nc_with_id << " )"
             << '\n';

        }

        auto ci = num_aggregated_consfeat_of_size_with_id.rbegin();
        auto fi = num_aggregated_feat_of_size_with_id.rbegin();
        for (; ci != num_aggregated_consfeat_of_size_with_id.rend(); ++ci, ++fi)
        {
          const Size csize = ci->first;
          const Size nconsfeatures = ci->second;
          const Size nfeatures = fi->second;

          os << "  peptides (with different mod. and charge) observed in " << setw(field_width) << csize << " maps: " << nconsfeatures 
             << "\t (features: " << nfeatures << " )"
             << '\n';

        }
        os << "  total consensus features:    "  << cons.size()
           << "  with at least one ID: " << number_cons_features_with_id << '\n'
           << "  total features:              " << number_features 
           << "  with at least one ID: " << string(field_width, ' ') << number_features_with_id
           << '\n';

        writeRangesHumanReadable_(cons, os);
        writeRangesMachineReadable_(cons, os_tsv);
      }

      // file descriptions
      const ConsensusMap::ColumnHeaders& descs = cons.getColumnHeaders();
      if (!descs.empty())
      {
        os << "File descriptions:"
           << '\n';
        for (ConsensusMap::ColumnHeaders::const_iterator it = descs.begin(); it != descs.end(); ++it)
        {
          os << "  " << it->second.filename << ":"
             << '\n'
             << "    identifier: " << it->first << '\n'
             << "    label:      " << it->second.label << '\n'
             << "    size:       " << it->second.size << '\n';
        }
        os << '\n';
      }
    }
    else if (in_type == FileTypes::IDXML || in_type == FileTypes::MZIDENTML) //identifications
    {
      UInt spectrum_count(0);
      Size peptide_hit_count(0);
      UInt runs_count(0);
      Size protein_hit_count(0);
      set<String> peptides;
      set<String> peptides_ignore_mods;
      set<String> proteins;
      Size modified_peptide_count(0);
      std::map<String, int> mod_counts;
      vector<uint16_t> peptide_length;

      // reading input
      SysInfo::MemUsage mu;
      if (in_type == FileTypes::MZIDENTML)
      {
        FileHandler().loadIdentifications(in, id_data.proteins, id_data.peptides, {FileTypes::MZIDENTML});
      }
      else
      {
        IdXMLFile().load(in, id_data.proteins, id_data.peptides, id_data.identifier);
      }
      std::cout << "\n\n" << mu.delta("loading idXML") << std::endl;

      // export metadata to second output stream
      os_tsv << "general: database"
             << '\t' << id_data.proteins[0].getSearchParameters().db << '\n'
             << "general: database version"
             << '\t' << id_data.proteins[0].getSearchParameters().db_version << '\n'
             << "general: taxonomy"
             << '\t' << id_data.proteins[0].getSearchParameters().taxonomy << '\n';

      // calculations
      Size average_peptide_hits{0}; // average number of hits per spectrum (ignoring the empty ones)
      for (Size i = 0; i < id_data.peptides.size(); ++i)
      {
        if (!id_data.peptides[i].empty())
        {
          ++spectrum_count;
          average_peptide_hits += id_data.peptides[i].getHits().size();
          peptide_hit_count += id_data.peptides[i].getHits().size();
          const vector<PeptideHit> &temp_hits = id_data.peptides[i].getHits();
          // collect stats about modifications from TOP HIT!
          if (temp_hits[0].getSequence().isModified())
          {
            ++modified_peptide_count;
            const AASequence& aa = temp_hits[0].getSequence();
            if (aa.hasCTerminalModification())
            {
              ++mod_counts[aa.getCTerminalModificationName()];
            }
            if (aa.hasNTerminalModification())
            {
              ++mod_counts[aa.getNTerminalModificationName()];
            }
            for (Size ia = 0; ia < aa.size(); ++ia)
            {
              if (aa[ia].isModified())
              {
                ++mod_counts[aa[ia].getModification()->getFullId()];
              }
            }
          }
          for (Size j = 0; j < temp_hits.size(); ++j)
          {
            peptides.insert(temp_hits[j].getSequence().toString());
            peptides_ignore_mods.insert(temp_hits[j].getSequence().toUnmodifiedString());
            peptide_length.push_back((uint16_t)temp_hits[j].getSequence().size());
          }
        }
      }
      set<pair<String, String>> search_engines;
      for (Size i = 0; i < id_data.proteins.size(); ++i)
      {
        ++runs_count;
        protein_hit_count += id_data.proteins[i].getHits().size();
        const vector<ProteinHit> &temp_hits = id_data.proteins[i].getHits();
        for (Size j = 0; j < temp_hits.size(); ++j)
        {
          proteins.insert(temp_hits[j].getAccession());
        }
        // collect all search engines which generated the data
        search_engines.emplace(id_data.proteins[i].getSearchEngine(), id_data.proteins[i].getSearchEngineVersion());
      }
      if (peptide_length.empty())
      { // avoid invalid-range exception when computing mean()
        peptide_length.push_back(0);
      }
      
      os << "Search Engine(s):\n";
      for (const auto& se : search_engines)
      {
        os << "  " << se.first << " (version: " << se.second << ")\n";
      }
      os << "Number of:"
         << '\n';
      os << "  runs:                       " << runs_count << '\n';
      os << "  protein hits:               " << protein_hit_count << '\n';
      os << "  non-redundant protein hits: " << proteins.size() << '\n';
      os << "  (only hits that differ in the accession)"
         << '\n';
      os << '\n';
      os << "  matched spectra:    " << spectrum_count << '\n';
      os << "  peptide sequences:  " << peptides_ignore_mods.size() << '\n';
      os << "  PSMs / spectrum (ignoring unidentified spectra):    " << average_peptide_hits / std::max(1, (Int)spectrum_count) << '\n';
      os << "  peptide hits:               " << peptide_hit_count << " (avg. length: " << Math::round(Math::mean(peptide_length.begin(), peptide_length.end())) << ")\n";
      os << "  modified top-hits:          " << modified_peptide_count << "/" << spectrum_count << (spectrum_count > 0 ? String(" (") + Math::round(modified_peptide_count * 1000.0 / spectrum_count) / 10 + "%)" : "") << '\n';
      os << "  non-redundant peptide hits: " << peptides.size() << '\n';
      os << "  (only hits that differ in sequence and/or modifications)"
         << '\n';
      for (std::map<String, int>::const_iterator it = mod_counts.begin(); it != mod_counts.end(); ++it)
      {
        if (it != mod_counts.begin())
        {
          os << ", ";
        }
        else
        {
          os << "  Modification count (top-hits only): ";
        }
        os << it->first << " " << it->second;
      }

      for (const auto& se : search_engines)
      {
        os_tsv << "general: search engine" << '\t' << se.first << '\t' << "(version: " << se.second << ")" << '\n';
      }
      os_tsv << "general: num. of runs" << '\t' << runs_count << '\n';
      os_tsv << "general: num. of protein hits" << '\t' << protein_hit_count << '\n';
      os_tsv << "general: num. of non-redundant protein hits (only hits that differ in the accession)"
             << '\t' << proteins.size() << '\n';
      os_tsv << "general: num. of matched spectra" << '\t' << spectrum_count << '\n';
      os_tsv << "general: num. of peptide hits" << '\t' << peptide_hit_count << '\n';
      os_tsv << "general: num. of modified top-hits" << '\t' << modified_peptide_count << '\n';
      os_tsv << "general: num. of non-redundant peptide hits (only hits that differ in sequence and/or modifications): "
             << '\t' << peptides.size() << '\n';
    }
    else if (in_type == FileTypes::PEPXML)
    {
      os << "\nFor pepXML files, only validation against the XML schema is implemented at this point."
         << '\n';
    }
    else if (in_type == FileTypes::MZTAB)
    {
      MzTab mztab;
      MzTabFile().load(in, mztab);
      os << "mzTab-version: " << mztab.getMetaData().mz_tab_version.get() << '\n'
         << "mzTab-mode: " << mztab.getMetaData().mz_tab_mode.get() << '\n'
         << "mzTab-type: " << mztab.getMetaData().mz_tab_type.get() << '\n'
         << "number of PSMs: " << mztab.getNumberOfPSMs() << '\n'
         << "number of peptides: " << mztab.getPeptideSectionRows().size() << '\n'
         << "number of proteins: " << mztab.getProteinSectionRows().size() << '\n'
         << "number of oligonucleotides: " << mztab.getOligonucleotideSectionRows().size() << '\n'
         << "number of OSMs: " << mztab.getOSMSectionRows().size() << '\n'
         << "number of small molecules: " << mztab.getSmallMoleculeSectionRows().size() << '\n'
         << "number of nucleic acids: " << mztab.getNucleicAcidSectionRows().size() << '\n';
    }
    else if (in_type == FileTypes::TRANSFORMATIONXML)
    {
      TransformationDescription trafo;
      FileHandler().loadTransformations(in, trafo, true, {FileTypes::TRANSFORMATIONXML});
      os << "\nTransformation model: " << trafo.getModelType() << '\n';
      trafo.printSummary(os);
    }
    else if (in_type == FileTypes::PQP)
    {
      TargetedExperiment targeted_exp; 
      TransitionPQPFile pqp_reader;
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertPQPToTargetedExperiment(in.c_str(), targeted_exp, true);
      os << targeted_exp.getSummary();
    }
    else // peaks
    {
      SysInfo::MemUsage mu;
      fh.loadExperiment(in, exp, {in_type}, log_type_, false, false);

      // update range information and retrieve which MS levels were recorded
      exp.updateRanges();
      // report memory consumption
      OPENMS_LOG_INFO << "\n\n" << mu.delta("loading MS data") << std::endl;

      os << '\n';

      os << "Instrument: " << exp.getInstrument().getName() << '\n';
      for (const auto& ma : exp.getInstrument().getMassAnalyzers())
      {
        os << "  Mass Analyzer: " << MassAnalyzer::NamesOfAnalyzerType[ma.getType()] << " (resolution: " << ma.getResolution() << ")\n";
      }
      os << '\n';

      const vector<UInt>& levels = exp.getMSLevels();
      os << "MS levels: " << ListUtils::concatenate(levels, ", ") << '\n';

      // basic info
      os << "Total number of peaks: " << exp.getSize() << '\n'; // count ALL peaks (also chromatographic)
      os << "Number of spectra: " << exp.size() << '\n'
        << '\n';
      os_tsv << "number of spectra"
        << '\t' << exp.size() << '\n'
        << "total number of peaks"
        << '\t' << exp.getSize() << '\n';

      writeRangesHumanReadable_(exp, os);
      writeRangesMachineReadable_(exp, os_tsv);
      // check if the meta data indicates that this is peak data
      // and count how many spectra per MS level there are
      map<Size, UInt> level_annotated_picked;
      map<Size, UInt> level_estimated_picked;
      struct ChAM
      {
        Size mslevel;
        Precursor::ActivationMethod am;
        bool operator<(const ChAM& rhs) const
        {
          return std::tie(mslevel, am) < std::tie(rhs.mslevel, rhs.am); 
        }
      };
      map<ChAM, Size> act_method_counts;
      map<Size, UInt> counts;
      for (const auto& spectrum : exp)
      {
        const Size level = spectrum.getMSLevel();
        ++counts[level];  // count MS level

        for (auto const& pc : spectrum.getPrecursors())
        {
          for (auto const& am : pc.getActivationMethods())
          {
            ++act_method_counts[{level, am}];
          }
        }

        // annotate peak type (profile / centroided) from meta data
        if (level_annotated_picked.count(level) == 0)
        {
          level_annotated_picked[level] = spectrum.getType(false);
        }

        // estimate peak type once for every level (take a spectrum with enough peaks for stable estimation)
        if (level_estimated_picked.count(level) == 0 && spectrum.size() > 10)
        {
          level_estimated_picked[level] = PeakTypeEstimator::estimateType(spectrum.begin(), spectrum.end());
        }
      }

      // output
      if (!counts.empty())
      {
        os << "Number of spectra per MS level:"
          << '\n';
        for (auto it = counts.begin(); it != counts.end(); ++it)
        {
          os << "  level " << it->first << ": " << it->second << '\n';
          os_tsv << "number of MS" << it->first << " spectra"
            << '\t' << it->second << '\n';
        }
        os << '\n';
      }

      // write peak types (centroided / profile mode)
      os << "Peak type from metadata (or estimated from data)\n";
      for (const auto& l : levels)
      {
        os << "  level " << l << ": "
           << SpectrumSettings::NamesOfSpectrumType[level_annotated_picked[l]] << " ("
           << SpectrumSettings::NamesOfSpectrumType[level_estimated_picked[l]] << ")\n";
        os_tsv << "peak type metadata [annotation, estimate]" << '\t' << SpectrumSettings::NamesOfSpectrumType[level_annotated_picked[l]] << '\t' << SpectrumSettings::NamesOfSpectrumType[level_estimated_picked[l]] << '\n';
      }
      os << '\n';

      os << "Activation methods\n";
      for (const auto& am : act_method_counts)
      {
        os << "    MS-Level " << am.first.mslevel << " & " << Precursor::NamesOfActivationMethodShort[am.first.am] << " (" << Precursor::NamesOfActivationMethod[am.first.am] << "): " << am.second << '\n';
        os_tsv << "activation methods (mslevel, method, count)" << '\t' << am.first.mslevel << '\t' << Precursor::NamesOfActivationMethodShort[am.first.am] << '\t' << am.second << '\n';
      }
      os << '\n';

      printChargeDistribution(exp,
                              [](const MSSpectrum& spec, int& q) { 
                                for (const auto& pc: spec.getPrecursors())
                                {
                                  q = pc.getCharge();
                                  return true;
                                }
                                return false;
                              },
                              os,
                              os_tsv,
                              "Precursor charge");


      // show meta data array names
      for (MSSpectrum& pm : exp)
      {
        for (Size i = 0; i < pm.getFloatDataArrays().size(); ++i)
        {
          ++meta_names[pm.getFloatDataArrays()[i].getName()];
        }
        for (Size i = 0; i < pm.getIntegerDataArrays().size(); ++i)
        {
          ++meta_names[pm.getIntegerDataArrays()[i].getName()];
        }
        for (Size i = 0; i < pm.getStringDataArrays().size(); ++i)
        {
          ++meta_names[pm.getStringDataArrays()[i].getName()];
        }
      }
      if (!meta_names.empty())
      {
        // nice formatting:
        Size max_length = 0;
        for (std::map<String, int>::const_iterator it = meta_names.begin(); it != meta_names.end(); ++it)
        {
          if (it->first.size() > max_length)
          {
            max_length = it->first.size();
          }
        }
        os << "Meta data array:"
           << '\n';
        for (std::map<String, int>::const_iterator it = meta_names.begin(); it != meta_names.end(); ++it)
        {
          String padding(max_length - it->first.size(), ' ');
          os << "  " << it->first << ": " << padding << it->second << " spectra"
             << '\n';
        }
        os << '\n';  
      }
      
      auto cvs = FAIMSHelper::getCompensationVoltages(exp);
      if (!cvs.empty())
      {        
        os << "IM (FAIMS_CV): ";
        StringList cvs_sl;
        for (double cv : cvs)
        {
          cvs_sl.push_back(String(cv));
        }
        os << cvs_sl << "\n\n";
      }

      // some chromatogram information
      if (!exp.getChromatograms().empty())
      {
        os << "Number of chromatograms: " << exp.getChromatograms().size() << '\n';
        os_tsv << "number of chromatograms"
               << '\t' << exp.getChromatograms().size() << '\n';

        Size num_chrom_peaks(0);
        std::map<ChromatogramSettings::ChromatogramType, Size> chrom_types;
        for (const MSChromatogram& ms : exp.getChromatograms())
        {
          num_chrom_peaks += ms.size();
          ++chrom_types[ms.getChromatogramType()];
        }
        os << "Number of chromatographic peaks: " << num_chrom_peaks << '\n'
           << '\n';
        os_tsv << "number of chromatographic peaks" << '\t' << num_chrom_peaks << '\n';

        os << "Number of chromatograms per type: "
           << '\n';
        for (std::map<ChromatogramSettings::ChromatogramType, Size>::const_iterator it = chrom_types.begin(); it != chrom_types.end(); ++it)
        {
          os << String("  ") + ChromatogramSettings::ChromatogramNames[it->first] + ":                         "
             << it->second << '\n';
        }
        if (getFlag_("d") && chrom_types.find(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM) != chrom_types.end())
        {
          os << '\n'
             << " -- Detailed chromatogram listing -- "
             << '\n';
          os << "\nSelected Reaction Monitoring Transitions:"
             << '\n';
          os << "Q1 Q3 RT_begin RT_end name comment"
             << '\n';
          for (const MSChromatogram& ms : exp.getChromatograms())
          {
            if (ms.getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
            {
              os << ms.getPrecursor().getMZ() << " " << ms.getProduct().getMZ() << " " << ms.front().getRT() << " " << ms.back().getRT() << " " << ms.getName() << " " << ms.getComment() << '\n';
            }
          }
        }
      }

      // Detailed listing of scans
      if (getFlag_("d") && !exp.empty())
      {
        os << '\n'
           << "-- Detailed spectrum listing --"
           << '\n';
        UInt count = 0;
        for (auto const& spectrum : exp)
        {
          ++count;
          os << '\n'
             << "Spectrum " << count << ":"
             << '\n'
             << "  mslevel:    " << spectrum.getMSLevel() << '\n'
             << "  scanMode:   " << InstrumentSettings::NamesOfScanMode[spectrum.getInstrumentSettings().getScanMode()] << '\n'
             << "  peaks:      " << spectrum.size() << '\n'
             << "  RT:         " << spectrum.getRT() << '\n'
             << "  m/z:        ";

          if (!spectrum.empty())
          {
            os << spectrum.begin()->getMZ() << " .. " << spectrum.rbegin()->getMZ() << '\n';
          }

          if (spectrum.getDriftTimeUnit() != DriftTimeUnit::NONE)
          {
            os << "  IM:         " <<  spectrum.getDriftTime() << ' '
                << spectrum.getDriftTimeUnitAsString()
                << '\n';
          }            

          os << "Precursors:  " << spectrum.getPrecursors().size() <<  '\n';

          auto pc_count = UInt{0};
          for (auto const& pc : spectrum.getPrecursors())
          {
            os << "Precursor[" << pc_count << "]\n"
               << "  charge: " << pc.getCharge() << '\n'
               << "  mz:     " << pc.getMZ() << '\n'
               << "  activation methods: \n";
            for (auto const& am : pc.getActivationMethods())
            {
              os << "    " << Precursor::NamesOfActivationMethodShort[am] << " (" << Precursor::NamesOfActivationMethod[am] << ")\n";
            }

            os << '\n';

            ++pc_count;
          }
        }
      }

      // Check for corrupt data
      if (getFlag_("c"))
      {
        os << '\n'
           << "-- Checking for corrupt data --"
           << '\n'
           << '\n';
        // RTs sorted?
        if (!exp.isSorted(false))
        {
          os << "Error: Spectrum retention times are not sorted in ascending order"
             << '\n';
        }
        vector<double> ms1_rts;
        ms1_rts.reserve(exp.size());
        for (Size s = 0; s < exp.size(); ++s)
        {
          // ms level = 0
          if (exp[s].getMSLevel() == 0)
          {
            os << "Error: MS-level 0 in spectrum (RT: " << exp[s].getRT() << ")"
               << '\n';
          }
          //scan size = 0
          if (exp[s].empty())
          {
            os << "Warning: No peaks in spectrum (RT: " << exp[s].getRT() << ")"
               << '\n';
          }
          //duplicate meta data array names
          std::map<String, int> names;
          for (Size m = 0; m < exp[s].getFloatDataArrays().size(); ++m)
          {
            String name = exp[s].getFloatDataArrays()[m].getName();
            if (names.find(name) != names.end())
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << '\n';
            }
            else
            {
              names[name] = 0;
            }
          }
          for (Size m = 0; m < exp[s].getIntegerDataArrays().size(); ++m)
          {
            String name = exp[s].getIntegerDataArrays()[m].getName();
            if (names.find(name) != names.end())
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << '\n';
            }
            else
            {
              names[name] = 0;
            }
          }
          for (Size m = 0; m < exp[s].getStringDataArrays().size(); ++m)
          {
            String name = exp[s].getStringDataArrays()[m].getName();
            if (names.find(name) != names.end())
            {
              os << "Error: Duplicate meta data array name '" << name << "' in spectrum (RT: " << exp[s].getRT() << ")"
                 << '\n';
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
            os << "Error: Duplicate spectrum retention time: " << ms1_rts[i] << '\n';
        }
        //check peaks
        for (Size s = 0; s < exp.size(); ++s)
        {
          //peaks sorted?
          if (!exp[s].isSorted())
          {
            os << "Error: Peak m/z positions are not sorted in ascending order in spectrum (RT: " << exp[s].getRT() << ")"
               << '\n';
          }
          vector<double> mzs;
          mzs.reserve(exp[s].size());
          for (Size p = 0; p < exp[s].size(); ++p)
          {
            //negative intensity
            if (exp[s][p].getIntensity() < 0.0)
            {
              os << "Warning: Negative peak intensity peak (RT: " << exp[s].getRT() << " MZ: " << exp[s][p].getMZ() << " intensity: " << exp[s][p].getIntensity() << ")"
                 << '\n';
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
                 << '\n';
          }
        }
      }
    }

    //-------------------------------------------------------------
    // meta information
    //-------------------------------------------------------------
    if (getFlag_("m") || !getStringOption_("out_tsv").empty())
    {
      //basic info
      os << '\n'
         << "-- Meta information --"
         << '\n'
         << '\n';

      if (in_type == FileTypes::FEATUREXML) //features
      {
        os << "Document ID: " << feat.getIdentifier() << '\n'
           << '\n';
        os_tsv << "meta: document ID" << '\t'
               << feat.getIdentifier() << '\n';
      }
      else if (in_type == FileTypes::CONSENSUSXML) //consensus features
      {
        os << "Document ID: " << cons.getIdentifier() << '\n'
           << '\n';
      }
      else if (in_type == FileTypes::IDXML) //identifications
      {
        os << "Document ID: " << id_data.identifier << '\n'
           << '\n';
        os_tsv << "meta: document ID" << '\t'
               << id_data.identifier;
      }
      else if (in_type == FileTypes::PEPXML)
      {
        // TODO
      }
      else if (in_type == FileTypes::FASTA)
      {
        // TODO
      }
      else //peaks
      {
        os << "Document ID:        " << exp.getIdentifier() << '\n'
           << "Date:               " << exp.getDateTime().get() << '\n';
        os_tsv << "document id"
               << '\t' << exp.getIdentifier() << '\n'
               << "date"
               << '\t' << exp.getDateTime().get() << '\n';

        //basic info
        os << '\n'
           << "Sample:"
           << '\n'
           << "  name:             " << exp.getSample().getName() << '\n'
           << "  organism:         " << exp.getSample().getOrganism() << '\n'
           << "  comment:          " << exp.getSample().getComment() << '\n';
        os_tsv << "sample name"
               << '\t' << exp.getSample().getName() << '\n'
               << "sample organism"
               << '\t' << exp.getSample().getOrganism() << '\n'
               << "sample comment"
               << '\t' << exp.getSample().getComment() << '\n';

        //instrument info
        os << '\n'
           << "Instrument:"
           << '\n'
           << "  name:             " << exp.getInstrument().getName() << '\n'
           << "  model:            " << exp.getInstrument().getModel() << '\n'
           << "  vendor:           " << exp.getInstrument().getVendor() << '\n'
           << "  ion source(s):    ";
        os_tsv << "instrument name"
               << '\t' << exp.getInstrument().getName() << '\n'
               << "instrument model"
               << '\t' << exp.getInstrument().getModel() << '\n'
               << "instrument vendor"
               << '\t' << exp.getInstrument().getVendor() << '\n';
        for (Size i = 0; i < exp.getInstrument().getIonSources().size(); ++i)
        {
          os << IonSource::NamesOfIonizationMethod[exp.getInstrument().getIonSources()[i].getIonizationMethod()];
          if (i != exp.getInstrument().getIonSources().size() - 1)
          {
            os << ", ";
          }
        }
        os << '\n'
           << "  mass analyzer(s): ";
        for (Size i = 0; i < exp.getInstrument().getMassAnalyzers().size(); ++i)
        {
          os << MassAnalyzer::NamesOfAnalyzerType[exp.getInstrument().getMassAnalyzers()[i].getType()];
          if (i != exp.getInstrument().getMassAnalyzers().size() - 1)
          {
            os << ", ";
          }
        }
        os << '\n'
           << "  detector(s):      ";
        for (Size i = 0; i < exp.getInstrument().getIonDetectors().size(); ++i)
        {
          os << IonDetector::NamesOfType[exp.getInstrument().getIonDetectors()[i].getType()];
          if (i != exp.getInstrument().getIonDetectors().size() - 1)
            os << ", ";
        }
        os << '\n'
           << '\n';

        //contact persons
        for (Size i = 0; i < exp.getContacts().size(); ++i)
        {
          os << "Contact person:"
             << '\n'
             << "  first name:     " << exp.getContacts()[i].getFirstName() << '\n'
             << "  last name:      " << exp.getContacts()[i].getLastName() << '\n'
             << "  email:          " << exp.getContacts()[i].getEmail() << '\n'
             << '\n';
        }
      }
    }

    //-------------------------------------------------------------
    // data processing
    //-------------------------------------------------------------
    if (getFlag_("p"))
    {
      //basic info
      os << '\n'
         << "-- Data processing information --"
         << '\n'
         << '\n';

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
             << '\n'
             << '\n';
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
           << '\n'
           << '\n';
      }
      else
      {
        for (Size i = 0; i < dp.size(); ++i)
        {
          os << "Processing " << (i + 1) << ":"
             << '\n';
          os << "  software name:    " << dp[i].getSoftware().getName() << '\n';
          os << "  software version: " << dp[i].getSoftware().getVersion() << '\n';
          os << "  completion time:  " << dp[i].getCompletionTime().get() << '\n';

          os_tsv << "data processing: " << (i + 1)
                 << ": software name" << '\t'
                 << dp[i].getSoftware().getName() << '\n';
          os_tsv << "data processing: " << (i + 1)
                 << ": software version" << '\t'
                 << dp[i].getSoftware().getVersion() << '\n';
          os_tsv << "data processing: " << (i + 1)
                 << ": completion time" << '\t'
                 << dp[i].getCompletionTime().get() << '\n';

          os << "  actions:          ";
          os_tsv << "data processing: " << (i + 1)
                 << ": actions" << '\t';
          for (set<DataProcessing::ProcessingAction>::const_iterator it = dp[i].getProcessingActions().begin();
               it != dp[i].getProcessingActions().end(); ++it)
          {
            if (it != dp[i].getProcessingActions().begin())
            {
              os << ", ";
              os_tsv << ", ";
            }
            os << DataProcessing::NamesOfProcessingAction[*it];
            os_tsv << DataProcessing::NamesOfProcessingAction[*it];
          }
          os << '\n'
             << '\n';
          os_tsv << '\n';
        }
      }
    }

    //-------------------------------------------------------------
    // statistics
    //-------------------------------------------------------------
    if (getFlag_("s"))
    {
      os << '\n'
         << "-- Statistics --"
         << '\n'
         << '\n';

      if (in_type == FileTypes::FEATUREXML) //features
      {
        Size size = feat.size();

        vector<double> intensities(size);
        vector<double> overall_qualities(size);
        vector<double> mz_qualities(size);
        vector<double> rt_qualities(size);
        vector<double> peak_widths(size);

        Size idx = 0;
        for (const Feature& fm : feat)
        {
          intensities[idx] = fm.getIntensity();
          overall_qualities[idx] = fm.getOverallQuality();
          rt_qualities[idx] = fm.getQuality(Feature::RT);
          mz_qualities[idx] = fm.getQuality(Feature::MZ);
          peak_widths[idx] = fm.getWidth();
          ++idx;
        }

        Math::SummaryStatistics<vector<double>> intensities_summary;
        intensities_summary = Math::SummaryStatistics<vector<double>>(intensities);
        os.precision(writtenDigits<>(Feature::IntensityType()));
        os << "Intensities:" << '\n' << intensities_summary << '\n';
        os_tsv.precision(writtenDigits<>(Feature::IntensityType()));
        writeSummaryStatisticsMachineReadable_(intensities_summary, os_tsv, "intensities");

        Math::SummaryStatistics<vector<double>> peak_widths_summary;
        peak_widths_summary = Math::SummaryStatistics<vector<double>>(peak_widths);
        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Feature FWHM in RT dimension:" << '\n' << peak_widths_summary << '\n';
        os_tsv.precision(writtenDigits<>(Feature::QualityType()));
        writeSummaryStatisticsMachineReadable_(peak_widths_summary, os_tsv, "feature FWHM in RT dimension");

        Math::SummaryStatistics<vector<double>> overall_qualities_summary;
        overall_qualities_summary = Math::SummaryStatistics<vector<double>>(overall_qualities);
        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Overall qualities:" << '\n' << overall_qualities_summary << '\n';
        os_tsv.precision(writtenDigits<>(Feature::QualityType()));
        writeSummaryStatisticsMachineReadable_(overall_qualities_summary, os_tsv, "overall qualities");

        Math::SummaryStatistics<vector<double>> rt_qualities_summary;
        rt_qualities_summary = Math::SummaryStatistics<vector<double>>(rt_qualities);
        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Qualities in retention time dimension:" << '\n' << rt_qualities_summary << '\n';
        os_tsv.precision(writtenDigits<>(Feature::QualityType()));
        writeSummaryStatisticsMachineReadable_(rt_qualities_summary, os_tsv, "qualities in retention time dimension");

        Math::SummaryStatistics<vector<double>> mz_qualities_summary;
        mz_qualities_summary = Math::SummaryStatistics<vector<double>>(mz_qualities);
        os.precision(writtenDigits<>(Feature::QualityType()));
        os << "Qualities in mass-to-charge dimension:" << '\n' << mz_qualities_summary << '\n';
        os_tsv.precision(writtenDigits<>(Feature::QualityType()));
        writeSummaryStatisticsMachineReadable_(mz_qualities_summary, os_tsv, "qualities in mass-to-charge dimension");
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

        for (const ConsensusFeature& cm : cons)
        {
          double rt_aad = 0;
          double mz_aad = 0;
          double it_aad = 0;
          intensities.push_back(cm.getIntensity());
          qualities.push_back(cm.getQuality());
          widths.push_back(cm.getWidth());
          for (ConsensusFeature::HandleSetType::const_iterator hs_iter = cm.begin();
               hs_iter != cm.end(); ++hs_iter)
          {
            double rt_diff = hs_iter->getRT() - cm.getRT();
            rt_delta_by_elems.push_back(rt_diff);
            if (rt_diff < 0)
            {
              rt_diff = -rt_diff;
            }
            rt_aad_by_elems.push_back(rt_diff);
            rt_aad += rt_diff;
            double mz_diff = hs_iter->getMZ() - cm.getMZ();
            mz_delta_by_elems.push_back(mz_diff);
            if (mz_diff < 0)
            {
              mz_diff = -mz_diff;
            }
            mz_aad_by_elems.push_back(mz_diff);
            mz_aad += mz_diff;
            double it_ratio = hs_iter->getIntensity() / (cm.getIntensity() > 0 ? cm.getIntensity() : 1.);
            it_delta_by_elems.push_back(it_ratio);
            if (it_ratio < 1.)
            {
              it_ratio = 1. / it_ratio;
            }
            it_aad_by_elems.push_back(it_ratio);
            it_aad += it_ratio;
          }
          if (!cm.empty())
          {
            rt_aad /= cm.size();
            mz_aad /= cm.size();
            it_aad /= cm.size();
          } // otherwise rt_aad etc. are 0 anyway
          rt_aad_by_cfs.push_back(rt_aad);
          mz_aad_by_cfs.push_back(mz_aad);
          it_aad_by_cfs.push_back(it_aad);
        }

        os.precision(writtenDigits(ConsensusFeature::IntensityType()));
        os << "Intensities of consensus features:"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(intensities) << '\n';

        os.precision(writtenDigits(ConsensusFeature::QualityType()));
        os << "Qualities of consensus features:"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(qualities) << '\n';

        os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
        os << "Retention time differences (\"element - center\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(rt_delta_by_elems) << '\n';
        os << "Absolute retention time differences (\"|element - center|\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(rt_aad_by_elems) << '\n';
        os << "Average absolute differences of retention time within consensus features (\"|element - center|\", weight 1 per consensus features):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(rt_aad_by_cfs) << '\n';

        os.precision(writtenDigits(ConsensusFeature::CoordinateType()));
        os << "Mass-to-charge differences (\"element - center\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(mz_delta_by_elems) << '\n';
        os << "Absolute differences of mass-to-charge (\"|element - center|\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(mz_aad_by_elems) << '\n';
        os << "Average absolute differences of mass-to-charge within consensus features (\"|element - center|\", weight 1 per consensus features):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(mz_aad_by_cfs) << '\n';

        os.precision(writtenDigits(ConsensusFeature::IntensityType()));
        os << "Intensity ratios (\"element / center\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(it_delta_by_elems) << '\n';
        os << "Relative intensity error (\"max{(element / center), (center / element)}\", weight 1 per element):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(it_aad_by_elems) << '\n';
        os << "Average relative intensity error within consensus features (\"max{(element / center), (center / element)}\", weight 1 per consensus features):"
           << '\n'
           << Math::SummaryStatistics<vector<double>>(it_aad_by_cfs) << '\n';
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
           << '\n'
           << Math::SummaryStatistics<vector<double>>(intensities) << '\n';

        //Statistics for meta information
        for (std::map<String, int>::const_iterator it = meta_names.begin(); it != meta_names.end(); ++it)
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
          os << "Meta data: " << name << '\n'
             << Math::SummaryStatistics<vector<double>>(m_values) << '\n';
        }
      }
    }

    os << '\n'
       << '\n';

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char **) override
  {
    String out = getStringOption_("out");
    String out_tsv = getStringOption_("out_tsv");

    ofstream os;
    ofstream os_tsv;
    boost::iostreams::filtering_ostream os_filt;
    boost::iostreams::filtering_ostream os_tsv_filt;


    if (out.empty())
    {
      os_filt.push(OpenMS_Log_info);
    }
    else
    {
      os.open(out);
      if (!os)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out);
      }
      os_filt.push(os);
    }

    if (out_tsv.empty())
    {
      os_tsv_filt.push(boost::iostreams::null_sink());
    }
    else
    {
      os_tsv.open(out_tsv);
      if (!os_tsv)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_tsv);
      }
      os_tsv_filt.push(os_tsv);
    }

    return outputTo_(os_filt, os_tsv_filt);
  }
};


int main(int argc, const char **argv)
{
  TOPPFileInfo tool;
  return tool.main(argc, argv);
}

/// @endcond
