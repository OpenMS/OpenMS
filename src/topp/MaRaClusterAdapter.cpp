// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Leon Bichmann $
// $Authors: Mathew The, Leon Bichmann $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <cmath>
#include <string>
#include <set>
#include <map>

#include <QtCore/QProcess>
#include <boost/algorithm/clamp.hpp>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MaRaClusterAdapter MaRaClusterAdapter

@brief MaRaClusterAdapter facilitates the input to, the call of and output integration of MaRaCluster.
MaRaCluster (https://github.com/statisticalbiotechnology/maracluster) is a tool to apply unsupervised clustering of ms2 spectra from shotgun proteomics datasets.

@experimental This tool is work in progress and usage and input requirements might change.

<center>
  <table>
      <tr>
          <th ALIGN = "center"> pot. predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=2> &rarr; MaRaClusterAdapter &rarr;</td>
          <th ALIGN = "center"> pot. successor tools </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1>any signal-/preprocessing tool @n (in mzML format) </td>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MSGFPlusAdapter </td>
      </tr>
  </table>
</center>
<p>MaRaCluster is dependent on the input parameter pcut, which is the logarithm of the pvalue cutoff.
The default value is -10, lower values will result in smaller but purer clusters. If specified peptide search results
can be provided as idXML files and the MaRaCluster Adapter will annotate cluster ids as attributes to each peptide
identification, which will be outputed as a merged idXML. Moreover the merged idXML containing only scan numbers,
cluster ids and file origin can be outputed without prior peptide identification searches. The assigned cluster ids in
the respective idXML are equal to the scanindex of the produced clustered mzML.
</p>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MaRaClusterAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MaRaClusterAdapter.html

MaRaCluster is written by Matthew The (https://github.com/statisticalbiotechnology/maracluster
Copyright Matthew The <matthew.the@scilifelab.se>)
Cite Publication:
MaRaCluster: A Fragment Rarity Metric for Clustering Fragment Spectra in Shotgun Proteomics
Journal of proteome research, 2016, 15(3), pp 713-720 DOI: 10.1021/acs.jproteome.5b00749
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class MaRaClusterAdapter :
  public TOPPBase
{
public:
  MaRaClusterAdapter() :
    TOPPBase("MaRaClusterAdapter", "Facilitate input to MaRaCluster and reintegrate.", true,
                { // citation(s), specific for this tool
                 { "The M and Käll L", "MaRaCluster: A Fragment Rarity Metric for Clustering Fragment Spectra in Shotgun Proteomics", "J Proteome Res 2016; 15: 3", "10.1021/acs.jproteome.5b00749"}
                }
            )
  {
  }

protected:
  struct MaRaClusterResult
  {
    Int file_idx;
    Int scan_nr;

    MaRaClusterResult(const Int f, const Int s):
        file_idx (f),
        scan_nr (s)
    {
    }

    explicit MaRaClusterResult(StringList& row)
    {
      file_idx = row[0].toInt();
      scan_nr = row[1].toInt();
    }

    bool operator!=(const MaRaClusterResult& rhs) const
    {
      if (file_idx != rhs.file_idx || scan_nr != rhs.scan_nr)
      {
        return true;
      }
      return false;
    }

    bool operator<(const MaRaClusterResult& rhs) const
    {
      if (file_idx < rhs.file_idx || (file_idx == rhs.file_idx && scan_nr < rhs.scan_nr))
      {
        return true;
      }
      return false;
    }

    bool operator==(const MaRaClusterResult& rhs) const
    {
      return !(operator !=(rhs));
    }
  };
  

  void registerOptionsAndFlags_() override
  {
    static const bool is_required(true);
    static const bool is_advanced_option(true);
   
    //input 
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", is_required);
    setValidFormats_("in", ListUtils::create<String>("mzML,mgf"));
    registerInputFileList_("id_in", "<files>", StringList(), "Optional idXML Input file(s) in the same order as mzML files - for Maracluster Cluster annotation", !is_required);
    setValidFormats_("id_in", ListUtils::create<String>("idXML"));

    //output
    registerOutputFile_("out", "<file>", "", "Output file in idXML format", !is_required);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("consensus_out", "<file>", "", "Consensus spectra in mzML format", !is_required);
    setValidFormats_("consensus_out", ListUtils::create<String>("mzML"));
    registerStringOption_("output_directory", "<directory>", "", "Output directory for MaRaCluster original consensus output", false);


    //pvalue cutoff
    registerDoubleOption_("pcut", "<value>", -10.0, "log(p-value) cutoff, has to be < 0.0. Default: -10.0.", !is_required);
    setMaxFloat_("pcut", 0.0);
    registerIntOption_("min_cluster_size", "<value>", 1, "minimum number of spectra in a cluster for consensus spectra", !is_required);

    // minimal cluster size
    setMinInt_("min_cluster_size", 1);

    // executable
    registerInputFile_("maracluster_executable", "<executable>",
        // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "maracluster.exe",
        #else
                       "maracluster",
        #endif
                       "The maracluster executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", is_required, !is_advanced_option, {"is_executable"}
    );

    // Advanced parameters
    registerIntOption_("verbose", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all.", !is_required, is_advanced_option);
    registerDoubleOption_("precursor_tolerance", "<tolerance>", 20.0, "Precursor monoisotopic mass tolerance", !is_required, is_advanced_option);
    registerStringOption_("precursor_tolerance_units", "<choice>", "ppm", "tolerance_mass_units 0=ppm, 1=Da", !is_required, is_advanced_option);
    setValidStrings_("precursor_tolerance_units", ListUtils::create<String>("ppm,Da"));


  }

  // read and parse clustering output csv to store specnumber and clusterid associations
  void readMClusterOutputAsMap_(String mcout_file, std::map<MaRaClusterResult, Int>& specid_to_clusterid_map, const std::map<String, Int>& filename_to_idx_map)
  {
    CsvFile csv_file(mcout_file, '\t');
    StringList row;
    Int clusterid = 0;

    for (Size i = 0; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      if (!row.empty())
      {
        row[0] = String(filename_to_idx_map.at(row[0]));

        MaRaClusterResult res(row);
        specid_to_clusterid_map[res] = clusterid;
      }
      else
      {
        ++clusterid;
      }
    }
  }

  //   replace with PercolatorAdapter function
  String getScanIdentifier_(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
  {
    // MSGF+ uses this field, is empty if not specified
    String scan_identifier = it->getSpectrumReference();
    if (scan_identifier.empty())
    {
      // XTandem uses this (integer) field
      // these ids are 1-based in contrast to the index which is 0-based. This might be problematic to use for merging
      if (it->metaValueExists("spectrum_id") && !it->getMetaValue("spectrum_id").toString().empty())
      {
        scan_identifier = "scan=" + it->getMetaValue("spectrum_id").toString();
      }
      else
      {
        scan_identifier = "index=" + String(it - start + 1);
        OPENMS_LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
      }
    }
    return scan_identifier.removeWhitespaces();
  }

  //   replace with PercolatorAdapter function
  Int getScanNumber_(String scan_identifier)
  {
    Int scan_number = 0;
    StringList fields = ListUtils::create<String>(scan_identifier);
    for (const String& st : fields)
    {
      // if scan number is not available, use the scan index
      Size idx = 0;
      if ((idx = st.find("scan=")) != string::npos)
      {
        scan_number = st.substr(idx + 5).toInt();
        break;
      }
      else if ((idx = st.find("index=")) != string::npos)
      {
        scan_number = st.substr(idx + 6).toInt();
        break;
      }
      else if ((idx = st.find("spectrum=")) != string::npos)
      {
        scan_number = st.substr(idx + 9).toInt();
      }
    }
    return scan_number;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList in_list = getStringList_("in");

    const String maracluster_executable(getStringOption_("maracluster_executable"));
    writeDebug_(String("Path to the maracluster executable: ") + maracluster_executable, 2);

    String maracluster_output_directory = getStringOption_("output_directory");   
    const String consensus_out(getStringOption_("consensus_out"));
    const String out(getStringOption_("out"));

    if (in_list.empty())
    {
      writeLogError_("Error:  no input file given (parameter 'in')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (consensus_out.empty() && out.empty())
    {
      writeLogError_("Error:  no output file given (parameter 'out' or 'consensus_out')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------

    // create temp directory to store maracluster temporary files
    File::TempDir tmp_dir(debug_level_ >= 2);

    double pcut = getDoubleOption_("pcut");

    String txt_designator = File::getUniqueName();
    String input_file_list(tmp_dir.getPath() + txt_designator + ".file_list.txt");
    String consensus_output_file(tmp_dir.getPath() + txt_designator + ".clusters_p" + String(Int(-1*pcut)) + ".tsv");

    // Create simple text file with one file path per line
    // TODO make a bit more exception safe
    ofstream os(input_file_list.c_str());
    map<String,Int> filename_to_file_idx;
    Int file_idx = 0;
    for (StringList::const_iterator fit = in_list.begin(); fit != in_list.end(); ++fit, ++file_idx) {
      filename_to_file_idx[*fit] = file_idx;
      os << *fit;
      if (fit + 1 != in_list.end())
      {
        os << endl;
      }
    }
    os.close();

    QStringList arguments;
    // Check all set parameters and get them into arguments StringList
    {
      arguments << "batch";
      arguments << "-b" << input_file_list.toQString();
      arguments << "-f" << tmp_dir.getPath().toQString();
      arguments << "-a" << txt_designator.toQString();

      map<String,int> precursor_tolerance_units;
      precursor_tolerance_units["ppm"] = 0;
      precursor_tolerance_units["Da"] = 1;

      arguments << "-p" << (String(getDoubleOption_("precursor_tolerance")) + precursor_tolerance_units[getStringOption_("precursor_tolerance_units")]).toQString();

      arguments << "-t" << String(pcut).toQString();
      arguments << "-c" << String(pcut).toQString();

      Int verbose_level = getIntOption_("verbose");
      if (verbose_level != 2)
      {
        arguments << "-v" << String(verbose_level).toQString();
      }
    }
    writeLogInfo_("Prepared maracluster command.");

    //-------------------------------------------------------------
    // run MaRaCluster for idXML output
    //-------------------------------------------------------------
    // MaRaCluster execution with the executable and the arguments StringList
    writeLogInfo_("Executing maracluster ...");
    auto exit_code = runExternalProcess_(maracluster_executable.toQString(), arguments);
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // reintegrate clustering results 
    //-------------------------------------------------------------
    std::map<MaRaClusterResult, Int> specid_to_clusterid_map;
    readMClusterOutputAsMap_(consensus_output_file, specid_to_clusterid_map, filename_to_file_idx);
    file_idx = 0;

    //if specified keep original output in designated directory
    if (!maracluster_output_directory.empty())
    {
      bool copy_status = File::copyDirRecursively(tmp_dir.getPath().toQString(), maracluster_output_directory.toQString());

      if (copy_status)
      { 
        OPENMS_LOG_INFO << "MaRaCluster original output was successfully copied to " << maracluster_output_directory << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "MaRaCluster original output could not be copied to " << maracluster_output_directory << ". Please run MaRaClusterAdapter with debug >= 2." << std::endl;
      }
    }

    //output idXML containing scannumber and cluster id annotation
    if (!out.empty())
    {
      const StringList id_in = getStringList_("id_in");
      vector<PeptideIdentification> all_peptide_ids;
      vector<ProteinIdentification> all_protein_ids;
      if (!id_in.empty())
      {
        for (const String& ss : id_in) {
          vector<PeptideIdentification> peptide_ids;
          vector<ProteinIdentification> protein_ids;
          FileHandler().loadIdentifications(ss, protein_ids, peptide_ids, {FileTypes::IDXML});
          for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it) {
            String scan_identifier = getScanIdentifier_(it, peptide_ids.begin());
            Int scan_number = getScanNumber_(scan_identifier);
            MaRaClusterResult res(file_idx, scan_number);
            // cluster index - 1 is equal to scan_number in consensus.mzML
            Int cluster_id = specid_to_clusterid_map[res] - 1;
            it->setMetaValue("cluster_id", cluster_id);
            String filename = in_list[file_idx];
            it->setMetaValue("file_origin", filename);
          }
          for (ProteinIdentification& prot : protein_ids) {
            String filename = in_list[file_idx];
            prot.setMetaValue("file_origin", filename);
          }
          all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
          all_protein_ids.insert(all_protein_ids.end(), protein_ids.begin(), protein_ids.end());
        }
        ++file_idx;
      }
      else
      {
        for (std::map<MaRaClusterResult,Int>::iterator sid = specid_to_clusterid_map.begin(); sid != specid_to_clusterid_map.end(); ++sid) {
          Int scan_nr = sid->first.scan_nr;
          Int file_id = sid->first.file_idx;
          Int cluster_id = sid->second;
          PeptideIdentification pid;
          PeptideHit pih;
          pid.insertHit(pih);
          pid.setSpectrumReference("scan=" + String(scan_nr));
          // cluster index - 1 is equal to scan_number in consensus.mzML
          pid.setMetaValue("cluster_id", cluster_id - 1);
          pid.setMetaValue("file_origin", in_list[file_id]);
          all_peptide_ids.push_back(pid);
        }
      }

      if (all_protein_ids.empty())
      {
        ProteinIdentification protid;
        all_protein_ids.push_back(protid);
      }

      all_protein_ids.back().setMetaValue("maracluster", "MaRaClusterAdapter");
      ProteinIdentification::SearchParameters search_parameters = all_protein_ids.back().getSearchParameters();

      search_parameters.setMetaValue("MaRaCluster:pvalue-cutoff", pcut);

      all_protein_ids.back().setSearchParameters(search_parameters);


      writeDebug_("write idXMLFile", 1);
      writeDebug_(out, 1);// As the maracluster output file is not needed anymore, the temporary directory is going to be deleted
      FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {FileTypes::IDXML});
    }

    //output consensus mzML
    if (!consensus_out.empty())
    {
      QStringList arguments_consensus;
      // Check all set parameters and get them into arguments StringList
      {
        arguments_consensus << "consensus";
        arguments_consensus << "-l" << consensus_output_file.toQString();
        arguments_consensus << "-f" << tmp_dir.getPath().toQString();
        arguments_consensus << "-o" << consensus_out.toQString();
        Int min_cluster_size = getIntOption_("min_cluster_size");
        arguments_consensus << "-M" << String(min_cluster_size).toQString();

        Int verbose_level = getIntOption_("verbose");
        if (verbose_level != 2) arguments_consensus << "-v" << String(verbose_level).toQString();
      }
      writeLogInfo_("Prepared maracluster-consensus command.");

      //-------------------------------------------------------------
      // run MaRaCluster for consensus output
      //-------------------------------------------------------------
      // MaRaCluster execution with the executable and the arguments StringList
      TOPPBase::ExitCodes exit_code = runExternalProcess_(maracluster_executable.toQString(), arguments_consensus);
      if (exit_code != EXECUTION_OK)
      {
        return exit_code;
      }

      // sort mzML
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(consensus_out);
      OPENMS_LOG_DEBUG << "Input type" << FileTypes::typeToName(in_type) << ". " << std::endl;
      PeakMap exp;
      fh.loadExperiment(FileHandler::stripExtension(consensus_out) + ".part1." + FileTypes::typeToName(in_type), exp, {in_type}, log_type_, true, true);
      exp.sortSpectra();
      fh.storeExperiment(consensus_out, exp, {FileTypes::MZML}, log_type_);
    }

    writeLogInfo_("MaRaClusterAdapter finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  MaRaClusterAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
