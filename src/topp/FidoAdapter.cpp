// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_FidoAdapter FidoAdapter

    @brief Computes a protein identification based on the target and decoy peptide IDs.
    <CENTER>
      <table>
      <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ FidoAdapter \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
      </tr>
      <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter (or other ID engines)</td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_PeptideIndexer </td>
      </tr>
      <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
      </tr>
      <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
      </tr>
      </table>
    </CENTER>
    @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

    Fido should be installed properly before running this tool.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_FidoAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FidoAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFidoAdapter :
  public TOPPBase
{
public:
  TOPPFidoAdapter() :
    TOPPBase("FidoAdapter", "A protein inference tool.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    addEmptyLine_();
    registerInputFile_("fido_executable", "<executable>", "/fido/bin/Fido", "bin directory e.g. '/usr/local/fido/bin/Fido' or '/usr/local/fido/bin/FidoChooseParameters'", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("default_input_file", "<file>", "", "Default parameters input file, if not given default parameters are used", false);
    addEmptyLine_();
    registerStringOption_("peptide_probability", "<string>", "", "Specify the MetaValue of peptide hits that will be used as the peptide probability, otherwise the score will be used", false);
    registerFlag_("fido_choose_parameter_off", "Disable FidoChooseParameters program. If set, three parameters: alpha, beta, gamma must be given mannually.");
    registerDoubleOption_("alpha", "<float>", 0, "Fido paramter, needed when not running FidoChooseParameter.", false);
    registerDoubleOption_("beta", "<float>", 0, "Fido paramter, needed when not running FidoChooseParameter.", false);
    registerDoubleOption_("gamma", "<float>", 0, "Fido paramter, needed when not running FidoChooseParameter.", false);
    registerStringOption_("decoy_string", "<string>", "_rev", "String that was appended (or prepended - see 'prefix' flag below) to the accession of the protein database to indicate a decoy protein.", false);
    registerFlag_("clean_peptide_names", "FidoChooseParameters option -p: omits cleaning the peptide names.");
    registerFlag_("use_all_psms", "FidoChooseParameters option -a: uses all PSM matches instead the best one.");
    registerFlag_("use_group_level", "FidoChooseParameters option -g: uses protein group level inference.");
    registerIntOption_("accurary_level", "<Int>", 1, "FidoChooseParameters option -c: sets start parameter's accurary level (1-3): 1 = best    / slower; 2 = relaxed / faster; 3 = sloppy  / very fast", false);
    registerFlag_("log2_size", "FidoChooseParameters option l1: l1 is the log2 of maximum number of subgraph connected states.");
    registerFlag_("log2_calculation", "FidoChooseParameters option l2: l2 is the log2 for the main calculation, and l1 is only used.");
  }

  // Fido psm_graph_file gamma alpha beta -options
  void writeFidoinput_(String& graph_filename, vector<PeptideIdentification>& pep_ids, String& peptide_probability)
  {
    /* psm graph file:
    e EEEMPEPK
    r SW:TRP6_HUMAN
    r GP:AJ271067_1
    r GP:AJ271068_1
    p 0.9849
    */
    ofstream graph(graph_filename.c_str());
    if (!graph)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, graph_filename);
    }

    for (vector<PeptideIdentification>::iterator id = pep_ids.begin(); id != pep_ids.end(); ++id)
    {
      //higher_score_better="false";
      for (vector<PeptideHit>::iterator hit_it = id->getHits().begin(); hit_it != id->getHits().end(); ++hit_it)
      {
        if (hit_it->getProteinAccessions().size() > 0)
        {
          graph << "e " << hit_it->getSequence() << "\n";
          StringList proteins = hit_it->getProteinAccessions();
          for (Size i = 0; i < hit_it->getProteinAccessions().size(); ++i)
          {
            graph << "r " << hit_it->getProteinAccessions()[i] << "\n";
          }
        }
        graph << "p ";
        if (peptide_probability.empty())
        {
          if (id->isHigherScoreBetter())
          {
            graph << hit_it->getScore() << "\n";
          }
          else
          {
            graph << 1.0 - hit_it->getScore() << "\n";
          }
        }
        else
        {
          graph << hit_it->getMetaValue(peptide_probability) << "\n";
        }
      }
    }
    graph.close();
  }

  void loadFidoOutput_(StringList& protein_score, vector<ProteinIdentification>& prot_ids, vector<PeptideIdentification>& pep_ids)
  {
    //1 { P16083|NQO2_HUMAN }
    //0.9097144245 { Q86XS8|GOLI_HUMAN }
    //0.9097144245 { Q9NVH1|DJC11_HUMAN , P05129|KPCG_HUMAN }
    StringList proteins_list;
    map<String, String> fido_results;
    for (Size i = 0; i < protein_score.size(); ++i)
    {
      StringList sp = ListUtils::create<String>(protein_score[i], ' ');
      for (Size j = 2; j < sp.size(); j = j + 2)
      {
        String protein_name = sp[j];
        proteins_list.push_back(protein_name);
        fido_results.insert(pair<String, String>(protein_name, sp[0]));
        if (sp[j + 1] == "}") break;
      }
    }

    for (vector<ProteinIdentification>::iterator pid = prot_ids.begin(); pid != prot_ids.end(); ++pid)
    {
      ProteinIdentification pit = *pid;
      pid->setScoreType("fido");
      vector<ProteinHit> hits;
      for (vector<ProteinHit>::iterator hit_it = pit.getHits().begin();
           hit_it != pit.getHits().end(); ++hit_it)
      {
        String protein_accession = hit_it->getAccession();
        if (fido_results[protein_accession].c_str())
        {
          double score = atof(fido_results[protein_accession].c_str());
          hit_it->setScore(score);
          hits.push_back(*hit_it);
        }
      }
      pid->setHits(hits);
      pid->assignRanks();

    }
    for (vector<PeptideIdentification>::iterator id = pep_ids.begin(); id != pep_ids.end(); ++id)
    {
      vector<PeptideHit> pep_hits;
      for (Size i = 0; i < proteins_list.size(); ++i)
      {
        vector<PeptideHit> peptides;
        id->getReferencingHits(proteins_list[i], peptides);
        for (Size j = 0; j < peptides.size(); ++j) pep_hits.push_back(peptides[j]);
      }
      id->setHits(pep_hits);
    }
  }

  // FidoChooseParameters -options psm_graph_file target_decoy_file
  void writeTargetdecoyfile_(String& targetdecoy_filename, vector<ProteinIdentification>& prot_ids, String& decoy_string)
  {
    vector<String> target, decoy;
    for (vector<ProteinIdentification>::iterator pid = prot_ids.begin(); pid != prot_ids.end(); ++pid)
    {
      ProteinIdentification pit = *pid;
      for (vector<ProteinHit>::const_iterator hit_it = pit.getHits().begin();
           hit_it != pit.getHits().end(); ++hit_it)
      {
        String protein_accession = hit_it->getAccession();
        if (protein_accession.hasSubstring(decoy_string))
        {
          decoy.push_back(protein_accession);
        }
        else
        {
          target.push_back(protein_accession);
        }
      }
    }
    ofstream f(targetdecoy_filename.c_str());
    if (!f)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, targetdecoy_filename);
    }
    f << "{ " << target[0];
    for (Size i = 1; i < target.size(); ++i) f << " , " << target[i];
    f << " }\n";
    f << "{ " << decoy[0];
    for (Size i = 1; i < decoy.size(); ++i) f << " , " << decoy[i];
    f << " }\n";
    f.close();
  }

  ExitCodes main_(int, const char**)
  {
    // path to the log file
    String logfile(getStringOption_("log"));
    String inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name == "")
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile().load(inputfile_name, prot_ids, pep_ids);

    String outputfile_name = getStringOption_("out");
    writeDebug_(String("Output file: ") + outputfile_name, 1);
    if (outputfile_name == "")
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String temp_directory = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }
    writeLog_("Write Fido psm graph file...");
    //-------------------------------------------------------------
    // write Fido graph input files
    //-------------------------------------------------------------
    String psm_graph_file(temp_directory + "psm_graph_file");
    String peptide_probability(getStringOption_("peptide_probability"));
    writeFidoinput_(psm_graph_file, pep_ids, peptide_probability);

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList paramlist;
    String exe_path;
    String fido_executable;
    if (File::isDirectory(getStringOption_("fido_executable")))
    {
      exe_path = getStringOption_("fido_executable");
    }
    else
    {
      exe_path = File::path(getStringOption_("fido_executable"));
    }
    if (getFlag_("fido_choose_parameter_off")) // only run fido: Fido psm_graph_file 0.9 0.04 0 > output is stdout
    {
      fido_executable = exe_path + "/Fido";
      paramlist << psm_graph_file;
      paramlist << String(getDoubleOption_("gamma"));
      paramlist << String(getDoubleOption_("alpha"));
      paramlist << String(getDoubleOption_("beta"));
    }
    else // run FidoChooseParameters -p -a -c 1 psm_graph_file target_decoy_file > stdout
    {
      writeLog_("Write Fido target_decoy_file...");
      String targetdecoy_filename(temp_directory + "targetdecoy_file");
      String decoy_string(getStringOption_("decoy_string"));
      writeTargetdecoyfile_(targetdecoy_filename, prot_ids, decoy_string);


      fido_executable = exe_path + "/FidoChooseParameters";
      if (getFlag_("clean_peptide_names"))
      {
        paramlist << "-p";
      }
      if (getFlag_("use_all_psms"))
      {
        paramlist << "-a";
      }
      if (getFlag_("use_group_level"))
      {
        paramlist << "-g";
      }
      paramlist << "-c " + String(getIntOption_("accurary_level"));
      paramlist << psm_graph_file;
      paramlist << targetdecoy_filename;
      writeDebug_("Run FidoChooseParameters...", 5);
    }

    //-------------------------------------------------------------
    // run Fido, Fido output is always stdout
    //-------------------------------------------------------------
    QStringList qparam;
    for (Size i = 0; i < paramlist.size(); ++i)
    {
      qparam << paramlist[i].toQString();
    }

    QProcess program;
    program.start(fido_executable.toQString(), qparam);

    bool success = program.waitForFinished();
    String fido_output(QString(program.readAllStandardOutput()));
    String fido_parameters(QString(program.readAllStandardError()));
    if (!success || program.exitStatus() != 0 || program.exitCode() != 0)
    {
      writeLog_("Fido problem. Aborting! Calling command was: '" + fido_executable + " \"" + inputfile_name + "\"'.\nDoes the Fido executable exist?");
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
      }
      return EXTERNAL_PROGRAM_ERROR;
    }
    else
    {
      if (!fido_parameters.empty())
      {
        writeLog_("Fido choose parameters process...");
        writeLog_(fido_parameters);
      }
      StringList protein_scores = ListUtils::create<String>(fido_output, '\n');
      writeLog_("Load Fido results...");
      loadFidoOutput_(protein_scores, prot_ids, pep_ids);
    }

    //-------------------------------------------------------------
    // writing results into idXML
    //-------------------------------------------------------------
    IdXMLFile().store(outputfile_name, prot_ids, pep_ids);

    // Deletion of temporary files
    if (this->debug_level_ < 2)
    {
      File::removeDirRecursively(temp_directory);
      LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
    }
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFidoAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
