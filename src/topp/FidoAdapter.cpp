// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_FidoAdapter FidoAdapter

    @brief Runs the protein inference engine Fido.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FidoAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability\n(with "prob_correct" option) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
    </table>
</CENTER>

    In order to prepare suitable input data for Fido, peptide/protein identification results (e.g. from a database search engine) should be processed with @ref TOPP_PeptideIndexer (with the "prob_correct" option) followed by @ref TOPP_IDPosteriorErrorProbability (with the "annotate_proteins" option).

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
    TOPPBase("FidoAdapter", "Runs the protein inference engine Fido.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output: identification results "
                        "with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("exe", "<file>", 
#if defined(OPENMS_WINDOWSPLATFORM)
                          "FidoChooseParameters.exe",
#else
                          "FidoChooseParameters",
#endif
                          "Executable for Fido with parameter estimation", 
                          false);
  }


  ExitCodes main_(int, const char**)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String fido_executable = getStringOption_("fido_executable");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    IdXMLFile().load(in, proteins, peptides);
    if (proteins.empty() || peptides.empty())
    {
      LOG_ERROR << "Error: Input file '" << in 
                << "' should contain both protein and peptide data." << endl;
      return INPUT_FILE_EMPTY;
    }
    if (proteins.size() > 1)
    {
      LOG_WARN << "Warning: Input contains more than one protein "
               << "identification run. Protein inference results will be "
               << "written to the first run only." << endl;
    }

    String temp_directory = 
      QDir::toNativeSeparators((File::getTempDirectory() + "/" + 
                                File::getUniqueName() + "/").toQString());
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }

    String fido_input_graph = temp_directory + "fido_input_graph.dat";
    String fido_input_proteins = temp_directory + "fido_input_proteins.dat";

    // write PSM graph:
    ofstream graph_out(fido_input_graph.c_str());
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      pep_it->sort();
      const PeptideHit& hit = pep_it->getHits()[0];
      if (hit.getSequence().empty() || 
          hit.getProteinAccessions().empty()) continue;

      double score = hit.getScore();

      String error_reason;
      if (!pep_it->isHigherScoreBetter())
      {
        error_reason = "lower scores are better";
      }
      else if (score < 0.0)
      {
        error_reason = "score < 0";
      }
      else if (score > 1.0)
      {
        error_reason = "score > 1";
      }
      if (!error_reason.empty())
      {
        String msg = "Error: Unsuitable score type for peptide-spectrum "
          "matches detected (problem: " + error_reason + ").\nFido requires "
          "probabilities as scores, e.g. as produced by "
          "IDPosteriorErrorProbability with the 'prob_correct' option.";
        LOG_ERROR << msg << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }

      graph_out << "e " << hit.getSequence() << endl;
      for (vector<String>::const_iterator acc_it = 
             hit.getProteinAccessions().begin(); acc_it != 
             hit.getProteinAccessions().end(); ++acc_it)
      {
        graph_out << "r " << *acc_it << endl;
      }
      graph_out << "p " << score << endl;
    }
    graph_out.close();

    // gather protein target/decoy data:
    set<String> targets, decoys;
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      for (vector<ProteinHit>::iterator hit_it = prot_it->getHits().begin();
           hit_it != prot_it->getHits().end(); ++hit_it)
      {
        String target_decoy = hit_it->getMetaValue("target_decoy").toString();
        if (target_decoy == "target")
        {
          targets.insert(hit_it->getAccession());
        }
        else if (target_decoy == "decoy")
        {
          decoys.insert(hit_it->getAccession());
        }
        else
        {
          String msg = "Error: All protein hits must be annotated with target/"
            "decoy meta data. Run PeptideIndexer with the 'annotate_proteins' "
            "option to accomplish this.";
          LOG_ERROR << msg << endl;
          return INCOMPATIBLE_INPUT_DATA;
        }
      }
    }

    // write target/decoy protein sets:
    ofstream proteins_out(fido_input_proteins.c_str());
    proteins_out << "{ ";
    for (set<String>::iterator it = targets.begin(); it != targets.end(); ++it)
    {
      if (it != targets.begin()) proteins_out << " , ";
      proteins_out << *it;
    }
    proteins_out << " }\n{ ";
    for (set<String>::iterator it = decoys.begin(); it != decoys.end(); ++it)
    {
      if (it != decoys.begin()) proteins_out << " , ";
      proteins_out << *it;
    }
    proteins_out << " }" << endl;
    proteins_out.close();
    
    // run program and read output:
    QProcess fido;
    QStringList inputs;
    inputs << fido_input_graph.toQString() << fido_input_proteins.toQString();
    fido.start(fido_executable.toQString(), inputs);

    ExitCodes exit_code = EXECUTION_OK;
    if (!fido.waitForFinished(-1))
    {
      String msg = "Fatal error running Fido (command: '" + fido_executable +
        " \"" + String(inputs.join("\" \"")) + "\"').\n"
        "Does the Fido executable exist?";
      LOG_ERROR << msg << endl;
      exit_code = EXTERNAL_PROGRAM_ERROR;
    }
    else // success!
    {
      String output = QString(fido.readAllStandardOutput());
      vector<String> lines;
      output.split("\n", lines);

      vector<ProteinIdentification::ProteinGroup> groups;
      for (vector<String>::iterator line_it = lines.begin(); 
           line_it != lines.end(); ++line_it)
      {
        // format of a line (example):
        // 0.6788 { SW:TRP6_HUMAN , GP:AJ271067_1 , GP:AJ271068_1 }
        istringstream line(*line_it);
        ProteinIdentification::ProteinGroup group;
        line >> group.probability;
        String accessions;
        line >> accessions;
        vector<String> parts;
        accessions.split(" ", parts);
        for (vector<String>::iterator part_it = parts.begin(); 
             part_it != parts.end(); ++part_it)
        {
          if (part_it->size() > 1) // skip braces and commas
          {
            group.accessions.push_back(*part_it);
          }
        }
      }
      proteins[0].getIndistinguishableProteins() = groups;

      // write output:
      IdXMLFile().store(out, proteins, peptides);
    }

    // clean up temporary files
    if (debug_level_ < 2)
    {
      File::removeDirRecursively(temp_directory);
      LOG_INFO << "Set debug level to 2 or higher to keep temporary files at '"
               << temp_directory << "'" << endl;
    }
    else
    {
      LOG_INFO << "Keeping temporary files at '" << temp_directory
               << "'. Set debug level to below 2 to remove them." << endl;
    }

    return exit_code;
  }

};


int main(int argc, const char** argv)
{
  TOPPFidoAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
