// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <OpenMS/QC/DBSuitability.h>
#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DatabaseSuitability DatabaseSuitability

@brief Calculates the suitability of a database which was used a for peptide identification search. Also reports the quality of LC-MS spectra.

@dot
digraph sample_workflow {
  node [ style="solid,filled", color=black, fillcolor=grey90, width=1.5, fixedsize=true, shape=square, fontname=Helvetica, fontsize=10 ];
  edge [ arrowhead="open", style="solid" ];
  rankdir="LR";
  splines=ortho;
  mzml [ label="mzML file(s)" shape=oval fillcolor=white group=1];
  novor [ label="NovorAdapter" URL="\ref OpenMS::NovorAdapter" group=2];
  id_filter [ label="IDFilter" URL="\ref OpenMS::IDFilter" group=2];
  id_convert [ label="IDFileConverter" URL="\ref OpenMS::IDFileConverter" group=2];
  decoy_db [ label="DecoyDatabase" URL="\ref OpenMS::DecoyDatabase" group=2];
  comet [ label="CometAdapter" URL="\ref OpenMS::CometAdapter" group=1];
  pep_ind [ label="PeptideIndexer" URL="\ref OpenMS::PeptideIndexer" group=1];
  db_suit [ label="DatabaseSuitability" fillcolor="#6F42C1" fontcolor=white group=3];
  tsv [ label="optional\ntsv output" shape=oval fillcolor=white group=3];
  {rank = same; db_suit; decoy_db;}
  mzml -> novor;
  mzml -> comet;
  comet -> pep_ind;
  pep_ind -> db_suit [ xlabel="in_id" ];
  novor -> id_filter;
  id_filter -> id_convert;
  id_convert -> decoy_db;
  decoy_db -> comet;
  mzml -> db_suit [ xlabel="in_spec" ];
  novor -> db_suit [ xlabel="in_novor" ];
  db_suit -> tsv;
}
@enddot

The metric this tool uses to determine the suitability of a database is based on a de novo model. Therefore it is crucial that your workflow is set up the right way. Above you can see an example.@n
Most importantly the peptide identification search needs to be done with a combination of the database in question and a de novo "database".@n
To generate the de novo "database":
      - @ref UTILS_NovorAdapter calculates de novo sequences.
      - @ref TOPP_IDFilter can filter out unwanted ones.
      - @ref TOPP_IDFileConverter generates the de novo fasta file.

For re-ranking all cases where a peptide hit only found in the de novo "database" scores above a peptide hit found in the actual database are checked. In all these cases the cross-correlation scores of those peptide hits are compared. If they are similar enough, the database hit will be re-ranked to be on top of the de novo hit. You can control how much of cases with similar scores will be re-ranked by using the @p reranking_cutoff_percentile.@n
For this to work it is important @ref TOPP_PeptideIndexer ran before. However it is also crucial that no FDR was performed. This tool does this itself and will crash if a q-value is found. You can still control the FDR that you want to establish using the corresponding flag.

@note For identification search the only supported search engine for the time being is Comet because the Comet cross-correlation score is needed for re-ranking.@n
You can still uses other search engines and disable the re-ranking via the @p no_rerank flag in this tool. This will probably result in an underestimated suitability though.@n


The results are written directly into the console. But you can provide an optional tsv output file where the most important results will be exported to.


This tool uses the metrics and algorithms first presented in:@n
<em>Assessing protein sequence database suitability using de novo sequencing. Molecular & Cellular Proteomics. January 1, 2020; 19, 1: 198-208. doi:10.1074/mcp.TIR119.001752.@n
Richard S. Johnson, Brian C. Searle, Brook L. Nunn, Jason M. Gilmore, Molly Phillips, Chris T. Amemiya, Michelle Heck, Michael J. MacCoss.</em>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DatabaseSuitability.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DatabaseSuitability.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

Citation c = { "Richard S. Johnson, Brian C. Searle, Brook L. Nunn, Jason M. Gilmore, Molly Phillips, Chris T. Amemiya, Michelle Heck, Michael J. MacCoss",
                    "Assessing protein sequence database suitability using de novo sequencing",
                    "Molecular & Cellular Proteomics. January 1, 2020; 19, 1: 198-208",
                    "10.1074/mcp.TIR119.001752" };

class DatabaseSuitability :
  public TOPPBase
{
public:
  DatabaseSuitability() :
    TOPPBase("DatabaseSuitability", "Computes a suitability score for a database which was used for a peptide identification search. Also reports the quality of LC-MS spectra.", false, {c})
  {
  }
protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return DBSuitability().getDefaults();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_id", "<file>", "", "Input idXML file from peptide search with combined database with added de novo peptide. PeptideIndexer is needed, FDR is forbidden.");
    setValidFormats_("in_id", { "idXML" });
    registerInputFile_("in_spec", "<file>", "", "Input MzML file used for the peptide identification");
    setValidFormats_("in_spec", { "mzML" });
    registerInputFile_("in_novo", "<file>", "", "Input idXML file containing de novo peptides (unfiltered)");
    setValidFormats_("in_novo", { "idXML" });
    registerOutputFile_("out", "<file>", "", "Optional tsv output containing database suitability information as well as spectral quality.", false);
    setValidFormats_("out", { "tsv" });

    registerSubsection_("algorithm", "Parameter section for the suitability calculation algorithm");
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in_id = getStringOption_("in_id");
    String in_spec = getStringOption_("in_spec");
    String in_novo = getStringOption_("in_novo");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    
    MzMLFile m;
    PeakFileOptions op;
    op.setMSLevels({ 2 }); // only ms2
    m.setOptions(op);
    PeakMap exp;
    m.load(in_spec, exp);
    
    IdXMLFile x;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    x.load(in_id, prot_ids, pep_ids);

    vector<ProteinIdentification> novo_prots;
    vector<PeptideIdentification> novo_peps;
    x.load(in_novo, novo_prots, novo_peps);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    set<AASequence> unique_novo;
    for (const auto& pep_id : novo_peps)
    {
      if (pep_id.getScoreType() != "novorscore")
      {
        OPENMS_LOG_ERROR << in_novo << " contains at least one identification without a novorscore! Make sure this file contains only deNovo sequences." << endl;
        return INPUT_FILE_CORRUPT;
      }
      if (pep_id.getHits().empty()) continue;
      unique_novo.insert(pep_id.getHits()[0].getSequence());
    }

    Ms2IdentificationRate q;
    q.compute(novo_peps, exp, true);

    Ms2IdentificationRate::IdentificationRateData spectral_quality = q.getResults()[0];

    DBSuitability s;
    Param p = getParam_().copy("algorithm:", true);
    s.setParameters(p);
    s.compute(pep_ids);

    DBSuitability::SuitabilityData suit = s.getResults()[0];

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << suit.num_top_db << " / " << (suit.num_top_db + suit.num_top_novo) << " top hits were found in the database." << endl;
    OPENMS_LOG_INFO << suit.num_top_novo << " / " << (suit.num_top_db + suit.num_top_novo) << " top hits were only found in the concatenated de novo peptide." << endl;
    OPENMS_LOG_INFO << suit.num_interest << " times scored a de novo hit above a database hit. Of those times " << suit.num_re_ranked << " top de novo hits where re-ranked." << endl;
    OPENMS_LOG_INFO << "database suitability [0, 1]: " << suit.suitability << endl << endl;
    OPENMS_LOG_INFO << unique_novo.size() << " / " << spectral_quality.num_peptide_identification << " de novo sequences are unique" << endl;
    OPENMS_LOG_INFO << spectral_quality.num_ms2_spectra << " ms2 spectra found" << endl;
    OPENMS_LOG_INFO << "spectral quality (id rate of de novo sequences) [0, 1]: " << spectral_quality.identification_rate << endl << endl;

    if (!out.empty())
    {
      OPENMS_LOG_INFO << "Writing output to: " << out << endl << endl;

      std::ofstream os(out);
      if (!os.is_open())
      {
        OPENMS_LOG_ERROR << "Output file given in 'out' isn't writable." << endl;
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      os.precision(writtenDigits(double()));
      os << "key\tvalue\n";
      os << "#top_db_hits\t" << suit.num_top_db << "\n";
      os << "#top_novo_hits\t" << suit.num_top_novo << "\n";
      os << "db_suitability\t" << suit.suitability << "\n";
      os << "#total_novo_seqs\t" << spectral_quality.num_peptide_identification << "\n";
      os << "#unique_novo_seqs\t" << unique_novo.size() << "\n";
      os << "#ms2_spectra\t" << spectral_quality.num_ms2_spectra << "\n";
      os << "spectral_quality\t" << spectral_quality.identification_rate << "\n";
      os.close();
    }

    return EXECUTION_OK;
  }
};


// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  DatabaseSuitability tool;
  return tool.main(argc, argv);
}

/// @endcond
