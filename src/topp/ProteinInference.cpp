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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <set>
#include <unordered_set>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_ProteinInference ProteinInference

    @brief Computes a protein identification based on the number of identified peptides.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ProteinInterference \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines)</td>
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

    This tool counts and aggregates the scores of peptide sequences that match a protein accession.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_ProteinInference.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_ProteinInference.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPProteinInference :
  public TOPPBase
{
public:
  TOPPProteinInference() :
    TOPPBase("ProteinInference", "Protein inference based on an aggregation of the scores of the identified peptides.")
  {
  }

  enum AggregationMethod {PROD, SUM, MAXIMUM};

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    addEmptyLine_();
    //registerIntOption_("min_peptides_per_protein", "<num>", 2, "Minimal number of peptides needed for a protein identification", false);
    //setMinInt_("min_peptides_per_protein", 1);

    registerFlag_("treat_charge_variants_separately", "If this flag is set, different charge variants of the same peptide sequence count as inidividual evidences.");
    registerFlag_("treat_modification_variants_separately", "If this flag is set, different modification variants of the same peptide sequence count as individual evidences.");
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    //Size min_peptides_per_protein = getIntOption_("min_peptides_per_protein");
    bool treat_charge_variants_separately(getFlag_("treat_charge_variants_separately"));
    bool treat_modification_variants_separately(getFlag_("treat_modification_variants_separately"));
    bool use_shared_peptides(getFlag_("use_shared_peptides"));

    AggregationMethod aggregation_method = AggregationMethod::MAXIMUM;

    String aggMethodString = "";

    switch (aggregation_method)
    {
      case AggregationMethod::PROD :
        aggMethodString = "Product";
        break;
      case AggregationMethod::SUM :
        aggMethodString = "Sum";
        break;
      case AggregationMethod::MAXIMUM :
        aggMethodString = "Maximum";
        break;
      default:
        break;
    }


    // load identifications
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile().load(in, prot_ids, pep_ids);

    unordered_map<std::string, map<Int,PeptideHit*>> best_pep{};

    // iterate over runs
    unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    for (auto& prot_run : prot_ids)
    {
      acc_to_protein_hitP_and_count.clear();
      best_pep.clear();

      ProteinIdentification::SearchParameters sp = prot_run.getSearchParameters();
      prot_run.setSearchEngine("TOPPProteinInference_" + aggMethodString);
      sp.setMetaValue("use_shared_peptides", use_shared_peptides);
      sp.setMetaValue("treat_charge_variants_separately", treat_charge_variants_separately);
      sp.setMetaValue("treat_modification_variants_separately", treat_modification_variants_separately);

      //create Accession to ProteinHit and peptide count map. To have quick access later.
      for (auto& phit : prot_run.getHits())
      {
        acc_to_protein_hitP_and_count[phit.getAccession()] = std::make_pair<ProteinHit*, Size>(&phit, 0);
      }

      for (auto& pep : pep_ids)
      {
        //skip if it does not belong to run
        if (pep.getIdentifier() != prot_run.getIdentifier()) continue;

        //skip if shared and skip_shared
        if (!use_shared_peptides && pep.metaValueExists("protein_references") && pep.getMetaValue("protein_references") == "non-unique") continue;

        for (auto& hit : pep.getHits())
        {
          String lookup_seq;
          if (!treat_modification_variants_separately)
          {
            lookup_seq = hit.getSequence().toUnmodifiedString();
          }
          else
          {
            lookup_seq = hit.getSequence().toString();
          }

          int lookup_charge = 0;
          if (treat_charge_variants_separately)
          {
            lookup_charge = hit.getCharge();
          }

          auto current_best_pep_it = best_pep.find(lookup_seq);
          if (current_best_pep_it == best_pep.end())
          {
            auto& new_entry = best_pep[lookup_seq];
            new_entry = std::map<Int, PeptideHit*>();
            new_entry[lookup_charge] = &hit;
          }
          else
          {
            auto current_best_pep_charge_it = current_best_pep_it->second.find(lookup_charge);
            if (current_best_pep_charge_it == current_best_pep_it->second.end())
            {
              current_best_pep_it->second[lookup_charge] = &hit;
            }
            else if ((prot_run.isHigherScoreBetter() && (hit.getScore() > current_best_pep_charge_it->second->getScore())) ||
                     (!prot_run.isHigherScoreBetter() && (hit.getScore() < current_best_pep_charge_it->second->getScore())))
            {
              current_best_pep_charge_it->second = &hit;
            }
          }
        }
      }

      // update protein scores
      for (const auto &charge_to_pep_hit_map : best_pep)
      {
        // The next line assumes that PeptideHits of different charge states necessarily share the same
        // protein accessions
        // TODO this could be done for mods later, too (first hashing AASeq, then the mods)
        for (const auto &acc : charge_to_pep_hit_map.second.begin()->second->extractProteinAccessionsSet())
        {
          for (const auto &pep_hit : charge_to_pep_hit_map.second)
          {
            auto prot_count_pair = acc_to_protein_hitP_and_count[acc];
            ProteinHit* protein = prot_count_pair.first;
            prot_count_pair.second++;

            double new_score = pep_hit.second->getScore();

            // Note: This requires/works only with Posterior (Error) Probabilities
            if (!prot_run.isHigherScoreBetter()) new_score = 1. - new_score;
            switch (aggregation_method)
            {
              //TODO for 0 probability peptides we could also multiply a minimum value
              case AggregationMethod::PROD :
                if (new_score > 0.0) protein->setScore(protein->getScore() * new_score);
                break;
              case AggregationMethod::SUM :
                protein->setScore(protein->getScore() + new_score);
                break;
              case AggregationMethod::MAXIMUM :
                protein->setScore(max(double(protein->getScore()), new_score));
                break;
              default:
                break;
            }
          }
        }
      }
    }

    //TODO Filtering? I think this should be done separate afterwards with IDFilter
    //for all protein hits for the id run, only accept proteins that have at least 'min_peptides_per_protein' peptides

    // write output
    IdXMLFile().store(out, prot_ids, pep_ids);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPProteinInference tool;
  return tool.main(argc, argv);
}

/// @endcond
