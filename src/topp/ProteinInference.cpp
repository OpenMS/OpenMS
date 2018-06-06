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

    @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

    This tool counts the peptide sequences that match a protein accession. From this count for all protein hits in the respective id run, only those proteins are accepted that have at least a given number of peptides sequences identified. The peptide identifications should be prefiltered with respect to false discovery rate and the score in general to remove bad identifications.

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
    TOPPBase("ProteinInference", "Protein inference based on the number of identified peptides.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    addEmptyLine_();
    registerIntOption_("min_peptides_per_protein", "<num>", 2, "Minimal number of peptides needed for a protein identification", false);
    setMinInt_("min_peptides_per_protein", 1);

    registerFlag_("treat_charge_variants_separately", "If this flag is set, different charge variants of the same peptide sequence count as inidividual evidences.");
    registerFlag_("treat_modification_variants_separately", "If this flag is set, different modification variants of the same peptide sequence count as individual evidences.");
    //registerSubsection_("algorithm","Consensus algorithm section");
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    Size min_peptides_per_protein = getIntOption_("min_peptides_per_protein");
    bool treat_charge_variants_separately(getFlag_("treat_charge_variants_separately"));
    bool treat_modification_variants_separately(getFlag_("treat_modification_variants_separately"));

    // load identifications
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile().load(in, prot_ids, pep_ids);

    // collect the different proteins (some of the protein hit copies are discarded)
    Map<String, ProteinHit> acc_to_protein_hit;
    for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        acc_to_protein_hit[pit->getAccession()] = *pit;
      }
    }

    writeDebug_(String(acc_to_protein_hit.size()) + " different protein accessions in the file.", 1);


    // count the sequences that match a protein accession
    // ProtAcc --> [charge, PepSeq]
    Map<String, Map<Size, set<String> > > acc_peptides;
    for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
      // for all peptide hits
      for (vector<PeptideHit>::const_iterator it2 = it1->getHits().begin(); it2 != it1->getHits().end(); ++it2)
      {
        String pep_seq;
        if (treat_modification_variants_separately)
        {
          pep_seq = it2->getSequence().toString();
        }
        else
        {
          pep_seq = it2->getSequence().toUnmodifiedString();
        }
        Size charge = 0;
        if (treat_charge_variants_separately)
        {
          charge = it2->getCharge();
        }

        // for all protein accessions
        set<String> protein_accessions = it2->extractProteinAccessionsSet();
        for (set<String>::const_iterator it3 = protein_accessions.begin(); it3 != protein_accessions.end(); ++it3)
        {
          acc_peptides[*it3][charge].insert(pep_seq);
        }
      }
    }

    writeDebug_("Peptides from " + String(acc_peptides.size()) + " proteins recorded.", 1);

    // for all protein hits for the id run, only accept proteins that have at least 'min_peptides_per_protein' peptides
    set<String> accepted_proteins;
    vector<ProteinHit> accepted_protein_hits;
    for (Map<String, ProteinHit>::ConstIterator it1 = acc_to_protein_hit.begin(); it1 != acc_to_protein_hit.end(); ++it1)
    {
      if (acc_peptides.has(it1->first))
      {
        Size num_peps(0);
        for (Map<Size, set<String> >::ConstIterator it2 = acc_peptides[it1->first].begin(); it2 != acc_peptides[it1->first].end(); ++it2)
        {
          num_peps += it2->second.size();
        }

        if (num_peps >= min_peptides_per_protein)
        {
          accepted_proteins.insert(it1->first);
          accepted_protein_hits.push_back(it1->second);
        }
      }
    }

    writeDebug_("Accepted " + String(accepted_protein_hits.size()) + " proteins.", 1);
    writeDebug_("Accepted " + String(accepted_proteins.size()) + " proteins.", 1);

    // remove peptides that are not accepted
    for (vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
      vector<PeptideHit> peptide_hits = it1->getHits();
      it1->setHits(vector<PeptideHit>());
      for (vector<PeptideHit>::const_iterator it2 = peptide_hits.begin(); it2 != peptide_hits.end(); ++it2)
      {
        set<String> protein_accessions = it2->extractProteinAccessionsSet();
        for (set<String>::const_iterator it3 = protein_accessions.begin(); it3 != protein_accessions.end(); ++it3)
        {
          if (accepted_proteins.find(*it3) != accepted_proteins.end())
          {
            it1->insertHit(*it2);
            break;
          }
        }
      }
    }

    // remove proteins that are not accepted
    prot_ids.resize(1);
    prot_ids[0].setHits(accepted_protein_hits);

    // fix wrong accessions of the peptides (to proteins that were removed)
    for (vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
      vector<PeptideHit> peptide_ids = it1->getHits();
      for (vector<PeptideHit>::iterator it2 = peptide_ids.begin(); it2 != peptide_ids.end(); ++it2)
      {
        vector<PeptideEvidence> filtered_evidence;
        vector<PeptideEvidence> old_evidence = it2->getPeptideEvidences();

        for (vector<PeptideEvidence>::const_iterator evidence_it = old_evidence.begin(); evidence_it != old_evidence.end(); ++evidence_it)
        {
          if (accepted_proteins.find(evidence_it->getProteinAccession()) != accepted_proteins.end())
          {
            filtered_evidence.push_back(*evidence_it);
          }
        }
        it2->setPeptideEvidences(filtered_evidence);
      }
      it1->setHits(peptide_ids);
    }

    DateTime now = DateTime::now();
    String identifier(now.get() + "_TOPPProteinInference");
    for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
    {
      it->setIdentifier(identifier);
    }

    prot_ids[0].setIdentifier(identifier);

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
