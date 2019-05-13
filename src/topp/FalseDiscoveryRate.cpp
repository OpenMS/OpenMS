// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Timo Sachsenberg, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

/**
    @page TOPP_FalseDiscoveryRate FalseDiscoveryRate

    @brief Tool to estimate the false discovery rate on peptide and protein level
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FalseDiscoveryRate \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDFilter </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
        </tr>
    </table>
</CENTER>

    This TOPP tool calculates the false discovery rate (FDR) for results of target-decoy searches. The FDR calculation can be performed for proteins and/or for peptides (more exactly, peptide spectrum matches).

    The false discovery rate is defined as the number of false discoveries (decoy hits) divided by the number of false and correct discoveries (both target and decoy hits) with a score better than a given threshold.

    @ref TOPP_PeptideIndexer must be applied to the search results (idXML file) to index the data and to annotate peptide and protein hits with their target/decoy status.

    @note When no decoy hits were found you will get a warning like this:<br>
    "FalseDiscoveryRate: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0!"<br>
    This should be a serious concern, since it indicates a possible problem with the target/decoy annotation step (@ref TOPP_PeptideIndexer), e.g. due to a misconfigured database.

    @note FalseDiscoveryRate only annotates peptides and proteins with their FDR. By setting FDR:PSM or FDR:protein the maximum q-value (e.g., 0.05 corresponds to an FDR of 5%) can be controlled on the PSM and protein level.
    Alternativly, FDR filtering can be performed in the @ref IDFilter tool by setting score:pep and score:prot to the maximum q-value.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_FalseDiscoveryRate.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FalseDiscoveryRate.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFalseDiscoveryRate :
  public TOPPBase
{
public:
  TOPPFalseDiscoveryRate() :
    TOPPBase("FalseDiscoveryRate", "Estimates the false discovery rate on peptide and protein level using decoy searches.")
  {
  }

protected:

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return FalseDiscoveryRate().getDefaults();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Identifications from searching a target-decoy database.");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identifications with annotated FDR");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("PSM", "<FDR level>", "true", "Perform FDR calculation on PSM level", false);
    setValidStrings_("PSM", ListUtils::create<String>("true,false"));
    registerStringOption_("protein", "<FDR level>", "true", "Perform FDR calculation on protein level", false);
    setValidStrings_("protein", ListUtils::create<String>("true,false"));

    registerTOPPSubsection_("FDR", "FDR control");
    registerDoubleOption_("FDR:PSM", "<fraction>", 1, "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)", false);
    setMinFloat_("FDR:PSM", 0);
    setMaxFloat_("FDR:PSM", 1);
    registerDoubleOption_("FDR:protein", "<fraction>", 1, "Filter proteins based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)", false);
    setMinFloat_("FDR:protein", 0);
    setMaxFloat_("FDR:protein", 1);
    registerSubsection_("algorithm", "Parameter section for the FDR calculation algorithm");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    Param alg_param = getParam_().copy("algorithm:", true);
    FalseDiscoveryRate fdr;

    fdr.setParameters(alg_param);
    writeDebug_("Parameters passed to FalseDiscoveryRate", alg_param, 3);

    // input/output files
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const double protein_fdr = getDoubleOption_("FDR:protein");
    const double psm_fdr = getDoubleOption_("FDR:PSM");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    vector<PeptideIdentification> pep_ids;
    vector<ProteinIdentification> prot_ids;

    IdXMLFile().load(in, prot_ids, pep_ids);

    Size n_prot_ids = prot_ids.size();
    Size n_prot_hits = IDFilter::countHits(prot_ids);
    Size n_pep_ids = pep_ids.size();
    Size n_pep_hits = IDFilter::countHits(pep_ids);

    bool filter_applied(false);

    try
    {
      if (getStringOption_("protein") == "true")
      {
        fdr.apply(prot_ids);
        filter_applied = true;

        if (protein_fdr < 1)
        {
          OPENMS_LOG_INFO << "FDR control: Filtering proteins..." << endl;
          IDFilter::filterHitsByScore(prot_ids, protein_fdr);
        }
      }

      if (getStringOption_("PSM") == "true")
      {
        fdr.apply(pep_ids);
        filter_applied = true;

        if (psm_fdr < 1)
        {
          OPENMS_LOG_INFO << "FDR control: Filtering PSMs..." << endl;
          IDFilter::filterHitsByScore(pep_ids, psm_fdr);
        }
      }
    }
    catch (Exception::MissingInformation&)
    {
      OPENMS_LOG_FATAL_ERROR << "FalseDiscoveryRate failed due to missing information (see above).\n";
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (filter_applied)
    {
      IDFilter::removeUnreferencedProteins(prot_ids, pep_ids);

      // keep decoy peptide hits without decoy protein references if flag is specified
      if (alg_param.getValue("add_decoy_peptides").toBool() == true)
      {
        IDFilter::updateProteinReferences(pep_ids, prot_ids, false);
      }
      else
      {
        IDFilter::updateProteinReferences(pep_ids, prot_ids, true);
      }    
      IDFilter::updateHitRanks(prot_ids);
      IDFilter::updateHitRanks(pep_ids);
      IDFilter::removeEmptyIdentifications(pep_ids);
      // we want to keep "empty" protein IDs because they contain search meta data
    }

    // update protein groupings if necessary:
    for (auto prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
    {
      bool valid = IDFilter::updateProteinGroups(prot_it->getProteinGroups(),
                                                 prot_it->getHits());
      if (!valid)
      {
        OPENMS_LOG_WARN << "Warning: While updating protein groups, some prot_ids were removed from groups that are still present. "
                 << "The new grouping (especially the group probabilities) may not be completely valid any more." 
                 << endl;
      }

      valid = IDFilter::updateProteinGroups(
        prot_it->getIndistinguishableProteins(), prot_it->getHits());

      if (!valid)
      {
        OPENMS_LOG_WARN << "Warning: While updating indistinguishable prot_ids, some prot_ids were removed from groups that are still present. "
                 << "The new grouping (especially the group probabilities) may not be completely valid any more." 
                 << endl;
      }
    }

    // some stats
    OPENMS_LOG_INFO << "Before filtering:\n"
             << n_prot_ids << " protein identification(s) with "
             << n_prot_hits << " protein hit(s),\n"
             << n_pep_ids << " peptide identification(s) with "
             << n_pep_hits << " pep_ids hit(s).\n"
             << "After filtering:\n"
             << prot_ids.size() << " protein identification(s) with "
             << IDFilter::countHits(prot_ids) << " protein hit(s),\n"
             << pep_ids.size() << " peptide identification(s) with "
             << IDFilter::countHits(pep_ids) << " pep_ids hit(s)." << endl;

    OPENMS_LOG_INFO << "Writing filtered output..." << endl;
    IdXMLFile().store(out, prot_ids, pep_ids);
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFalseDiscoveryRate tool;
  return tool.main(argc, argv);
}

/// @endcond
