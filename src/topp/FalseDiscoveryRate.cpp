// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
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

    This TOPP tool calculates the false discovery rate (FDR) for results of target-decoy searches. It can handle separate searches of a target database (e.g. forward sequences) and a decoy database (e.g. reversed sequences), or a combined search of a concatenated target-decoy database. The FDR calculation can be performed for proteins and/or for peptides (more exactly, peptide spectrum matches).

    The false discovery rate is defined as the number of false discoveries (decoy hits) divided by the number of false and correct discoveries (both target and decoy hits) with a score better than a given threshold.

    When using a combined database of target and decoy sequences (thus only running one search per ID engine), @ref TOPP_PeptideIndexer must be applied to the search results (idXML file) to index the data and to annotate peptide and protein hits with their target/decoy status.

    @note When no decoy hits were found you will get a warning like this:<br>
    "FalseDiscoveryRate: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0!"<br>
    This should be a serious concern, since it indicates a possible problem with the target/decoy annotation step (@ref TOPP_PeptideIndexer), e.g. due to a misconfigured database.

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

  Param getSubsectionDefaults_(const String& /*section*/) const
  {
    return FalseDiscoveryRate().getDefaults();
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Identification input file containing a search against "
                                           "a concatenated sequence database. "
                                           "Either specify '-in' alone or 'fwd_in' together with 'rev_in' as input.", false);
    setValidFormats_("in", ListUtils::create<String>("idXML"));

    registerInputFile_("fwd_in", "<file>", "", "Identification input to estimate FDR, forward run.", false);
    setValidFormats_("fwd_in", ListUtils::create<String>("idXML"));
    registerInputFile_("rev_in", "<file>", "", "Identification input to estimate FDR, decoy run.", false);
    setValidFormats_("rev_in", ListUtils::create<String>("idXML"));

    registerOutputFile_("out", "<file>", "", "Identification output with annotated FDR");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerFlag_("proteins_only", "If set, the FDR of the proteins only is calculated");
    registerFlag_("peptides_only", "If set, the FDR of the peptides only is calculated");

    registerSubsection_("algorithm", "Parameter section for the FDR calculation algorithm");
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    Param alg_param = getParam_().copy("algorithm:", true);
    FalseDiscoveryRate fdr;

    if (!alg_param.empty())
    {
      fdr.setParameters(alg_param);
      writeDebug_("Parameters passed to FalseDiscoveryRate", alg_param, 3);
    }

    // input/output files
    // either "fwd_in" and "rev_in" must be given, or just "in" (which contains results of a search against a concatenated target-decoy sequence db):
    String fwd_in = getStringOption_("fwd_in"),
      rev_in = getStringOption_("rev_in"), in = getStringOption_("in");
    bool combined = false;
    if (!fwd_in.empty() && !rev_in.empty())
    {
      if (!in.empty())
      {
        writeLog_("Error, either 'fwd_in' and 'rev_in' must be given or 'in', but not both");
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      if (!in.empty())
      {
        combined = true;
      }
      else
      {
        writeLog_("Error, at least 'fwd_in' and 'rev_in' or 'in' must be given");
        return ILLEGAL_PARAMETERS;
      }
    }
    String out = getStringOption_("out");
    bool proteins_only = getFlag_("proteins_only");
    bool peptides_only = getFlag_("peptides_only");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    if (combined) // -in was given
    {
      vector<PeptideIdentification> pep_ids;
      vector<ProteinIdentification> prot_ids;

      IdXMLFile().load(in, prot_ids, pep_ids);

      try
      {
        if (!proteins_only)
        {
          fdr.apply(pep_ids);
        }
        if (!peptides_only)
        {
          fdr.apply(prot_ids);
        }
      }
      catch (Exception::MissingInformation)
      {
        LOG_FATAL_ERROR << "FalseDiscoveryRate failed due to missing information (see above).\n";
        return INCOMPATIBLE_INPUT_DATA;
      }

      for (vector<ProteinIdentification>::iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
      {
        it->assignRanks();
      }
      for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
      {
        it->assignRanks();
      }

      IdXMLFile().store(out, prot_ids, pep_ids);
    }
    else         // -fw_in & rev_in given
    {
      vector<PeptideIdentification> fwd_pep, rev_pep;
      vector<ProteinIdentification> fwd_prot, rev_prot;

      IdXMLFile().load(fwd_in, fwd_prot, fwd_pep);
      IdXMLFile().load(rev_in, rev_prot, rev_pep);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      writeDebug_("Starting calculations with " + String(fwd_pep.size()) + "/" + String(rev_pep.size()) + " read peptide IDs", 1);

      if (!proteins_only)
      {
        fdr.apply(fwd_pep, rev_pep);
      }
      if (!peptides_only)
      {
        fdr.apply(fwd_prot, rev_prot);
      }

      // TODO @all shouldn't ranks be assigned here as well?

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      IdXMLFile().store(out, fwd_prot, fwd_pep);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFalseDiscoveryRate tool;
  return tool.main(argc, argv);
}

/// @endcond
