// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

    @note FalseDiscoveryRate only annotates peptides and proteins with their FDR. A subsequent FDR filtering step needs to be conducted downstream via @ref IDFilter or after exporting the data to a text-based format.

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
    registerInputFile_("in", "<file>", "", "Identifications from searching a target-decoy database.");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identifications with annotated FDR");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerFlag_("proteins_only", "If set only the FDR on protein level is calculated");
    registerFlag_("peptides_only", "If set only the FDR on peptide (PSM) level is calculated");

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
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    bool proteins_only = getFlag_("proteins_only");
    bool peptides_only = getFlag_("peptides_only");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

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
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFalseDiscoveryRate tool;
  return tool.main(argc, argv);
}

/// @endcond
