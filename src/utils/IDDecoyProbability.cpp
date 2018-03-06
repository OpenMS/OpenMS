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
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

/**
    @page UTILS_IDDecoyProbability IDDecoyProbability

    @brief Util to estimate probability of peptide hits

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDDecoyProbability \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
        </tr>
    </table>
</CENTER>

    @experimental This util is deprecated and might behave not as expected!

    So far an estimation of the false score distribution with a gamma distribution
    and the correct score distribution with a gaussian distribution is performed.
    The probabilities are calculated using bayes law, similar to PeptideProphet.
    This implementation is much simpler than that of PeptideProphet.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDDecoyProbability.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDDecoyProbability.html

    For the parameters of the algorithm section see the algorithms documentation: @n
        @ref OpenMS::IDDecoyProbability "decoy_algorithm" @n

*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDDecoyProbability :
  public TOPPBase
{
public:
  TOPPIDDecoyProbability() :
    TOPPBase("IDDecoyProbability", "Estimates peptide probabilities using a decoy search strategy.\nWARNING: This util is deprecated.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Identification input of combined forward decoy search (reindex with PeptideIndexer first)", false);
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerInputFile_("fwd_in", "<file>", "", "Identification input of forward run", false);
    setValidFormats_("fwd_in", ListUtils::create<String>("idXML"));
    registerInputFile_("rev_in", "<file>", "", "Identification input of decoy run", false);
    setValidFormats_("rev_in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Identification output with forward scores converted to probabilities");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerSubsection_("decoy_algorithm", "Algorithm parameter subsection");
    addEmptyLine_();
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    IDDecoyProbability decoy_prob;
    return decoy_prob.getParameters();
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    // either fwd_in and rev_in must be given or just the in which contains results of a search against a concatenated target decoy sequence db
    String fwd_in(getStringOption_("fwd_in")), rev_in(getStringOption_("rev_in")), in(getStringOption_("in"));
    bool combined(false);
    if (fwd_in != "" && rev_in != "")
    {
      if (in != "")
      {
        writeLog_("Error, either 'fwd_in' and 'rev_in' must be given or 'in', but not both");
        return ILLEGAL_PARAMETERS;
      }
    }
    else
    {
      if (in != "")
      {
        combined = true;
      }
      else
      {
        writeLog_("Error, at least 'fwd_in' and 'rev_in' or 'in' must be given");
        return ILLEGAL_PARAMETERS;
      }
    }

    String out(getStringOption_("out"));

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    IDDecoyProbability decoy_prob;
    Param decoy_param = getParam_().copy("decoy_algorithm:", true);
    decoy_prob.setParameters(decoy_param);

    if (!combined)
    {
      vector<PeptideIdentification> fwd_pep, rev_pep, out_pep;
      vector<ProteinIdentification> fwd_prot, rev_prot;
      String document_id;
      IdXMLFile().load(fwd_in, fwd_prot, fwd_pep, document_id);
      IdXMLFile().load(rev_in, rev_prot, rev_pep, document_id);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      writeDebug_("Starting calculations", 1);
      decoy_prob.apply(out_pep, fwd_pep, rev_pep);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      IdXMLFile().store(out, fwd_prot, out_pep);
    }
    else
    {
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      String document_id;
      IdXMLFile().load(in, prot_ids, pep_ids, document_id);

      decoy_prob.apply(pep_ids);
      IdXMLFile().store(out, prot_ids, pep_ids);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPIDDecoyProbability tool;
  return tool.main(argc, argv);
}

/// @endcond

