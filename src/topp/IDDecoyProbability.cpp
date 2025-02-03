// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/ID/IDDecoyProbability.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_IDDecoyProbability IDDecoyProbability

@brief Util to estimate probability of peptide hits

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; IDDecoyProbability &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
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
The probabilities are calculated using Bayes law, similar to PeptideProphet.
This implementation is much simpler than that of PeptideProphet.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDDecoyProbability.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDDecoyProbability.html

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
    TOPPBase("IDDecoyProbability", "Estimates peptide probabilities using a decoy search strategy.\nWARNING: This util is deprecated.")
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
    if (!fwd_in.empty() && !rev_in.empty())
    {
      if (!in.empty())
      {
        writeLogError_("Error: either 'fwd_in' and 'rev_in' must be given or 'in', but not both");
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
        writeLogError_("Error: at least 'fwd_in' and 'rev_in' or 'in' must be given");
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
      FileHandler().loadIdentifications(fwd_in, fwd_prot, fwd_pep, {FileTypes::IDXML});
      FileHandler().loadIdentifications(rev_in, rev_prot, rev_pep, {FileTypes::IDXML});

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

      writeDebug_("Starting calculations", 1);
      decoy_prob.apply(out_pep, fwd_pep, rev_pep);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      FileHandler().storeIdentifications(out, fwd_prot, out_pep, {FileTypes::IDXML});
    }
    else
    {
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      String document_id;
      FileHandler().loadIdentifications(in, prot_ids, pep_ids, {FileTypes::IDXML});

      decoy_prob.apply(pep_ids);
      FileHandler().storeIdentifications(out, prot_ids, pep_ids, {FileTypes::IDXML});
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

