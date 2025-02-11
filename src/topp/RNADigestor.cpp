// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Nico Pfeiffer, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/FORMAT/FASTAFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_RNADigestor RNADigestor

@brief Digests an RNA sequence database in-silico.
<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; RNADigestor &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (so far)</td>
        </tr>
    </table>
</CENTER>

This application is used to digest an RNA sequence database to get all fragments given a cleavage enzyme.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_RNADigestor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_RNADigestor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRNADigestor :
  public TOPPBase
{
public:
  TOPPRNADigestor() :
    TOPPBase("RNADigestor", "Digests an RNA sequence database in-silico.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file containing RNA sequences");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "Output file containing sequence fragments");
    setValidFormats_("out", ListUtils::create<String>("fasta"));

    registerIntOption_("missed_cleavages", "<number>", 1, "The number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<number>", 3, "Minimum length of a fragment", false);
    registerIntOption_("max_length", "<number>", 30, "Maximum length of a fragment", false);
    vector<String> all_enzymes;
    RNaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<string>", "RNase_T1", "Digestion enzyme (RNase)", false);
    setValidStrings_("enzyme", all_enzymes);
    registerFlag_("unique", "Report each unique sequence fragment only once");
    registerFlag_("cdna", "Input file contains cDNA sequences - replace 'T' with 'U')");
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    Size min_size = getIntOption_("min_length");
    Size max_size = getIntOption_("max_length");
    Size missed_cleavages = getIntOption_("missed_cleavages");

    bool unique = getFlag_("unique");
    bool cdna = getFlag_("cdna");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    vector<FASTAFile::FASTAEntry> seq_data;
    FASTAFile().load(in, seq_data);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String enzyme = getStringOption_("enzyme");
    RNaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(missed_cleavages);

    std::vector<FASTAFile::FASTAEntry> all_fragments;
    set<NASequence> unique_fragments;

    for (FASTAFile::FASTAEntry& entry : seq_data)
    {
      vector<NASequence> fragments;
      if (cdna) entry.sequence.toUpper().substitute('T', 'U');
      NASequence seq = NASequence::fromString(entry.sequence);
      digestor.digest(seq, fragments, min_size, max_size);
      Size counter = 1;
      for (vector<NASequence>::const_iterator frag_it = fragments.begin();
           frag_it != fragments.end(); ++frag_it)
      {
        if (!unique || !unique_fragments.count(*frag_it))
        {
          String id = entry.identifier + "_" + String(counter);
          String desc;
          if (!entry.description.empty()) desc = entry.description + " ";
          desc += "(fragment " + String(counter) + ")";
          FASTAFile::FASTAEntry fragment(id, desc, frag_it->toString());
          all_fragments.push_back(fragment);
          unique_fragments.insert(*frag_it);
          counter++;
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    FASTAFile().store(out, all_fragments);

    OPENMS_LOG_INFO << "Digested " << seq_data.size() << " sequence(s) into "
             << all_fragments.size()
             << " fragments meeting the length restrictions." << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPRNADigestor tool;
  return tool.main(argc, argv);
}

/// @endcond
