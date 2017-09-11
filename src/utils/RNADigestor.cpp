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
   @page UTILS_RNADigestor RNADigestor

    @brief Digests an RNA sequence database in-silico.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ RNADigestor \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (so far)</td>
        </tr>
    </table>
</CENTER>

    This application is used to digest an RNA sequence database to get all fragments given a cleavage enzyme.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RNADigestor.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RNADigestor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRNADigestor :
  public TOPPBase
{
public:
  TOPPRNADigestor() :
    TOPPBase("RNADigestor", "Digests an RNA sequence database in-silico.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
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
  }

  ExitCodes main_(int, const char**)
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
    set<String> unique_fragments;

    for (vector<FASTAFile::FASTAEntry>::const_iterator fa_it = seq_data.begin();
         fa_it != seq_data.end(); ++fa_it)
    {
      vector<String> fragments;
      digestor.digest(fa_it->sequence, fragments, min_size, max_size);
      Size counter = 1;
      for (vector<String>::const_iterator frag_it = fragments.begin();
           frag_it != fragments.end(); ++frag_it)
      {
        if (!unique || !unique_fragments.count(*frag_it))
        {
          String id = fa_it->identifier + "_" + String(counter);
          String desc;
          if (!fa_it->description.empty()) desc = fa_it->description + " ";
          desc += "(fragment " + String(counter) + ")";
          FASTAFile::FASTAEntry fragment(id, desc, *frag_it);
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

    LOG_INFO << "Digested " << seq_data.size() << " sequence(s) into "
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
