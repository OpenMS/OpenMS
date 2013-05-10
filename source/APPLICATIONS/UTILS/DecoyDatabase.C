// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_DecoyDatabase DecoyDatabase

    @brief Create a decoy peptide database from standard FASTA databases.

  Decoy databases are useful to control false discovery rates and thus estimate score cutoffs for identified spectra.

  The decoy can either be generated from reversed or shuffled sequences.

  To get a 'contaminants' database have a look at http://www.thegpm.org/crap/index.html or find/create your own contaminant database.

  Multiple databases can be provided as input, which will internally be concatenated before being used for decoy generation.
  This allows you to specify your target database plus a contaminant file and (upon using the @p append flag) obtain a concatenated
  target-decoy database using a single call, e.g., DecoyDatabase -in human.fasta crap.fasta -out human_TD.fasta -append


    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_DecoyDatabase.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_DecoyDatabase.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDecoyDatabase :
  public TOPPBase
{
public:
  TOPPDecoyDatabase() :
    TOPPBase("DecoyDatabase", "Create decoy peptide databases from normal ones.", false)
  {

  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<file(s)>", StringList::create(""), "Input FASTA file(s), each containing a database. It is recommended to include a contaminant database as well.");
    setValidFormats_("in", StringList::create("fasta"));
    registerOutputFile_("out", "<file>", "", "Output FASTA file where the decoy database will be written to.");
    setValidFormats_("out", StringList::create("fasta"));
    registerStringOption_("decoy_string", "<string>", "_rev", "String that is appended to the accession of the protein database to indicate a decoy protein.", false);
    registerStringOption_("decoy_string_position", "<enum>", "suffix", "Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?", false);
    setValidStrings_("decoy_string_position", StringList::create("prefix,suffix"));
    registerFlag_("append", "If this flag is used, the decoy database is appended to the target database, allowing combined target decoy searches.");
    registerFlag_("shuffle", "If 'true' then the decoy hit are shuffled from the target sequences, otherwise they are reversed");
  }

  String getIdentifier_(const String & identifier, const String & decoy_string, const bool as_prefix)
  {
    if (as_prefix) return decoy_string + identifier;
    else return identifier + decoy_string;
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in(getStringList_("in"));
    String out(getStringOption_("out"));
    bool append = getFlag_("append");
    bool shuffle = getFlag_("shuffle");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    vector<FASTAFile::FASTAEntry> proteins;
    for (Size i=0; i<in.size(); ++i)
    {
      vector<FASTAFile::FASTAEntry> single_proteins;
      // this is a little inefficient since it requires copying; appending during load() would be better...
      FASTAFile().load(in[i], single_proteins);  
      // append
      proteins.insert(proteins.end(), single_proteins.begin(), single_proteins.end());
    }

    if (in.size() == 1)
    {
      LOG_WARN << "Warning: Only one FASTA input file was provided, which might not contain contaminants. You probably want to have them! Just add the contaminant file to the input file list 'in'." << endl;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    String decoy_string(getStringOption_("decoy_string"));
    bool decoy_string_position_prefix =   (String(getStringOption_("decoy_string_position")) == "prefix" ? true : false);
    Size num_proteins = proteins.size();
    set<String> identifiers;
    if (shuffle)
    {
      for (Size i = 0; i < num_proteins; ++i)
      {
        if (identifiers.find(proteins[i].identifier) != identifiers.end())
        {
          LOG_WARN << "DecoyDatabase: Warning, identifier is not unique to sequence file: '" << proteins[i].identifier << "'!" << endl;
        }
        identifiers.insert(proteins[i].identifier);

        FASTAFile::FASTAEntry entry = proteins[i];

        String pro_seq, temp;
        pro_seq = entry.sequence;
        Size x = pro_seq.size();
        srand(time(0));
        while (x != 0)
        {
          Size y = rand() % x;
          temp += pro_seq[y];
          pro_seq[y] = pro_seq[x - 1];
          --x;
        }
        entry.sequence = temp;

        if (append)
        {
          entry.identifier = getIdentifier_(entry.identifier, decoy_string, decoy_string_position_prefix);
          proteins.push_back(entry);
        }
        else
        {
          proteins[i].sequence = entry.sequence;
          proteins[i].identifier = getIdentifier_(proteins[i].identifier, decoy_string, decoy_string_position_prefix);
        }
      }
    }
    else         // !shuffle
    {
      for (Size i = 0; i < num_proteins; ++i)
      {
        if (identifiers.find(proteins[i].identifier) != identifiers.end())
        {
          LOG_WARN << "DecoyDatabase: Warning, identifier is not unique to sequence file: '" << proteins[i].identifier << "'!" << endl;
        }
        identifiers.insert(proteins[i].identifier);

        if (append)
        {
          FASTAFile::FASTAEntry entry = proteins[i];
          entry.sequence.reverse();
          entry.identifier = getIdentifier_(entry.identifier, decoy_string, decoy_string_position_prefix);
          proteins.push_back(entry);
        }
        else
        {
          proteins[i].sequence.reverse();
          proteins[i].identifier = getIdentifier_(proteins[i].identifier, decoy_string, decoy_string_position_prefix);
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    FASTAFile().store(out, proteins);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPDecoyDatabase tool;
  return tool.main(argc, argv);
}

/// @endcond
