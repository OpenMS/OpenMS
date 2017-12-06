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
  This allows you to specify your target database plus a contaminant file and obtain a concatenated
  target-decoy database using a single call, e.g., DecoyDatabase -in human.fasta crap.fasta -out human_TD.fasta

  By default, a combined database is created where target and decoy sequences are written interleaved
  (i.e., target1, decoy1, target2, decoy2,...).
  If you need all targets before the decoys for some reason, use @p only_decoy and concatenate the files
  externally.

  The tool will keep track of all protein identifiers and report duplicates.

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
    TOPPBase("DecoyDatabase", "Create decoy protein DB from forward protein DB.", false)
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file(s)>", ListUtils::create<String>(""), "Input FASTA file(s), each containing a database. It is recommended to include a contaminant database as well.");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "Output FASTA file where the decoy database will be written to.");
    setValidFormats_("out", ListUtils::create<String>("fasta"));
    registerStringOption_("decoy_string", "<string>", "DECOY_", "String that is combined with the accession of the protein identifier to indicate a decoy protein.", false);
    registerStringOption_("decoy_string_position", "<enum>", "prefix", "Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?", false);
    setValidStrings_("decoy_string_position", ListUtils::create<String>("prefix,suffix"));
    registerFlag_("only_decoy", "Write only decoy proteins to the output database instead of a combined database.", false);
    registerStringOption_("method", "<enum>", "reverse", "Method by which decoy sequences are generated from target sequences.", false);
    setValidStrings_("method", ListUtils::create<String>("reverse,shuffle"));
  }

  String getIdentifier_(const String & identifier, const String & decoy_string, const bool as_prefix)
  {
    if (as_prefix) return decoy_string + identifier;
    else return identifier + decoy_string;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in(getStringList_("in"));
    String out(getStringOption_("out"));
    bool append = (!getFlag_("only_decoy"));
    bool shuffle = (getStringOption_("method") == "shuffle");
    String decoy_string(getStringOption_("decoy_string"));
    bool decoy_string_position_prefix = (String(getStringOption_("decoy_string_position")) == "prefix" ? true : false);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    if (in.size() == 1)
    {
      LOG_WARN << "Warning: Only one FASTA input file was provided, which might not contain contaminants. You probably want to have them! Just add the contaminant file to the input file list 'in'." << endl;
    }

    set<String> identifiers; // spot duplicate identifiers  // std::unordered_set<string> has slightly more RAM, but slightly less CPU

    FASTAFile f;
    f.writeStart(out);
    FASTAFile::FASTAEntry protein;
      
    for (Size i = 0; i < in.size(); ++i)
    {
      f.readStart(in[i]);  

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
      while (f.readNext(protein))
      {
        if (identifiers.find(protein.identifier) != identifiers.end())
        {
          LOG_WARN << "DecoyDatabase: Warning, identifier '" << protein.identifier << "' occurs more than once!" << endl;
        }
        identifiers.insert(protein.identifier);

        if (append)
        {
          f.writeNext(protein);
        }
      
        // identifier
        protein.identifier = getIdentifier_(protein.identifier, decoy_string, decoy_string_position_prefix);
      
        // sequence
        if (shuffle)
        {
          String temp;
          Size x = protein.sequence.size();
          srand(time(nullptr));
          while (x != 0)
          {
            Size y = rand() % x;
            temp += protein.sequence[y];
            --x;
            protein.sequence[y] = protein.sequence[x]; // overwrite consumed position with last position (about to go out of scope for next dice roll)
          }
        }
        else // reverse
        {
          protein.sequence.reverse();
        }
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        f.writeNext(protein);
      
      } // next protein
    } // input files

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPDecoyDatabase tool;
  return tool.main(argc, argv);
}

/// @endcond
