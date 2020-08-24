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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <boost/regex.hpp>

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
    TOPPBase("DecoyDatabase", "Create decoy sequence database from forward sequence database.", false)
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
    registerStringOption_("decoy_string_position", "<choice>", "prefix", "Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?", false);
    setValidStrings_("decoy_string_position", ListUtils::create<String>("prefix,suffix"));
    registerFlag_("only_decoy", "Write only decoy proteins to the output database instead of a combined database.", false);

    registerStringOption_("type", "<choice>", "protein", "Type of sequence. RNA sequences may contain modification codes, which will be handled correctly if this is set to 'RNA'.", false);
    setValidStrings_("type", ListUtils::create<String>("protein,RNA"));

    registerStringOption_("method", "<choice>", "reverse", "Method by which decoy sequences are generated from target sequences. Note that all sequences are shuffled using the same random seed, ensuring that identical sequences produce the same shuffled decoy sequences. Shuffled sequences that produce highly similar output sequences are shuffled again (see shuffle_sequence_identity_threshold).", false);
    setValidStrings_("method", ListUtils::create<String>("reverse,shuffle"));
    registerIntOption_("shuffle_max_attempts", "<int>", 30, "shuffle: maximum attempts to lower the amino acid sequence identity between target and decoy for the shuffle algorithm", false, true);
    registerDoubleOption_("shuffle_sequence_identity_threshold", "<double>", 0.5, "shuffle: target-decoy amino acid sequence identity threshold for the shuffle algorithm. If the sequence identity is above this threshold, shuffling is repeated. In case of repeated failure, individual amino acids are 'mutated' to produce a different amino acid sequence.", false, true);

    registerStringOption_("seed", "<int>", '1', "Random number seed (use 'time' for system time)", false, true);

    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<enzyme>", "Trypsin", "Enzyme used for the digestion of the sample. Only applicable if parameter 'type' is 'protein'.",false);
    setValidStrings_("enzyme", all_enzymes);

    registerSubsection_("Decoy", "Decoy parameters section");
  }

  Param getSubsectionDefaults_(const String& /* name */) const override
  {
    Param p = MRMDecoy().getDefaults();
    // change the default to also work with other proteases
    p.setValue("non_shuffle_pattern", "", "Residues to not shuffle (keep at a constant position when shuffling). Separate by comma, e.g. use 'K,P,R' here.");
    return p;
  }

  String getIdentifier_(const String& identifier, const String& decoy_string, const bool as_prefix)
  {
    if (as_prefix) return decoy_string + identifier;
    else return identifier + decoy_string;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    enum SeqType {protein, RNA};
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    bool append = !getFlag_("only_decoy");
    bool shuffle = (getStringOption_("method") == "shuffle");
    String decoy_string = getStringOption_("decoy_string");
    bool decoy_string_position_prefix =
      (getStringOption_("decoy_string_position") == "prefix");
    SeqType input_type = SeqType::protein; //default to protein
    if (getStringOption_("type") == "RNA")
    {
      input_type = SeqType::RNA;
    }

    Param decoy_param = getParam_().copy("Decoy:", true);
    bool keepN = decoy_param.getValue("keepPeptideNTerm").toBool();
    bool keepC = decoy_param.getValue("keepPeptideCTerm").toBool();

    String keep_const_pattern = decoy_param.getValue("non_shuffle_pattern");
    Int max_attempts = getIntOption_("shuffle_max_attempts");
    double identity_threshold = getDoubleOption_("shuffle_sequence_identity_threshold");

    // Set the seed for shuffling always to the same number, this
    // ensures that identical peptides get shuffled the same way
    // every time (without keeping track of them explicitly). This
    // will ensure that the total number of unique tryptic peptides
    // is identical in both databases.
    int seed;
    String seed_option(getStringOption_("seed"));
    if (seed_option == "time") seed = time(nullptr);
    else seed = seed_option.toInt();

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    if (in.size() == 1)
    {
      OPENMS_LOG_WARN << "Warning: Only one FASTA input file was provided, which might not contain contaminants. "
               << "You probably want to have them! Just add the contaminant file to the input file list 'in'." << endl;
    }

    set<String> identifiers; // spot duplicate identifiers  // std::unordered_set<string> has slightly more RAM, but slightly less CPU

    FASTAFile f;
    f.writeStart(out);
    FASTAFile::FASTAEntry entry;

    // Configure Enzymatic digestion
    // TODO: allow user-specified regex
    ProteaseDigestion digestion;
    String enzyme = getStringOption_("enzyme").trim();
    if ((input_type == SeqType::protein) && !enzyme.empty())
    {
      digestion.setEnzyme(enzyme);
    }

    MRMDecoy m;
    m.setParameters(decoy_param);

    for (Size i = 0; i < in.size(); ++i)
    {
      f.readStart(in[i]);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
      while (f.readNext(entry))
      {
        if (identifiers.find(entry.identifier) != identifiers.end())
        {
          OPENMS_LOG_WARN << "DecoyDatabase: Warning, identifier '" << entry.identifier << "' occurs more than once!" << endl;
        }
        identifiers.insert(entry.identifier);

        if (append)
        {
          f.writeNext(entry);
        }

        // identifier
        entry.identifier = getIdentifier_(entry.identifier, decoy_string, decoy_string_position_prefix);

        // sequence
        if (input_type == SeqType::RNA)
        {
          string quick_seq = entry.sequence;
          bool five_p = (entry.sequence.front() == 'p');
          bool three_p = (entry.sequence.back() == 'p');
          if (five_p) //we don't want to reverse terminal phosphates
          {
            quick_seq.erase(0, 1);
          }
          if (three_p)
          {
            quick_seq.pop_back();
          }
          vector<String> tokenized;
          boost::smatch m;
          while (boost::regex_search(quick_seq, m, boost::regex("[^\\[]|(\\[[^\\[\\]]*\\])")))
          {
            tokenized.push_back(m.str(0));
            quick_seq = m.suffix();
          }

          if (shuffle)
          {
            srand(seed);
            random_shuffle(tokenized.begin(), tokenized.end());
          }
          else  // reverse
          {
            reverse(tokenized.begin(), tokenized.end()); //reverse the tokens
          }
          if (five_p)  //add back 5'
          {
            tokenized.insert(tokenized.begin(), String("p"));
          }
          if (three_p) //add back 3'
          {
            tokenized.push_back(String("p"));
          }
          entry.sequence = ListUtils::concatenate(tokenized, "");
        }
        else // protein input
        {
          // if (terminal_aminos != "none")
          if (enzyme != "no cleavage" && (keepN || keepC))
          {
            std::vector<AASequence> peptides;
            digestion.digest(AASequence::fromString(entry.sequence), peptides);
            String new_sequence = "";
            for (auto const& peptide : peptides)
            {
              if (shuffle)
              {
                OpenMS::TargetedExperiment::Peptide p;
                p.sequence = peptide.toString();
                OpenMS::TargetedExperiment::Peptide decoy_p = m.shufflePeptide(p, identity_threshold, seed, max_attempts);
                new_sequence += decoy_p.sequence;
              }
              else
              {
                OpenMS::TargetedExperiment::Peptide p;
                p.sequence = peptide.toString();
                OpenMS::TargetedExperiment::Peptide decoy_p = MRMDecoy::reversePeptide(p, keepN, keepC, keep_const_pattern);
                new_sequence += decoy_p.sequence;
              }
            }
            entry.sequence = new_sequence;
          }
          else
          {
            // sequence
            if (shuffle)
            {
              String temp;
              Size x = entry.sequence.size();
              srand(seed); // identical proteins are shuffled the same way
              while (x != 0)
              {
                Size y = rand() % x;
                temp += entry.sequence[y];
                --x;
                entry.sequence[y] = entry.sequence[x]; // overwrite consumed position with last position (about to go out of scope for next dice roll)
              }
            }
            else // reverse
            {
              entry.sequence.reverse();
            }
          }
        }

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        f.writeNext(entry);
      } // next protein
    } // input files

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPDecoyDatabase tool;
  return tool.main(argc, argv);
}

/// @endcond
