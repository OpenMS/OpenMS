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
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_SequenceCoverageCalculator SequenceCoverageCalculator

    @brief Prints information about idXML files.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_SequenceCoverageCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SequenceCoverageCalculator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceCoverageCalculator :
  public TOPPBase
{
public:
  TOPPSequenceCoverageCalculator() :
    TOPPBase("SequenceCoverageCalculator", "Prints information about idXML files.", false)
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_database", "<file>", "", "input file containing the database in FASTA format");
    setValidFormats_("in_database", ListUtils::create<String>("fasta"));
    registerInputFile_("in_peptides", "<file>", "", "input file containing the identified peptides");
    setValidFormats_("in_peptides", ListUtils::create<String>("idXML"), true);
    registerOutputFile_("out", "<file>", "", "Optional text output file. If left out, the output is written to the command line.", false);
    setValidFormats_("out", ListUtils::create<String>("txt"));
  }

  void getStartAndEndIndex(const String& sequence, const String& substring, pair<Size, Size>& indices)
  {
    indices.first = 0;
    indices.second = 0;

    if (sequence.hasSubstring(substring))
    {
      for (Size i = 0; i <= sequence.size() - substring.size(); ++i)
      {
        Size temp_index = i;
        Size temp_count = 0;
        while (temp_index < sequence.size()
              && temp_count < substring.size()
              && sequence.at(temp_index) == substring.at(temp_index - i))
        {
          ++temp_index;
          ++temp_count;
        }
        if (temp_count == substring.size())
        {
          indices.first = i;
          indices.second = temp_index;
          i = sequence.size();
        }
      }
    }
  }

  ExitCodes outputTo_(ostream& os)
  {
    IdXMLFile idXML_file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<FASTAFile::FASTAEntry> proteins;
    vector<double> statistics;
    vector<Size> counts;
    vector<Size> mod_counts;
    vector<PeptideHit> temp_hits;
    vector<Size> coverage;
    Size spectrum_count = 0;
    map<String, Size> unique_peptides;
    map<String, Size> temp_unique_peptides;
    map<String, Size> temp_modified_unique_peptides;

    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in_peptides");
    String database_name = getStringOption_("in_database");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    String document_id;
    idXML_file.load(inputfile_name, protein_identifications, identifications, document_id);
    FASTAFile().load(database_name, proteins);

    statistics.resize(proteins.size(), 0.);
    counts.resize(proteins.size(), 0);
    mod_counts.resize(proteins.size(), 0);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    for (Size j = 0; j < proteins.size(); ++j)
    {
      coverage.clear();
      coverage.resize(proteins[j].sequence.size(), 0);
      temp_unique_peptides.clear();
      temp_modified_unique_peptides.clear();

      for (Size i = 0; i < identifications.size(); ++i)
      {
        if (!identifications[i].empty())
        {
          if (identifications[i].getHits().size() > 1)
          {
            LOG_ERROR << "Spectrum with more than one identification found, which is not allowed.\n"
                      << "Use the IDFilter with the -best_hits option to filter for best hits." << endl;
            return ILLEGAL_PARAMETERS;
          }

          set<String> accession;
          accession.insert(proteins[j].identifier);
          temp_hits = PeptideIdentification::getReferencingHits(identifications[i].getHits(), accession);

          if (temp_hits.size() == 1)
          {
            pair<Size, Size> indices;
            getStartAndEndIndex(proteins[j].sequence, temp_hits[0].getSequence().toUnmodifiedString(), indices);
            for (Size k = indices.first; k < indices.second; ++k)
            {
              coverage[k] = 1;
            }
            if (indices.first != indices.second)
            {
              // os <<  temp_hits[0].getSequence().toUnmodifiedString() << endl;
            }
            ++spectrum_count;
            if (unique_peptides.find(temp_hits[0].getSequence().toString()) == unique_peptides.end())
            {
              unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
            }
            if (temp_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_unique_peptides.end())
            {
              temp_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toUnmodifiedString(), 0));
            }
            if (temp_modified_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_modified_unique_peptides.end())
            {
              temp_modified_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
            }
          }
        }
      }
/* << proteins[j].sequence << endl;
                for (Size k = 0; k < coverage.size(); ++k)
                {
                    os << coverage[k];
                }
                os << endl;
*/
      // statistics[j] = make_pair(,
      // accumulate(coverage.begin(), coverage.end(), 0) / proteins[j].sequence.size());
      statistics[j] = ((double) accumulate(coverage.begin(), coverage.end(), Size(0))) / proteins[j].sequence.size();
      counts[j] = temp_unique_peptides.size();
      mod_counts[j] = temp_modified_unique_peptides.size();

      // details for this protein
      if (counts[j] > 0)
      {
        os << proteins[j].identifier << "(coverage%, #unique hits): " <<  statistics[j] * 100 << "%, " << counts[j] << "\n";
      }

// os << statistics[j] << endl;
    }

// os << "Sum of coverage is " << accumulate(statistics.begin(), statistics.end(), 0.) << endl;
    os << "Average coverage per protein is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
    os << "Average number of peptides per protein is " << (((double) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
    os << "Average number of un/modified peptides per protein is " << (((double) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;
    os << "Number of identified spectra: " << spectrum_count << endl;
    os << "Number of unique identified peptides: " << unique_peptides.size() << endl;

    vector<double>::iterator it = statistics.begin();
    vector<Size>::iterator it2 = counts.begin();
    vector<Size>::iterator it3 = mod_counts.begin();
    while (it != statistics.end())
    {
      if (*it == 0.)
      {
        it = statistics.erase(it);
        it2 = counts.erase(it2);
        it3 = mod_counts.erase(it3);
      }
      else
      {
        ++it;
        ++it2;
        ++it3;
      }
    }
    os << "Average coverage per found protein (" << statistics.size() << ") is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
    os << "Average number of peptides per found protein is " << (((double) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
    os << "Average number of un/modified peptides per protein is " << (((double) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;

    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    String out = getStringOption_("out");

    TOPPBase::ExitCodes ret;
    if (out != "")
    {
      ofstream os(out.c_str());
      ret = outputTo_(os);
      os.close();
    }
    else
    {
      ret = outputTo_(LOG_INFO);
    }

    return ret;
  }

};


int main(int argc, const char** argv)
{
  TOPPSequenceCoverageCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
