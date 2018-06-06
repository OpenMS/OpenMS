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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <map>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IDExtractor IDExtractor

    @brief Extracts 'n' peptides randomly or best 'n' from idXML files.

  Input and output format are 'idXML'. The tools allows you to extract subsets of peptides
  from idXML files.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDExtractor.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDExtractor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDExtractor :
  public TOPPBase
{
public:
  TOPPIDExtractor() :
    TOPPBase("IDExtractor", "Extracts 'n' peptides randomly or best 'n' from idXML files.", false)
  {

  }

  static bool compareIDsWithScores(pair<double, PeptideIdentification> a, pair<double, PeptideIdentification> b)
  {
    if (a.second.isHigherScoreBetter())
    {
      return a.first > b.first;
    }
    else
    {
      return a.first < b.first;
    }
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerIntOption_("number_of_peptides", "<int>", 10, "Number of randomly chosen peptides", false);
    setMinInt_("number_of_peptides", 1);
    registerIntOption_("number_of_rand_invokations", "<int>", 0, "Number of rand invocations before random draw", false);
    setMinInt_("number_of_rand_invokations", 0);
    registerFlag_("best_hits", "If this flag is set the best n peptides are chosen.");
  }

  ExitCodes main_(int, const char**) override
  {
    IdXMLFile idXML_file;
    vector<ProteinIdentification> protein_identifications;
    vector<ProteinIdentification> chosen_protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<PeptideIdentification> chosen_identifications;
    vector<Size> indices;
    vector<PeptideHit> temp_peptide_hits;
    vector<ProteinHit> temp_protein_hits;
    vector<ProteinHit> chosen_protein_hits;
    map<String, vector<PeptideIdentification> > identifiers;
    PeptideIdentification temp_identification;
    vector<String> chosen_ids;
    vector<pair<double, PeptideIdentification> > identifications_with_scores;
    vector<pair<double, PeptideIdentification> >::iterator it = identifications_with_scores.begin();
    vector<PeptideIdentification> temp_identifications;


    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");
    Size number_of_peptides = getIntOption_("number_of_peptides");
    Size number_of_rand_invokations = getIntOption_("number_of_rand_invokations");
    bool best_hits = getFlag_("best_hits");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    String document_id;
    idXML_file.load(inputfile_name, protein_identifications, identifications, document_id);

    if (number_of_peptides > identifications.size())
    {
      writeLog_("Number of existing peptides smaller than number of chosen peptides. Aborting!");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    if (best_hits)
    {
      for (Size i = 0; i < identifications.size(); ++i)
      {
        identifications_with_scores.push_back(make_pair(identifications[i].getHits()[0].getScore(), identifications[i]));
      }
      sort(identifications_with_scores.begin(), identifications_with_scores.end(), TOPPIDExtractor::compareIDsWithScores);
      it = identifications_with_scores.begin();
      while (it != identifications_with_scores.end() && chosen_ids.size() < number_of_peptides)
      {
        if (find(chosen_ids.begin(), chosen_ids.end(), it->second.getHits()[0].getSequence().toString()) == chosen_ids.end())
        {
          chosen_ids.push_back(it->second.getHits()[0].getSequence().toString());
          chosen_identifications.push_back(it->second);
          if (identifiers.find(it->second.getIdentifier()) == identifiers.end())
          {
            temp_identifications.clear();
          }
          else
          {
            temp_identifications = identifiers[it->second.getIdentifier()];
            identifiers.erase(it->second.getIdentifier());
          }
          temp_identifications.push_back(it->second);
          identifiers.insert(make_pair(it->second.getIdentifier(), temp_identifications));
        }
        ++it;
      }
    }
    else
    {
      indices.resize(identifications.size(), 0);
      for (Size i = 0; i < identifications.size(); ++i)
      {
        indices[i] = i;
      }
      for (Size i = 0; i < number_of_rand_invokations; ++i)
      {
        rand();
      }
      random_shuffle(indices.begin(), indices.end());

      Size index = 0;
      while (chosen_ids.size() < number_of_peptides && index < indices.size())
      {
        if (identifications[indices[index]].getHits().size() > 0 && find(chosen_ids.begin(), chosen_ids.end(), identifications[indices[index]].getHits()[0].getSequence().toString()) == chosen_ids.end())
        {
          chosen_ids.push_back(identifications[indices[index]].getHits()[0].getSequence().toString());
          chosen_identifications.push_back(identifications[indices[index]]);
          if (identifiers.find(identifications[indices[index]].getIdentifier()) == identifiers.end())
          {
            temp_identifications.clear();
          }
          else
          {
            temp_identifications = identifiers[identifications[indices[index]].getIdentifier()];
            identifiers.erase(identifications[indices[index]].getIdentifier());
          }
          temp_identifications.push_back(identifications[indices[index]]);
          identifiers.insert(make_pair(identifications[indices[index]].getIdentifier(), temp_identifications));
        }
        ++index;
      }
    }

    if (chosen_ids.size() < number_of_peptides)
    {
      writeLog_("Number of existing unique peptides (" + String(chosen_ids.size()) + ") smaller than number of chosen peptides. Aborting!");
      return ILLEGAL_PARAMETERS;
    }

    for (Size i = 0; i < protein_identifications.size(); ++i)
    {
      temp_protein_hits = protein_identifications[i].getHits();
      chosen_protein_hits.clear();
      if (identifiers.find(protein_identifications[i].getIdentifier()) != identifiers.end())
      {
        temp_identifications = identifiers[protein_identifications[i].getIdentifier()];
        for (Size j = 0; j < temp_protein_hits.size(); ++j)
        {
          bool already_chosen = false;
          for (Size k = 0; k < temp_identifications.size(); ++k)
          {
            temp_peptide_hits.clear();
            set<String> accession;
            accession.insert(temp_protein_hits[j].getAccession());
            temp_peptide_hits = PeptideIdentification::getReferencingHits(temp_identifications[k].getHits(), accession);
            if (!temp_peptide_hits.empty() && !already_chosen)
            {
              chosen_protein_hits.push_back(temp_protein_hits[j]);
              already_chosen = true;
            }
          }
        }
        if (chosen_protein_hits.empty())
        {
          cout << "No protein hits found for " << protein_identifications[i].getIdentifier()
               << " although having " << temp_identifications.size() << " ids" << endl;
        }
        protein_identifications[i].setHits(chosen_protein_hits);
        chosen_protein_identifications.push_back(protein_identifications[i]);
      }
    }

    idXML_file.store(outputfile_name,
                     chosen_protein_identifications,
                     chosen_identifications);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
