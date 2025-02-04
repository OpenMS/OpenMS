// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <iostream>
#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDExtractor IDExtractor

@brief Extracts 'n' peptides randomly or best 'n' from idXML files.

Input and output format are 'idXML'. The tools allows you to extract subsets of peptides
from idXML files.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDExtractor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDExtractor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDExtractor :
  public TOPPBase
{
public:
  TOPPIDExtractor() :
    TOPPBase("IDExtractor", "Extracts 'n' peptides randomly or best 'n' from idXML files.")
  {

  }

  static bool compareIDsWithScores(pair<double, PeptideIdentification>& a, pair<double, PeptideIdentification>& b)
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
    registerIntOption_("number_of_rand_invokations", "<int>", 0, "Number of rand invocations before random draw (basically a seed)", false);
    setMinInt_("number_of_rand_invokations", 0);
    registerFlag_("best_hits", "If this flag is set the best n peptides are chosen.");
  }

  ExitCodes main_(int, const char**) override
  {
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
    FileHandler().loadIdentifications(inputfile_name, protein_identifications, identifications, {FileTypes::IDXML});

    if (number_of_peptides > identifications.size())
    {
      writeLogError_("Number of existing peptides smaller than number of chosen peptides. Aborting!");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    if (best_hits)
    {
      for (Size i = 0; i < identifications.size(); ++i)
      {
        identifications_with_scores.emplace_back(identifications[i].getHits()[0].getScore(), identifications[i]);
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
      Math::RandomShuffler r(number_of_rand_invokations);
      r.portable_random_shuffle(indices.begin(), indices.end());

      Size index = 0;
      while (chosen_ids.size() < number_of_peptides && index < indices.size())
      {
        if (!identifications[indices[index]].getHits().empty() && find(chosen_ids.begin(), chosen_ids.end(), identifications[indices[index]].getHits()[0].getSequence().toString()) == chosen_ids.end())
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
      writeLogError_("Number of existing unique peptides (" + String(chosen_ids.size()) + ") smaller than number of chosen peptides. Aborting!");
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

    FileHandler().storeIdentifications(outputfile_name,
                     chosen_protein_identifications,
                     chosen_identifications,
                     {FileTypes::IDXML});

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
