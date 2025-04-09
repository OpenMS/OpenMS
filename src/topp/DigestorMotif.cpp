// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
// $Author: Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <map>

using namespace OpenMS;
using namespace std;

#define SEP "\t"

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DigestorMotif DigestorMotif

@brief This application is used to digest a protein database to get all peptides given a cleavage enzyme. It will also produce peptide statistics given the mass
accuracy of the instrument. You can extract peptides with specific motifs,e.g. onyl cysteine containing peptides for ICAT experiments. At the moment only trypsin is supported.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DigestorMotif.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DigestorMotif.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDigestorMotif :
  public TOPPBase
{
public:
  TOPPDigestorMotif() :
    TOPPBase("DigestorMotif", "Digests a protein database in-silico")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "FASTA input file");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "output file (peptides)\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerIntOption_("missed_cleavages", "<number>", 1, "the number of allowed missed cleavages", false);
    registerIntOption_("mass_accuracy", "<number>", 1000, "give your mass accuracy in ppb", false);
    registerIntOption_("min_length", "<number>", 6, "minimum length of peptide", false);
    registerIntOption_("out_option", "<number>", 1, "indicate 1 (peptide table only), 2 (statistics only) or (both peptide table + statistics)", false);
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("enzyme", all_enzymes);
    registerStringOption_("motif", "<string>", "M", "the motif for the restricted peptidome", false);
    setMinInt_("missed_cleavages", 0);
  }

  ExitCodes main_(int, const char**) override
  {
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    std::vector<FASTAFile::FASTAEntry> protein_data;
    FASTAFile file;
    ProteaseDigestion digestor;
    vector<AASequence> temp_peptides;
    PeptideIdentification peptide_identification;
    ProteinIdentification protein_identification;
    PeptideHit temp_peptide_hit;
    ProteinHit temp_protein_hit;
    vector<String> protein_accessions;
    String inputfile_name;
    String outputfile_name;
    UInt min_size, counter = 0;
    UInt missed_cleavages;
    double accurate_mass, min_mass, max_mass;
    UInt mass_acc, out_opt;
    EmpiricalFormula EF;
    UInt zero_count = 0;
    ProteinIdentification::SearchParameters search_parameters;

    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    inputfile_name = getStringOption_("in");
    outputfile_name = getStringOption_("out");
    min_size = getIntOption_("min_length");
    mass_acc = getIntOption_("mass_accuracy");
    out_opt = getIntOption_("out_option");
    missed_cleavages = getIntOption_("missed_cleavages");
    AASequence M = AASequence::fromString(getStringOption_("motif"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    file.load(inputfile_name, protein_data);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // This should be updated if more cleavage enzymes are available
    String enzyme_name = getStringOption_("enzyme");
    digestor.setEnzyme(enzyme_name);
    search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_name));
    digestor.setMissedCleavages(missed_cleavages);

    for (UInt i = 0; i < protein_data.size(); ++i)
    {
      PeptideEvidence pe;
      temp_protein_hit.setSequence(protein_data[i].sequence);
      temp_protein_hit.setAccession(protein_data[i].identifier);
      pe.setProteinAccession(protein_data[i].identifier);
      digestor.digest(AASequence::fromString(protein_data[i].sequence), temp_peptides);
      temp_peptide_hit.setPeptideEvidences(vector<PeptideEvidence>(1, pe));
      for (UInt j = 0; j < temp_peptides.size(); ++j)
      {
        if (temp_peptides[j].size() >= min_size)
        {
          if (temp_peptides[j].hasSubsequence(M))
          {
            temp_peptide_hit.setSequence(temp_peptides[j]);
            peptide_identification.insertHit(temp_peptide_hit);
          }
        }
      }
      protein_identifications[0].insertHit(temp_protein_hit);
    }
    DateTime date_time;
    String date_time_string;
    date_time.now();

    date_time_string = date_time.get();
    protein_identifications[0].setSearchParameters(search_parameters);
    protein_identifications[0].setDateTime(date_time);
    protein_identifications[0].setSearchEngine("In-silico digestion");
    protein_identifications[0].setIdentifier("In-silico_digestion" + date_time_string);
    peptide_identification.setIdentifier("In-silico_digestion" + date_time_string);
    identifications.push_back(peptide_identification);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    ofstream fp_out(outputfile_name.c_str());
    if (out_opt == 2)
    {
      fp_out << "mass_error" << SEP << "#proteins in database" << SEP << "# tryptic peptides" << SEP << "# unique peptide weights" << SEP << "# identifiable proteins" << SEP << "average window_size" << "\n";
    }
    UInt mass_iter = mass_acc;
    while (mass_iter > 0)
    {
      vector<double> MIN, MAX;
      vector<String> protein_names, PROTEINS;
      vector<vector<double> > Y;
      vector<UInt> OVER;
      UInt total = 0;
      if (out_opt == 1 || out_opt == 3)
      {
        fp_out << "counter" << SEP << "ProteinID" << SEP << "PeptideLocation" << SEP << "PeptideSequence" << SEP << "C" << SEP << "H" << SEP << "N" << SEP << "O" << SEP << "S" << SEP << "length" << SEP << "weight" << SEP << "min_weight" << SEP << "max_weight" << SEP << "Formula" << SEP << "D" << SEP << "E" << SEP << "K" << SEP << "R" << SEP << "H" << SEP << "Y" << SEP << "W" << SEP << "F" << SEP << "C" << SEP << "M" << SEP << "S" << SEP << "T" << SEP << "N" << SEP << "Q" << SEP << "G" << SEP << "A" << SEP << "V" << SEP << "L" << SEP << "I" << SEP << "P" << SEP << "hydrophobicity" << "\n";
      }

      for (UInt i = 0; i < protein_data.size(); ++i)
      {
        PeptideEvidence pe;
        pe.setProteinAccession(protein_data[i].identifier);
        temp_protein_hit.setAccession(protein_data[i].identifier);
        digestor.digest(AASequence::fromString(protein_data[i].sequence), temp_peptides);
        temp_peptide_hit.setPeptideEvidences(vector<PeptideEvidence>(1, pe));
        for (UInt j = 0; j < temp_peptides.size(); ++j)
        {
          //vector<double> B_peptide, Y_peptide;
          vector<double> peptide_ions;
          accurate_mass = temp_peptides[j].getMonoWeight();
          min_mass = accurate_mass - mass_iter * accurate_mass / 1000000000;
          max_mass = accurate_mass + mass_iter * accurate_mass / 1000000000;
          EF = temp_peptides[j].getFormula();
          for (UInt r = 1; r <= temp_peptides[j].size(); ++r)
          {
            //B_peptide.push_back(temp_peptides[j].getPrefix(r).getMonoWeight());
            peptide_ions.push_back(temp_peptides[j].getPrefix(r).getMonoWeight());
            peptide_ions.push_back(temp_peptides[j].getSuffix(r).getMonoWeight());
            //Y_peptide.push_back(temp_peptides[j].getSuffix(r).getMonoWeight());
          }
          if (temp_peptides[j].size() >= min_size)
          {
            if (temp_peptides[j].hasSubsequence(M))
            {
              OVER.push_back((-1)); //because the increment of the first will always be counted;
              //IonCounter.push_back(0);
              MIN.push_back(min_mass);
              MAX.push_back(max_mass);
              Y.push_back(peptide_ions);
              //B.push_back(B_peptide);
              protein_names.push_back(protein_accessions[0]);
              temp_peptide_hit.setSequence(temp_peptides[j]);
              peptide_identification.insertHit(temp_peptide_hit);
              if (out_opt == 1 || out_opt == 3)
              {
                const String unmodified_peptide = temp_peptides[j].toUnmodifiedString();
                const Size nK = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'K');
                const Size nD = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'D');
                const Size nR = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'R');
                const Size nW = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'W');
                const Size nM = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'M');
                const Size nN = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'N');
                const Size nA = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'A');
                const Size nI = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'I');
                const Size nE = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'E');
                const Size nH = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'H');
                const Size nF = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'F');
                const Size nS = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'S');
                const Size nQ = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'Q');
                const Size nV = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'V');
                const Size nP = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'P');
                const Size nY = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'Y');
                const Size nC = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'C');
                const Size nT = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'T');
                const Size nG = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'G');
                const Size nL = std::count(unmodified_peptide.begin(), unmodified_peptide.end(), 'L');

                const ElementDB* db = ElementDB::getInstance();
                const Element* C = db->getElement("C");
                const Element* H = db->getElement("H");
                const Element* N = db->getElement("N");
                const Element* O = db->getElement("O");
                const Element* S = db->getElement("S");

                fp_out << counter << SEP << ">" << protein_accessions[0] << SEP << j << SEP << temp_peptides[j] << SEP
                       << EF.getNumberOf(C) << SEP << EF.getNumberOf(H) << SEP << EF.getNumberOf(N) << SEP << EF.getNumberOf(O) << SEP
                       << EF.getNumberOf(S) << SEP << temp_peptides[j].size() << SEP << precisionWrapper(temp_peptides[j].getMonoWeight()) << SEP
                       << precisionWrapper(min_mass) << SEP << precisionWrapper(max_mass) << SEP << temp_peptides[j].getFormula() << SEP
                       << nD << SEP << nE << SEP << nK << SEP << nR << SEP << nH << SEP << nY << SEP << nW << SEP << nF << SEP << nC << SEP
                       << nM << SEP << nS << SEP << nT << SEP << nN << SEP << nQ << SEP << nG << SEP << nA << SEP << nV << SEP << nL << SEP
                       << nI << SEP << nP << SEP << nD * (-3.5) + nE * (-3.5) + nK * (-3.9) + nR * (-4.5) + nH * (-3.2) + nY * (-1.3) + nW * (-0.9) + nF * (2.8) + nC * (2.5) + nM * (1.9) + nS * (-0.8) + nT * (-0.7) + nN * (-3.5) + nQ * (-3.5) + nG * (-0.4) + nA * (1.8) + nV * (4.2) + nL * (4.5) + nI * (4.5) + nP * (-1.6) << "\n";
              }
              counter++;
            }
          }
        }
        protein_identifications[0].insertHit(temp_protein_hit);
      }
      if (out_opt != 2)
      {
        fp_out << "MW_count" << SEP;
      }

      for (UInt r = 1; r <= 100; ++r)
      {
        fp_out << "y" << r << SEP << "b" << r << SEP;
      }
      fp_out << "\n";
      fp_out << "MW_count" << SEP << "Overlapping ions in search space" << "\n";
      for (UInt x = 0; x < MAX.size(); ++x)
      {
        /*for(UInt it = 0; it < ions.size(); ++it)
        {
            ions[it] = -1; //all ions from the same peptide all be counted
        }
        */
        cout << "2nd loop" << SEP << MAX.size() - x << endl;
        vector<UInt> IonCounter;
        for (UInt y = 0; y < MAX.size(); ++y)
        {
          if ((MIN[y] < MIN[x] && MAX[y] > MIN[x]) || (MAX[y] > MAX[x] && MIN[y] < MAX[x]) || (MIN[x] == MIN[y]))
          {
            OVER[x] = OVER[x] + 1;
            //find overlapping tandem ions
            vector<double> X_temp, Y_temp;
            X_temp = Y[x];
            Y_temp = Y[y];
            UInt ions = 0;
            for (UInt xx = 0; xx < X_temp.size(); ++xx)
            {
              for (UInt yy = 0; yy < Y_temp.size(); ++yy)
              {
                if (fabs(X_temp[xx] - Y_temp[yy]) <= 1)
                {
                  ions += 1;
                }
              }
            }
            IonCounter.push_back(ions);
          }
        }
        if (out_opt == 3)
        {
          fp_out << OVER[x] << SEP;
          if (MAX[x] < 3500)
          {
            for (UInt it = 0; it < IonCounter.size(); ++it)
            {
              fp_out << IonCounter[it] << SEP;
            }
          }
          fp_out << "\n";
          cout << OVER[x];
        }
        total = total + OVER[x];
        if (OVER[x] == 0)
        {
          ++zero_count;
          PROTEINS.push_back(protein_names[x]);
        }
      }
      UInt pro_count = 0;
      for (UInt a = 0; a + 1 < PROTEINS.size(); ++a)
      {
        if (PROTEINS[a] == PROTEINS[a + 1])
        {
          ++pro_count;
        }
        cout << PROTEINS.size() << endl << pro_count << endl;
      }

      if (out_opt != 2)
      {
        mass_iter = 0;
      }
      else
      {
        mass_iter = mass_iter - 1;
      }

      if (out_opt > 1)
      {
        fp_out << mass_iter << SEP << protein_data.size() << SEP << MAX.size() << SEP << zero_count << SEP << PROTEINS.size() - pro_count << SEP << total << endl;
      }
      pro_count = 0;
      zero_count = 0;
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPDigestorMotif tool;
  return tool.main(argc, argv);
}

/// @endcond
