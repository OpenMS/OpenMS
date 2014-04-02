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
// $Author: Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <map>

using namespace OpenMS;
using namespace std;

#define SEP "\t"

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_DigestorMotif DigestorMotif

    @brief This application is used to digest a protein database to get all peptides given a cleavage enzyme. It will also produce peptide statistics given the mass
    accuracy of the instrument. You can extract peptides with specific motifs,e.g. onyl cysteine containing peptides for ICAT experiments. At the moment only trypsin is supported.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_DigestorMotif.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_DigestorMotif.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDigestorMotif :
  public TOPPBase
{
public:
  TOPPDigestorMotif() :
    TOPPBase("DigestorMotif", "digests a protein database in-silico", false)
  {

  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "FASTA input file");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "output file (peptides)\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerIntOption_("missed_cleavages", "<number>", 1, "the number of allowed missed cleavages", false);
    registerIntOption_("mass_accuracy", "<number>", 1000, "give your mass accuracy in ppb", false);
    registerIntOption_("min_length", "<number>", 6, "minimum length of peptide", false);
    registerIntOption_("out_option", "<number>", 1, "indicate 1 (peptide table only), 2 (statistics only) or (both peptide table + statistics)", false);
    registerStringOption_("enzyme", "<string>", "Trypsin", "the digestion enzyme", false);
    registerStringOption_("motif", "<string>", "M", "the motif for the restricted peptidome", false);
    setMinInt_("missed_cleavages", 0);
  }

  ExitCodes main_(int, const char**)
  {
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    std::vector<FASTAFile::FASTAEntry> protein_data;
    FASTAFile file;
    EnzymaticDigestion digestor;
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
    AASequence M = AASequence(getStringOption_("motif"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    file.load(inputfile_name, protein_data);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // This should be updated if more cleavage enzymes are available
    digestor.setEnzyme(EnzymaticDigestion::ENZYME_TRYPSIN);
    search_parameters.enzyme = ProteinIdentification::TRYPSIN;
    digestor.setMissedCleavages(missed_cleavages);

    protein_accessions.resize(1, String(""));
    for (UInt i = 0; i < protein_data.size(); ++i)
    {
      protein_accessions[0] = protein_data[i].identifier;
      temp_protein_hit.setSequence(protein_data[i].sequence);
      temp_protein_hit.setAccession(protein_accessions[0]);

      digestor.digest(AASequence(protein_data[i].sequence), temp_peptides);
      temp_peptide_hit.setProteinAccessions(protein_accessions);
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
        protein_accessions[0] = protein_data[i].identifier;
        temp_protein_hit.setAccession(protein_accessions[0]);
        digestor.digest(AASequence(protein_data[i].sequence), temp_peptides);
        temp_peptide_hit.setProteinAccessions(protein_accessions);
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
                fp_out << counter << SEP << ">" << protein_accessions[0] << SEP << j << SEP << temp_peptides[j] << SEP << EF.getNumberOf("C") << SEP << EF.getNumberOf("H") << SEP << EF.getNumberOf("N") << SEP << EF.getNumberOf("O") << SEP << EF.getNumberOf("S") << SEP << temp_peptides[j].size() << SEP << precisionWrapper(temp_peptides[j].getMonoWeight()) << SEP << precisionWrapper(min_mass) << SEP << precisionWrapper(max_mass) << SEP << temp_peptides[j].getFormula() << SEP << temp_peptides[j].getNumberOf("D") << SEP << temp_peptides[j].getNumberOf("E") << SEP << temp_peptides[j].getNumberOf("K") << SEP << temp_peptides[j].getNumberOf("R") << SEP << temp_peptides[j].getNumberOf("H") << SEP << temp_peptides[j].getNumberOf("Y") << SEP << temp_peptides[j].getNumberOf("W") << SEP << temp_peptides[j].getNumberOf("F") << SEP << temp_peptides[j].getNumberOf("C") << SEP << temp_peptides[j].getNumberOf("M") << SEP << temp_peptides[j].getNumberOf("S") << SEP << temp_peptides[j].getNumberOf("T") << SEP << temp_peptides[j].getNumberOf("N") << SEP << temp_peptides[j].getNumberOf("Q") << SEP << temp_peptides[j].getNumberOf("G") << SEP << temp_peptides[j].getNumberOf("A") << SEP << temp_peptides[j].getNumberOf("V") << SEP << temp_peptides[j].getNumberOf("L") << SEP << temp_peptides[j].getNumberOf("I") << SEP << temp_peptides[j].getNumberOf("P") << SEP << temp_peptides[j].getNumberOf("D") * (-3.5) + temp_peptides[j].getNumberOf("E") * (-3.5) + temp_peptides[j].getNumberOf("K") * (-3.9) + temp_peptides[j].getNumberOf("R") * (-4.5) + temp_peptides[j].getNumberOf("H") * (-3.2) + temp_peptides[j].getNumberOf("Y") * (-1.3) + temp_peptides[j].getNumberOf("W") * (-0.9) + temp_peptides[j].getNumberOf("F") * (2.8) + temp_peptides[j].getNumberOf("C") * (2.5) + temp_peptides[j].getNumberOf("M") * (1.9) + temp_peptides[j].getNumberOf("S") * (-0.8) + temp_peptides[j].getNumberOf("T") * (-0.7) + temp_peptides[j].getNumberOf("N") * (-3.5) + temp_peptides[j].getNumberOf("Q") * (-3.5) + temp_peptides[j].getNumberOf("G") * (-0.4) + temp_peptides[j].getNumberOf("A") * (1.8) + temp_peptides[j].getNumberOf("V") * (4.2) + temp_peptides[j].getNumberOf("L") * (4.5) + temp_peptides[j].getNumberOf("I") * (4.5) + temp_peptides[j].getNumberOf("P") * (-1.6) << "\n";
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
      for (UInt a = 0; a < PROTEINS.size() - 1; ++a)
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
