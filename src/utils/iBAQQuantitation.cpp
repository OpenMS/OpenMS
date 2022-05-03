 
// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Maintainer $
// $Authors: Author1, Author2 $
// --------------------------------------------------------------------------
 
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>

using namespace OpenMS;
using namespace std;
 
//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
 
// We do not want this class to show up in the docu:
 
class TOPPiBAQQuantitation :
    public TOPPBase
{
public:
  TOPPiBAQQuantitation() :
      TOPPBase("iBAQQuantitation", "Intensity based absolute quantification (iBAQ) of proteins.", false)
  {
 
  }
 
protected:
 
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input FASTA file, containing a database.");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerInputFile_("id", "<file>", "", "Input file containing identified peptides and proteins.");
    setValidFormats_("id", ListUtils::create<String>("mzTab"));
    registerOutputFile_("out", "<file>", "", "Output mzTab including iBAQ quantities.");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));
  }
 
 
  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String id(getStringOption_("id"));
    String out(getStringOption_("out"));
 
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
//    FileTypes::Type filetype = FileHandler::getType(in);
//    if (filetype == FileTypes::MZTAB)

    vector<FASTAFile::FASTAEntry> db;
    FASTAFile().load(in, db);

    MzTab mztab;
    MzTabFile mztab_file;
    mztab_file.load(id, mztab);
 
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // Enzymatic digestion
    map<StringView, double> protein_digests;

    for (FASTAFile::FASTAEntry const& db_entry : db)
    {
      vector<StringView> digested;

      StringView identifier{ db_entry.identifier };
      StringView sequence{ db_entry.sequence };

      // TODO: Should single amino acids be excluded?
      EnzymaticDigestion().digestUnmodified(sequence, digested);
      protein_digests[identifier] = static_cast<double>(digested.size());
    }

    // Sum of intensities
    // TODO: StringView would probably be better, but doesn't work somehow. Also, replace double with long double (potentially large numbers)?
    map<string, map<const unsigned long, double>> protein_abundances;

    for (MzTabPeptideSectionRow &peptide_row : mztab.getPeptideSectionRows())
    {
      const string &accession{ peptide_row.accession.get() };

      if (protein_abundances.find(accession) == protein_abundances.end())
      {
        protein_abundances[accession] = map<const unsigned long, double>{};
      }

      for (pair<const unsigned long, MzTabDouble> &peptide_abundance_entry : peptide_row.peptide_abundance_study_variable)
      {
        const unsigned long &study_variable{ peptide_abundance_entry.first };
        double peptide_abundance{ (!peptide_abundance_entry.second.isNull()) ? peptide_abundance_entry.second.get() : 0. };

        protein_abundances[accession][study_variable] += peptide_abundance;
      }
    }

    // Divide by theoretical number of peptides
    for (pair<const string, map<const unsigned long, double>> &accession_abundances : protein_abundances)
    {
      for (pair<const unsigned long, double> &study_abundance : accession_abundances.second)
      {
        study_abundance.second /= protein_digests[accession_abundances.first];
      }
    }


//    for (auto &pair1 : protein_abundances)
//    {
//      cout << pair1.first << endl;
//      for (auto &pair2 : pair1.second)
//      {
//        cout << pair2.first << " " << pair2.second << endl;
//      }
//    }

 
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    mztab_file.store(out, mztab);


    return EXECUTION_OK;
  }
 
};
 
 
// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPiBAQQuantitation tool;
  return tool.main(argc, argv);
}