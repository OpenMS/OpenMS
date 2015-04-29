// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_NucleotideID NucleotideID

    @brief Generate a csv file containing masses    sequences   Identifiers from a fasta file

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NucleotideID.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NucleotideID.html
*/

class TOPPNucleotideID :
  public TOPPBase
{

public:
  TOPPNucleotideID() :
    TOPPBase("NucleotideID", "Tool for nucleotide chain identification.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_fasta", "<file>", "", "Fasta file containing the nucleotide library\n");
    setValidFormats_("in_fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_db_mapping", "<file>", "", "mapping file used in AccurateMassSearch\n");
    setValidFormats_("out_db_mapping", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_struct_mapping", "<file>", "", "structure file used in AccurateMassSearch\n");
    setValidFormats_("out_struct_mapping", ListUtils::create<String>("tsv"));

// TODO: other AccurateMassSearch files
  }


  ExitCodes main_(int, const char**)
  {
    // TODO: CODE
    //
    //
    String input_filepath(getStringOption_("in_fasta"));

    std::vector<FASTAFile::FASTAEntry> input_file;
    FASTAFile fasta_file;
    fasta_file.load(input_filepath,input_file);

    TextFile db_mapping_file;
    TextFile struct_mapping_file;

    //generate monophosphate mass list
    map<char, EmpiricalFormula> monophosphate_to_formula;
    monophosphate_to_formula['A']=EmpiricalFormula("C10H14N5O7P");
    monophosphate_to_formula['C']=EmpiricalFormula("C9H14N3O8P");
    monophosphate_to_formula['G']=EmpiricalFormula("C10H14N5O8P");
    monophosphate_to_formula['U']=EmpiricalFormula("C9H13N2O9P");

    db_mapping_file.addLine("database_name\tmirna\ndatabase_version\t1.0"); //header for db_mapping file

    for (size_t i=0;i<input_file.size();i++) //for each sequence
    {
        EmpiricalFormula entryformula;
        string seq= input_file.at(i).sequence;
        int seqlen= seq.length();
        for (int j=0;j<seqlen-1;j++){ //for each residue
            entryformula=entryformula+monophosphate_to_formula[seq[j]]-EmpiricalFormula("H2O"); //Add the mass of the next nucleotide - the condensed water
        }
        entryformula=entryformula+monophosphate_to_formula[seq[seqlen-1]]-EmpiricalFormula("HPO3"); //Add the mass of the next nucleotide - the condensed water, drop the phosphate


        db_mapping_file.addLine(String(entryformula.getMonoWeight()) + "\t" + entryformula.toString() + "\t" + "mirna:" + input_file.at(i).identifier);
        struct_mapping_file.addLine( "mirna:" +String(input_file.at(i).identifier) + "\t" + entryformula.toString() + "\t" + entryformula.toString() + "\t" + entryformula.toString());
    }

    db_mapping_file.store(getStringOption_("out_db_mapping"));
    struct_mapping_file.store(getStringOption_("out_struct_mapping"));
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPNucleotideID tool;
  return tool.main(argc, argv);
}
