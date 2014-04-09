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
// $Maintainer: Mathias Walzer $
// $Authors: Dilek Dere, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <OpenMS/FORMAT/TextFile.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <map>
#include <cstddef>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MSGFPlusAdapter MSGFPlusAdapter

   @brief 
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MSGFplusAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em MSGF+ must be installed before this wrapper can be used. This wrapper
    has been successfully tested. Please be sure that java and MSGF+ are working.

    This adapter supports relative database filenames, which (when not found in the 
    current working directory) is looked up in the directories specified by 
    'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).
		
    The adapter has four input parameters. The input file is a spectrum file, the database
    is the used database file for the search usually as fasta file,
    the name of the output file as idXML and the MSGF+ executable.
    Other parameters are set in the implementation and cannot be changed by the user
    in the moment. 
    First using the input spectrum file and the default parameters MSGF+ is started.
    MSGF+ muss be installed in the same folder as the executable MSGF+Adapter. 
    The output of the MSGF+Adapter is stored in a temporary directory using the file name
    "msgfplus_output_file.mzid".
    This file is then used and converted in a second QProcess run into a tsv file.
    For this part the process "java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv" is used 
    and the parameters are again fixed in the implementation yet.
    In the last step the created tsv file is parsed in and an idXML file is created.  

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_MSGFplusAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MSGFplusAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;
using namespace std;

class MSGFPlusAdapter :
  public TOPPBase
{
public:
  MSGFPlusAdapter() :
    TOPPBase("MSGFPlusAdapter", "MS/MS database search using MSGFPlus.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML,mgf,ms2"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerInputFile_("database", "<file>", "", "FASTA file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));

    registerInputFile_("MSGFPlus_executable", "<executable>", "", "MSGFPlus executable of the installation e.g. 'java - jar MSGFPlus.jar'");

    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 20, "Precursor mono mass tolerance.", false);
    registerStringOption_("precursor_mass_tolerance_unit", "<unit>", "ppm", "Unit to be used for precursor mass tolerance.", false);
    setValidStrings_("precursor_mass_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerStringOption_("isotope_error_range", "<range>", "0,1", "Range of allowed isotope peak errors. Takes into account of the error introduced by chooosing a non-monoisotopic peak for fragmentation. The combination of -t and -ti determins the precursor mass tolerance. E.g. -precursor_mass_tolerance 20 -precursor_mass_tolerance_unit ppm -isotope_error_range -1,2 tests abs(exp-calc-n*1.00335Da)<20ppm for n=-1, 0, 1, 2", false);

    registerIntOption_("decoy", "<0/1>", 0, "0: don't search decoy database, 1: search decoy database", false);
    registerIntOption_("fragment_method", "<method>", 0, "0: As written in the spectrum or CID if no info, 1: CID, 2: ETD, 3: HCD", false);
    registerIntOption_("instrument", "<instrument>", 0, "0: Low-res LCQ/LTQ, 1: High-res LTQ, 2: TOF, 3: Q-Exactive", false);
    registerIntOption_("enzyme", "<enzyme>", 1, "0: unspecific cleavage, 1: Trypsin, 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: glutamyl endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage", false);
    registerIntOption_("protocol", "<protocol>", 0, "0: NoProtocol, 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho, 4: TMT", false);    
    registerIntOption_("tolerable_termini", "<num>", 2, "For trypsin, 0: non-tryptic, 1: semi-tryptic, 2: fully-tryptic peptides only.", false);
    
    registerInputFile_("mod", "<file>", "", "Modification configuration file", false);

    registerIntOption_("min_precursor_charge", "<charge>", 2, "Minimum precursor ion charge", false);
    registerIntOption_("max_precursor_charge", "<charge>", 3, "Maximum precursor ion charge", false);

    registerIntOption_("min_peptide_length", "<length>", 6, "Minimum peptide length to consider", false);
    registerIntOption_("max_peptide_length", "<length>", 40, "Maximum peptide length to consider", false);   

    registerIntOption_("matches_per_spec", "<num>", 1, "Number of matches per spectrum to be reported", false);
    registerIntOption_("add_features", "<num>", 0, "0: output basic scores only, 1: output additional features", false);
  }

  //The following sequence modification methods are used to modify the sequence stored in the tsv such that it can be used by AASequence

  //Method to cut the first and last character of the sequence.
  //The sequences in the tsv file has the form K.AAAA.R (AAAA stands for any amino acid sequence.
  //After this method is used the sequence AAAA results
  String cutSequence (String sequence)
  {
    String modifiedSequence = sequence;

    //search for the first and last occurence of .
    std::size_t findFirst = sequence.find_first_of(".");
    std::size_t findLast = sequence.find_last_of(".");
       
    //used the found positions and cut the sequence 
    if (findFirst!=std::string::npos && findLast!=std::string::npos && findFirst != findLast)
    {
      modifiedSequence = sequence.substr(findFirst+1, findLast-2);
    }
		
    return modifiedSequence;
  }

  //Method to replace comma by point.
  //This is used as point should be used as separator of decimals instead of comma
  String changeKomma (String seq)
  {
    std::size_t found = seq.find_first_of(".,");
    while (found!=std::string::npos) 
    {
      seq[found]='.';
      found=seq.find_first_of(".,",found+1);
    }
    return seq;
  }
	
  //Method to replace the mass representation of modifications.
  //Modifications in the tsv file has the form M+15.999 e.g.
  //After using this method the sequence should look like this: M[+15.999] 
  String modifySequence (String seq)
  {
    String modifiedSequence = seq;
	std::size_t found = modifiedSequence.find_first_of("+-");
    while (found!=std::string::npos)
    {
      modifiedSequence = modifiedSequence.insert(found, "[");
      std::size_t found1 = modifiedSequence.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", found);
      if (found1!=std::string::npos)
      {
        modifiedSequence.insert(found1, "]");
        found = modifiedSequence.find_first_of("+-", found1+2);
      } 
      else 
      { //if last amino acid is modified
        modifiedSequence = modifiedSequence + "]";
        return modifiedSequence;
      }
    }
    return modifiedSequence;
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name == "")
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String outputfile_name = getStringOption_("out");
    writeDebug_(String("Output file: ") + outputfile_name, 1);
    if (outputfile_name == "")
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String db_name(getStringOption_("database"));
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }

    // write the msgf output file in the temporary directory
    String temp_directory = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }

    String msgfplus_output_filename(temp_directory + "msgfplus_output_file.mzid");
    String  parameters = "";
    parameters += "-s " + inputfile_name;
    parameters += " -o " + msgfplus_output_filename;
    parameters += " -d " + db_name;
    parameters += " -t " + String(getDoubleOption_("precursor_mass_tolerance")) + getStringOption_("precursor_mass_tolerance_unit");
    parameters += " -ti " + getStringOption_("isotope_error_range");
    parameters += " -tda " + String(getIntOption_("decoy"));
    parameters += " -m " + String(getIntOption_("fragment_method"));
    parameters += " -inst " + String(getIntOption_("instrument"));
    parameters += " -e " + String(getIntOption_("enzyme"));
    parameters += " -protocol " + String(getIntOption_("protocol"));
    parameters += " -ntt " + String(getIntOption_("tolerable_termini"));
    parameters += " -minLength " + String(getIntOption_("min_peptide_length"));
    parameters += " -maxLength " + String(getIntOption_("max_peptide_length"));
    parameters += " -minCharge " + String(getIntOption_("min_precursor_charge"));
    parameters += " -maxCharge " + String(getIntOption_("max_precursor_charge"));
    parameters += " -n " + String(getIntOption_("matches_per_spec"));
    parameters += " -addFeatures " + String(getIntOption_("add_features"));

    String modfile_name = getStringOption_("mod");
    if(modfile_name != "") {
      parameters += " -mod " + getStringOption_("mod");
    }

    cout << "Parameters: " << parameters << endl;

    //-------------------------------------------------------------
    // execute MSGFPlus
    //-------------------------------------------------------------
   
    //run MSGFPlus process and create the mzid file  
    String msgf_executable("java -Xmx3500M -jar " + getStringOption_("MSGFPlus_executable"));

    QProcess process;
    int status = process.execute((msgf_executable + " " + parameters).toQString());
    if (status != 0)
    {
      writeLog_("MSGFPlus problem. Aborting! Calling command was: '" + msgf_executable + " \"" + inputfile_name + "\"'.\nDoes the MSGFPlus executable exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // execute tsv converter
    //------------------------------------------------------------- 
    String mzidtotsv_output_filename(temp_directory + "svFile.tsv");
    String converter_executable("java -cp " + getStringOption_("MSGFPlus_executable") + " edu.ucsd.msjava.ui.MzIDToTsv ");

    parameters = "-i " +  msgfplus_output_filename;
    parameters += " -o " + mzidtotsv_output_filename;
    parameters += " -showQValue 1";
    parameters += " -showDecoy 0";
    parameters += " -unroll 1";

    status = process.execute((converter_executable + " " + parameters).toQString());
    if (status != 0)
    {
      writeLog_("MzIDToTSVConverter problem. Aborting! \nDoes the MzIDToTSVConverter executable exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // create idXML
    //------------------------------------------------------------- 
    CsvFile tsvfile;
    tsvfile.load(mzidtotsv_output_filename, "\t");

    //create idXML file
    vector<ProteinIdentification> protein_ids(1);
    ProteinIdentification protein_id;

    //store all peptide hits in a map, the key is the scannumber
    map<int,PeptideIdentification> peptide_hits;

    DoubleReal score; //use SpecEValue from the tsv file
    UInt rank; 
    Int charge;
    AASequence sequence;
    int scanNumber;
    double precursor_mz;

    //iterate over the rows of the tsv file
    for (CsvFile::Iterator it = tsvfile.begin() + 1 ; it != tsvfile.end(); ++it)
    {
      vector<String> elements; 
      it->split("\t", elements);
      scanNumber = elements[2].toInt();
      score = elements[12].toDouble();
      rank = 0; //set to 0 in the moment
      charge = elements[7].toInt();
      sequence = AASequence(modifySequence(changeKomma(cutSequence(elements[8])))); //sequence must be cutted and modified
      PeptideHit p_hit(score, rank, charge, sequence);
      precursor_mz = elements[4].toDouble(); 
      peptide_hits[scanNumber].insertHit(p_hit);
      peptide_hits[scanNumber].setMetaValue("MZ", precursor_mz);
      peptide_hits[scanNumber].setMetaValue("ScanNumber", scanNumber);
      peptide_hits[scanNumber].setScoreType("SpecEValue");
      peptide_hits[scanNumber].setHigherScoreBetter(false);
      peptide_hits[scanNumber].setMetaValue("RT", 0);   //ToDo, retention time is not given in the tsv file yet is set to 0 at the moment
    }

    //iterate over map and create a vector<PeptideIdentification>
    map<int,PeptideIdentification>::iterator it;
    vector<PeptideIdentification> peptide_ids;
    PeptideIdentification pep;
    for(map<int,PeptideIdentification>::iterator it=peptide_hits.begin(); it!=peptide_hits.end(); ++it)
    {
      pep = it->second;
      pep.sort();
      peptide_ids.push_back(pep);
    }
    
    cout << "outputfile_name: " << outputfile_name << endl;
    IdXMLFile().store(outputfile_name, protein_ids, peptide_ids);
    return EXECUTION_OK;
  }	
};


int main(int argc, const char** argv)
{
  MSGFPlusAdapter tool;

  return tool.main(argc, argv);
}

  
  
