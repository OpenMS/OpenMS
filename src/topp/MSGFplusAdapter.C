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

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <map>
#include <cstddef>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MSGFplusAdapter MSGFplusAdapter

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

class TOPPMSGFPlusAdapter :
  public TOPPBase
{
public:
  TOPPMSGFPlusAdapter() :
    TOPPBase("MSGF+Adapter", "MS/MS database search using MSGF+.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {

		registerInputFile_("in", "<file>", "", "Input file");
		setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("database", "<file>", "", "FASTA file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
		registerInputFile_("d", "<file>", "", "Database file");
		registerInputFile_("MSGFplus_executable", "<executable>", "java -jar MSGFPlus.jar", "MSGF+ executable of the installation e.g. 'java - jar MSGFPlus.jar'", false);

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

    String msgf_executable("java -Xmx3500M -jar " + getStringOption_("MSGFPlus_executable"));
   
    //~ input_filename = "yeast_mitochondrium.mzML";

    //-------------------------------------------------------------
    // add hardcoded parameters
    //-------------------------------------------------------------
    String parameterMSGF = "-s " + inputfile_name + " -d" + db_name + " -t 20 -ti 0,1 -thread 2 -tda 0 -m 0 -inst 0 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 3 -n 1 -addFeatures 0 ";
    

    // write the msgf output file in the temporary directory
    String temp_directory = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }

    //-------------------------------------------------------------
    // execute MSGF+
    //-------------------------------------------------------------
    String msgfplus_output_filename(temp_directory + "msgfplus_output_file.mzid");
   
    //run MSGFPlus process and create the mzid file  
    int status = QProcess::execute((msgf_executable + " " + parameterMSGF + "-o " + msgfplus_output_filename).toQString());
    
    if (status != 0)
    {
      writeLog_("MSGF+ problem. Aborting! Calling command was: '" + msgf_executable + " \"" + inputfile_name + "\"'.\nDoes the MSGF+ executable exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // execute tsv converter
    //------------------------------------------------------------- 
		String mzidtotsv_output_filename(temp_directory + "svFile.tsv");
    String converter = "java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i " + msgfplus_output_filename + " -o "+mzidtotsv_output_filename+" -showQValue 1 -showDecoy 0 -unroll 1";
		
		status = QProcess::execute((msgf_executable + " " + parameterMSGF + "-o " + msgfplus_output_filename).toQString());
    
		if (status != 0)
    {
      writeLog_("MzIDToTSVConverter problem. Aborting! \nDoes the MzIDToTSVConverter executable exist?");
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

		IdXMLFile idxmlFile;
		idxmlFile.store(outputfile_name, protein_ids, peptide_ids);
		
	}	

};


int main(int argc, const char** argv)
{
  TOPPMSGFPlusAdapter tool;

  return tool.main(argc, argv);
}

  
  
