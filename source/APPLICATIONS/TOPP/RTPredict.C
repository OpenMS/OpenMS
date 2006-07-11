// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/MATH/STATISTICS/EvaluationFunctions.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include "TOPPBase.h"

#include <qfileinfo.h>
#include <qfile.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page RTPredict RTPredict
	
	@brief This application is used to predict retention times for peptides
	
	The input of this component is a svm model and an analysisXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the RTModel component. 
	<br>
	The peptide sequences are extracted from the analysisXML inputfile 
	and passed to the svm. The svm then predicts retention times
	according to the trained model. The predicted retention times
	are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" /> 
	@endcode inside the peptide entities in the analysisXML output file.

	@todo Fix --help and --help-opt output (Nico)
	
	@ingroup TOPP
*/



// We do not want this class to show up in the docu -> @cond
/// @cond 

class TOPPRTPredict
	: public TOPPBase
{
	public:
		TOPPRTPredict()
			: TOPPBase("RTPredict")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			cerr << endl
       << tool_name_ << " -- Predicts retention times for peptides"
       << " via the svm_model that is trained by RTModel."
       << endl
       << "Usage:" << endl
			 << " " << tool_name_ << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>   			 input file in analysisXML format (default read from INI file)" << endl
			 << "  -svm_model <file>   		 svm model in libsvm format (can be produced by RTModel) "
			 << "  -total_gradient_time <file> the time (in seconds) of the gradient "
			 << "(default read from INI file)" << endl
			 << "  -out <file>  			 output file in analysisXML format (default read from INI file)" << endl
			 << endl ;
		}

	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-svm_model"] = "svm_model";
			options_["-total_gradient_time"] = "total_gradient_time";
			options_["-ini"] = "ini";
			options_["-log"] = "log";
			options_["-n"] = "instance";
			options_["-d"] = "debug";
			options_["--help"] = "help";
		}

		void printToolHelpOpt_()
		{
			cerr << endl
       << tool_name_ << " -- Predicts retention times for peptides"
       << " via the svm_model that is trained by RTModel."
       << endl
       << "Usage:" << endl
			 << " " << tool_name_ << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>   			 input file in analysisXML format (default read from INI file)" << endl
			 << "  -svm_model <file>   		 svm model in libsvm format (can be produced by RTModel) "
			 << "  -total_gradient_time <file> the time (in seconds) of the gradient "
			 << "(default read from INI file)" << endl
			 << "  -out <file>  			 output file in analysisXML format (default read from INI file)" << endl
			 << endl
			 << "Common TOPP options are:" << endl
			 << "  -ini <file>       TOPP INI file (default: TOPP.ini)" << endl
			 << "  -log <file>       log file (default: TOPP.log)" << endl
			 << "  -n <int>          instance number (default: 1)" << endl
			 << "  -d <level>        sets debug level (default: 0)" << endl
			 << "  --help            shows this help" << endl
       << "  --help-opt        shows help on the INI options accepted" << endl
			 << endl ;
		}	

		ExitCodes main_(int , char**)
		{
			QFileInfo file_info;
			QFile file;
			String inputfile_name;
			String svmfile_name;
			String outputfile_name;
			AnalysisXMLFile analysisXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<Identification> identifications;
			vector<float> precursor_retention_times;
			vector<float> precursor_mz_values;
			ContactPerson contact_person;
			vector< String > peptides;
			vector< DoubleReal > training_retention_times;
			vector<PeptideHit> temp_peptide_hits;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			vector< vector< pair<UnsignedInt, DoubleReal> > >* encoded_composition_vectors;
			vector<svm_node*>* encoded_LIBSVM_vectors;
			vector<DoubleReal>* predicted_retention_times;
			map< String, DoubleReal > predicted_data;
			Real total_gradient_time = 1;
  

		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
		
			//input file names and types
			inputfile_name = getParamAsString_("in");			
			writeDebug_(String("Input file: ") + inputfile_name, 1);
			if (inputfile_name == "")
			{
				writeLog_("No input file specified. Aborting!");
				cout << "No input file specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			
			//output file names and types
			outputfile_name = getParamAsString_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);
			if (outputfile_name == "")
			{
				writeLog_("No output file specified. Aborting!");
				cout << "No output file specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}				
			
			total_gradient_time = getParamAsString_("total_gradient_time", "0.f").toFloat();
			writeDebug_(String("Total gradient time: ") + String(total_gradient_time), 1);
			if (total_gradient_time == 0.f)
			{
				writeLog_(String("Total gradient time has to") + 
						String(" be specified. Aborting!"));
				cout << "Total gradient time has to"
						<< " be specified. Aborting!" << endl;
				return ILLEGAL_PARAMETERS;
			}				

			//svm model file
			svmfile_name = getParamAsString_("svm_model");			
			writeDebug_(String("SVM model file: ") + svmfile_name, 1);
			if (svmfile_name == "")
			{
				writeLog_("No svm model file specified. Aborting!");
				cout << "No svm model file specified. Aborting!" << endl;
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			//-------------------------------------------------------------
			// testing whether input and output files are accessible
			//-------------------------------------------------------------
			
			file_info.setFile(inputfile_name.c_str());
			if (!file_info.exists())
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
			}
			if (!file_info.isReadable())
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);			
			}
			if (file_info.size() == 0)
			{
				throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
			}		
			file_info.setFile(svmfile_name.c_str());
			if (!file_info.exists())
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, svmfile_name);
			}
			if (!file_info.isReadable())
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, svmfile_name);			
			}
			file.setName(outputfile_name.c_str());
			file.open( IO_WriteOnly );
			if (!file.isWritable())
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, outputfile_name);
			}
			file.close();				


			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			analysisXML_file.load(inputfile_name,
														&protein_identifications,
														&identifications,
														&precursor_retention_times,
														&precursor_mz_values,
														&contact_person);
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for(UnsignedInt i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].getPeptideHits();
				for(UnsignedInt j = 0; j < temp_peptide_hits.size(); j++)
				{
					peptides.push_back(temp_peptide_hits[j].getSequence());
				}
			}
			
			encoded_composition_vectors = 
				encoder.encodeCompositionVectors(peptides,
																				 allowed_amino_acid_characters);
			encoded_LIBSVM_vectors = 
				encoder.encodeLIBSVMVectors(*encoded_composition_vectors);
		
			svm.loadModel(svmfile_name);	
			predicted_retention_times = 
				svm.predict(*encoded_LIBSVM_vectors);
		
			delete encoded_composition_vectors;
			delete encoded_LIBSVM_vectors;
		
			for(UnsignedInt i = 0; i < peptides.size(); i++)
			{
				predicted_data.insert(make_pair(peptides[i],
																				((*predicted_retention_times)[i] * total_gradient_time)));
			}
		
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			analysisXML_file.store(outputfile_name,
														 protein_identifications,
														 identifications,
														 precursor_retention_times,
														 precursor_mz_values,
														 contact_person,
														 predicted_data,
														 svm.getSVRProbability());
			return OK;
		}
};

///@endcond


int main( int argc, char ** argv )
{
	TOPPRTPredict tool;
	return tool.main(argc,argv);
}
  




