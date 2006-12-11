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
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

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
	
	The input of this application is a svm model and an analysisXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the RTModel application. 
	<br>
	The peptide sequences are extracted from the analysisXML inputfile 
	and passed to the svm. The svm then predicts retention times
	according to the trained model. The predicted retention times
	are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" /> 
	@endcode inside the peptide entities in the analysisXML output file.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTPredict
	: public TOPPBase
{
	public:
		TOPPRTPredict()
			: TOPPBase("RTPredict")
		{
			
		}
	
	protected:
		void printToolUsage_() const
		{
			cerr << endl
       << getToolName() << " -- Predicts retention times for peptides via the svm_model that is trained by RTModel." << endl
       << "Version: " << VersionInfo::getVersion() << endl
       << endl
       << "Usage:" << endl
			 << " " << getToolName() << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>              input file in analysisXML format (default read from INI file)" << endl
			 << "  -svm_model <file>       svm model in libsvm format (can be produced by RTModel) " << endl
			 << "  -total_gradient_time    the time (in seconds) of the gradient "
			 << "(default read from INI file)" << endl
			 << "  -out <file>             output file in analysisXML format (default read from INI file)" << endl
			 << endl ;
		}

	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-svm_model"] = "svm_model";
			options_["-total_gradient_time"] = "total_gradient_time";
		}

		void printToolHelpOpt_() const
		{
			cerr << endl
		       << getToolName() << endl
		       << endl
		       << "INI options:" << endl
					 << "  in                        input file" << endl
					 << "  out                       output file" << endl
					 << "  svm_model                 svm model in libsvm format (can be produced by RTModel) " << endl
					 << "  total_gradient_time       the time (in seconds) of the gradient "
					 << endl << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.analysisXML\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"output.analysisXML\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"svm_model\" value=\"model.svm\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"total_gradient_time\" value=\"3000\" type=\"float\"/>" << endl;
 					 
		}	

		ExitCodes main_(int , char**)
		{
			String inputfile_name;
			String svmfile_name;
			String outputfile_name;
			AnalysisXMLFile analysisXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<IdentificationData> identifications;
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
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}
			
			//output file names and types
			outputfile_name = getParamAsString_("out");
			writeDebug_(String("Output file: ") + outputfile_name, 1);
			if (outputfile_name == "")
			{
				writeLog_("No output file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}				
			
			total_gradient_time = getParamAsString_("total_gradient_time", "0.f").toFloat();
			writeDebug_(String("Total gradient time: ") + String(total_gradient_time), 1);
			if (total_gradient_time == 0.f)
			{
				writeLog_("Total gradient time has to be specified. Aborting!");
				return ILLEGAL_PARAMETERS;
			}				

			//svm model file
			svmfile_name = getParamAsString_("svm_model");			
			writeDebug_(String("SVM model file: ") + svmfile_name, 1);
			if (svmfile_name == "")
			{
				writeLog_("No svm model file specified. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;
			}

			//-------------------------------------------------------------
			// testing whether input and output files are accessible
			//-------------------------------------------------------------
			
			if (!File::exists(inputfile_name))
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
			}
			if (!File::readable(inputfile_name))
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);			
			}
			if (File::empty(inputfile_name))
			{
				throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
			}		
			

			if (!File::exists(svmfile_name))
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, svmfile_name);
			}
			if (!File::readable(svmfile_name))
			{
				throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, svmfile_name);			
			}
			if (!File::readable(outputfile_name))
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, outputfile_name);
			}
			

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			analysisXML_file.load(inputfile_name,
														protein_identifications,
														identifications);
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for(UnsignedInt i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].id.getPeptideHits();
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
														 predicted_data,
														 svm.getSVRProbability());
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPRTPredict tool;
	return tool.main(argc,argv);
}
  
/// @endcond





