// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <map>
#include <iterator>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_PTPredict PTPredict
	
	@brief This application is used to predict the likelihood 
				 of peptides to be proteotypic.
	
	This method has been described in the publication 
	
	Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gr&ouml;pl, Oliver Kohlbacher and Knut Reinert LC-MSsim - a simulation software for Liquid ChromatographyMass Spectrometry data
  BMC Bioinformatics 2008, 9:423.

	The input of this application is an svm model and an IdXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the @ref TOPP_PTModel application.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PTPredict.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_PTPredict.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPTPredict
	: public TOPPBase
{
	public:
		TOPPPTPredict()
			: TOPPBase("PTPredict","predicts the likelihood of peptides to be proteotypic via svm_model which is trained by PTModel")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file ");
			setValidFormats_("in",StringList::create("idXML"));
			registerOutputFile_("out","<file>","","output file\n", false);
			setValidFormats_("out",StringList::create("idXML"));
			registerInputFile_("svm_model","<file>","","svm model in libsvm format (can be produced by PTModel)");
			registerIntOption_("max_number_of_peptides", "<int>",100000,"the maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).\n",false);
		}

		ExitCodes main_(int , const char**)
		{
			IdXMLFile IdXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			vector< String > peptides;
			vector<PeptideHit> temp_peptide_hits;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			vector<DoubleReal> predicted_likelihoods;
			vector<DoubleReal> predicted_labels;
			map< String, DoubleReal > predicted_data;
			svm_problem* training_data = NULL;
			svm_problem* prediction_data = NULL;
			UInt border_length = 0;
			UInt k_mer_length = 0;
			DoubleReal sigma = 0;
			String temp_string = "";
			UInt maximum_length = 50;
			String inputfile_name = "";
			String outputfile_name = "";
			String svmfile_name = "";
			UInt max_number_of_peptides = getIntOption_("max_number_of_peptides");

			
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
			inputfile_name = getStringOption_("in");			
			outputfile_name = getStringOption_("out");	

			svmfile_name = getStringOption_("svm_model");				

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			svm.loadModel(svmfile_name);

			// Since the POBK is not included in the libsvm we have to load
			// additional parameters from additional files.
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
			{
				inputFileReadable_(svmfile_name + "_additional_parameters", "svm_model (derived)");
	
				Param additional_parameters;
				
				additional_parameters.load(svmfile_name + "_additional_parameters");
				if (additional_parameters.getValue("kernel_type") != DataValue::EMPTY)
				{
					svm.setParameter(SVMWrapper::KERNEL_TYPE, ((String) additional_parameters.getValue("kernel_type")).toInt());
				}
								
				if (additional_parameters.getValue("border_length") == DataValue::EMPTY
						&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					writeLog_("No border length saved in additional parameters file. Aborting!");
					cout << "No border length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				border_length = ((String)additional_parameters.getValue("border_length")).toInt();
				if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
						&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					writeLog_("No k-mer length saved in additional parameters file. Aborting!");
					cout << "No k-mer length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
				if (additional_parameters.getValue("sigma") == DataValue::EMPTY
						&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					writeLog_("No sigma saved in additional parameters file. Aborting!");
					cout << "No sigma saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
				
			}				
			String document_id;			
			IdXML_file.load(inputfile_name, protein_identifications, identifications, document_id);
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for (Size i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].getHits();
				for (Size j = 0; j < temp_peptide_hits.size(); j++)
				{
					peptides.push_back(temp_peptide_hits[j].getSequence().toUnmodifiedString());
				}
			}
			
			vector<DoubleReal> labels;
			labels.resize(peptides.size(), 0);
			
			vector<String>::iterator it_from = peptides.begin();
			vector<String>::iterator it_to = peptides.begin();
			while(it_from != peptides.end())
			{
				vector<String> temp_peptides;
				vector<DoubleReal> temp_labels;
				UInt i = 0;
				while(i <= max_number_of_peptides && it_to != peptides.end())
				{
					++it_to;
					++i;
				}
								
				temp_peptides.insert(temp_peptides.end(), it_from, it_to);
				temp_labels.resize(temp_peptides.size(), 0);
				
				if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
				{
					prediction_data = 
						encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(temp_peptides,
																																				temp_labels,
																															 					allowed_amino_acid_characters,
																															 					maximum_length);
				}
				else if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(temp_peptides, 
																																							temp_labels, 
																																							k_mer_length, 
																																							allowed_amino_acid_characters, 
																																							border_length);				
				}
						
				if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					inputFileReadable_((svmfile_name + "_samples").c_str(), "svm_model (derived)");
	
					training_data = encoder.loadLibSVMProblem(svmfile_name + "_samples");
					svm.setTrainingSample(training_data);
	
					svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
					svm.setParameter(SVMWrapper::SIGMA, sigma);
				}
		    svm.getSVCProbabilities(prediction_data, predicted_likelihoods, predicted_labels);
	
				for (Size i = 0; i < temp_peptides.size(); i++)
				{
					predicted_data.insert(make_pair(temp_peptides[i],
																					(predicted_likelihoods[i])));
				}
				predicted_likelihoods.clear();
				predicted_labels.clear();
				LibSVMEncoder::destroyProblem(prediction_data);			

				it_from = it_to;
			}
			
			for (Size i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].getHits();
				for (Size j = 0; j < temp_peptide_hits.size(); j++)
				{
					DoubleReal temp_likelihood = predicted_data[temp_peptide_hits[j].getSequence().toUnmodifiedString()];

					temp_peptide_hits[j].setMetaValue("predicted_PT", temp_likelihood);
				}
				identifications[i].setHits(temp_peptide_hits);				
			}
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			IdXML_file.store(outputfile_name,
											 protein_identifications,
											 identifications);
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPPTPredict tool;
	return tool.main(argc,argv);
}
  
/// @endcond





