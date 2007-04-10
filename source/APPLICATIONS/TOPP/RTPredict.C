// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <map>

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
	
	@todo Store RT prediction p-values as identifications (Nico)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTPredict
	: public TOPPBase
{
	public:
		TOPPRTPredict()
			: TOPPBase("RTPredict","predicts retention times for peptides via the svm_model that is trained by RTModel")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>",""," input file in analysisXML format");
			registerStringOption_("out","<file>","","output file in analysisXML format");
			registerStringOption_("svm_model","<file>","","svm model in libsvm format (can be produced by RTModel)");
			registerDoubleOption_("total_gradient_time","<time>",1.0,"the time (in seconds) of the gradient");
		}

		ExitCodes main_(int , char**)
		{
			AnalysisXMLFile analysisXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<IdentificationData> identifications;
			vector< String > peptides;
			vector< DoubleReal > training_retention_times;
			vector<PeptideHit> temp_peptide_hits;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			vector<DoubleReal> predicted_retention_times;
			map< String, DoubleReal > predicted_data;
			svm_problem* training_data = NULL;
			svm_problem* prediction_data = NULL;
			UInt border_length = 0;
			UInt k_mer_length = 0;
			DoubleReal sigma = 0;
			DoubleReal sigma_0 = 0;
			DoubleReal sigma_max = 0;
			String temp_string = "";
			UInt maximum_length = 50;
			pair<DoubleReal, DoubleReal> temp_point;
			vector<Real> performance_retention_times;

			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
			String inputfile_name = getStringOption_("in");			
			inputFileReadable_(inputfile_name);
			String outputfile_name = getStringOption_("out");	
			outputFileWritable_(outputfile_name);
			String svmfile_name = getStringOption_("svm_model");
			inputFileReadable_(svmfile_name);			
			Real total_gradient_time = getDoubleOption_("total_gradient_time");

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			svm.loadModel(svmfile_name);

			// Since the oligo border kernel is not included in the libsvm we have to load
			// additional parameters from additional files.
			if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				inputFileReadable_(svmfile_name + "_additional_parameters");
	
				Param additional_parameters;
				
				additional_parameters.load(svmfile_name + "_additional_parameters");
				if (additional_parameters.getValue("kernel_type") != DataValue::EMPTY)
				{
					svm.setParameter(KERNEL_TYPE, ((String) additional_parameters.getValue("kernel_type")).toInt());
					cout << "Kernel type = " << svm.getIntParameter(KERNEL_TYPE) << endl;					
				}
								
				if (additional_parameters.getValue("border_length") == DataValue::EMPTY
						&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
				{
					writeLog_("No border length saved in additional parameters file. Aborting!");
					cout << "No border length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				border_length = ((String)additional_parameters.getValue("border_length")).toInt();
				if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
						&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
				{
					writeLog_("No k-mer length saved in additional parameters file. Aborting!");
					cout << "No k-mer length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
				if (additional_parameters.getValue("sigma") == DataValue::EMPTY
						&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
				{
					writeLog_("No sigma saved in additional parameters file. Aborting!");
					cout << "No sigma length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
				
				if (additional_parameters.getValue("sigma_0") == DataValue::EMPTY)
				{
					writeLog_("No sigma_0 saved in additional parameters file. Aborting!");
					cout << "No sigma_0 length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				sigma_0 = ((String)additional_parameters.getValue("sigma_0")).toFloat();
				if (additional_parameters.getValue("sigma_max") == DataValue::EMPTY)
				{
					writeLog_("No sigma_max saved in additional parameters file. Aborting!");
					cout << "No sigma_max length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				sigma_max = ((String)additional_parameters.getValue("sigma_max")).toFloat();
			}				
			
			analysisXML_file.load(inputfile_name, protein_identifications, identifications);
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			for(UInt i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].id.getPeptideHits();
				for(UInt j = 0; j < temp_peptide_hits.size(); j++)
				{
					peptides.push_back(temp_peptide_hits[j].getSequence());
				}
			}
			
			vector<DoubleReal> rts;
			rts.resize(peptides.size(), 0);
			if (svm.getIntParameter(KERNEL_TYPE) != OLIGO)
			{
				prediction_data = 
					encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(peptides,
																																			rts,
																														 					allowed_amino_acid_characters,
																														 					maximum_length);
			}
			else if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptides, 
																																						rts, 
																																						k_mer_length, 
																																						allowed_amino_acid_characters, 
																																						border_length);				
			}
					
			if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				inputFileReadable_((svmfile_name + "_samples").c_str());

				training_data = encoder.loadLibSVMProblem(svmfile_name + "_samples");
				cout << "Loading training_data" << endl;
				svm.setTrainingSample(training_data);

				svm.setParameter(BORDER_LENGTH, (Int) border_length);
				svm.setParameter(SIGMA, sigma);
				svm.predict(prediction_data, predicted_retention_times);
			}
			else
			{
				svm.predict(prediction_data, predicted_retention_times);
			}

			for(UInt i = 0; i < peptides.size(); i++)
			{
				predicted_data.insert(make_pair(peptides[i],
																				(predicted_retention_times[i] * total_gradient_time)));
			}
		
			for(UInt i = 0; i < identifications.size(); i++)
			{
				temp_peptide_hits = identifications[i].id.getPeptideHits();
				vector<ProteinHit> temp_protein_hits = identifications[i].id.getProteinHits();
				for(UInt j = 0; j < temp_peptide_hits.size(); j++)
				{
					DoubleReal temp_rt = predicted_data[temp_peptide_hits[j].getSequence()];

					temp_point.first = identifications[i].rt;
					temp_point.second = temp_rt;
					DoubleReal temp_p_value = svm.getPValue(sigma_0, sigma_max, temp_point);
					temp_peptide_hits[j].setPredictedRTPValue(temp_p_value);

					performance_retention_times.push_back(identifications[i].rt);					
				}
				identifications[i].id.setPeptideAndProteinHits(temp_peptide_hits,
																											temp_protein_hits);				
			}
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			analysisXML_file.store(outputfile_name,
														 protein_identifications,
														 identifications,
														 predicted_data);
			writeDebug_("Linear correlation between predicted and measured rt is: "
									+ String(Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(predicted_retention_times.begin(), 
																		predicted_retention_times.end(), 
																		performance_retention_times.begin(), 
																		performance_retention_times.end())), 1);														 
			writeDebug_("MSE between predicted and measured rt is: "
									+ String(Math::BasicStatistics<Real>::meanSquareError(predicted_retention_times.begin(), 
																		predicted_retention_times.end(), 
																		performance_retention_times.begin(), 
																		performance_retention_times.end())), 1);														 
			writeDebug_("Linear correlation between predicted and measured rt is: "
									+ String(Math::BasicStatistics<Real>::pearsonCorrelationCoefficient(predicted_retention_times.begin(), 
																		predicted_retention_times.end(), 
																		performance_retention_times.begin(), 
																		performance_retention_times.end())), 1);														 
			writeDebug_("MSE between predicted and measured rt is: "
									+ String(Math::BasicStatistics<Real>::meanSquareError(predicted_retention_times.begin(), 
																		predicted_retention_times.end(), 
																		performance_retention_times.begin(), 
																		performance_retention_times.end())), 1);														 
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPRTPredict tool;
	return tool.main(argc,argv);
}
  
/// @endcond





