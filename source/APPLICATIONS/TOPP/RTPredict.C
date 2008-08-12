// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page RTPredict RTPredict
	
	@brief This application is used to predict retention times 
				 for peptides or peptide separation.
	
	The input of this application 
	is an svm model and an IdXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the RTModel application. 
	<br>
	For retention time prediction the peptide sequences are extracted 
	from the IdXML inputfile 
	and passed to the svm. The svm then predicts retention times
	according to the trained model. The predicted retention times
	are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" /> 
	@endcode inside the peptide entities in the IdXML output file.
	For separation prediction you have to specify two output file names.
	'out_positive' is the filename of the peptides which are predicted
	to be collected by the column and 'out_negative' is the file
	of the predicted flowthrough peptides.
	
	@todo Reactivate TOPPtest after fixing issue below on MinGW (Niko, Chris)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTPredict
	: public TOPPBase
{
	public:
		TOPPRTPredict()
			: TOPPBase("RTPredict","Predicts retention times for peptides using a model trained by RTModel.")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file ");
			setValidFormats_("in",StringList::create("IdXML"));
			registerOutputFile_("out","<file>","","output file (peptide RT prediction)\n", false);
			setValidFormats_("out",StringList::create("IdXML"));
			registerFlag_("textfile_input", "if this flag is set, RTPredict expects a textfile instead of an IdXML file as input which contains one peptide sequence per line output as a textfile is switched on as well");
			registerFlag_("textfile_output", "if this flag is set, RTPredict just writes a peptide sequence with the corresponding predicted retention time per line");
			registerOutputFile_("out_positive","<file>","","output file in IdXML format containing positive predictions (peptide separation prediction)\n", false);
			setValidFormats_("out_positive",StringList::create("IdXML"));
			registerOutputFile_("out_negative","<file>","","output file in IdXML format containing negative predictions (peptide separation prediction)\n", false);
			setValidFormats_("out_negative",StringList::create("IdXML"));
			registerInputFile_("svm_model","<file>","","svm model in libsvm format (can be produced by RTModel)");
			registerDoubleOption_("total_gradient_time","<time>",1.0,"the time (in seconds) of the gradient (peptide RT prediction)", false);
			setMinFloat_("total_gradient_time", 0.00001);
		}

		void loadStrings_(String filename, std::vector<String>& sequences)
		{
		    TextFile text_file(filename.c_str(), true);
		    TextFile::iterator it;
		    
		    sequences.clear();
		    
		    it = text_file.begin();	
		    while(it != text_file.end())
		    {
			  	sequences.push_back((*it).trim());
			  	++it;
		    } 	
		}
		
		void writeStringLabelLines_(String filename, map< String, DoubleReal > predicted_data)
		{
				ofstream os;
				map< String, DoubleReal >::iterator it;
				
				os.open(filename.c_str(), ofstream::out);
				
				for(it = predicted_data.begin(); it != predicted_data.end(); ++it)
				{
					os << it->first << " " << it->second << endl;
				}
				os.flush();
				os.close();																					
		}

		ExitCodes main_(int , const char**)
		{
			IdXMLFile IdXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			vector< String > peptides;
			vector< AASequence > modified_peptides;
			vector< DoubleReal > training_retention_times;
			vector<PeptideHit> temp_peptide_hits;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			vector<DoubleReal> predicted_retention_times;
			map< String, DoubleReal > predicted_data;
			map< AASequence, DoubleReal > predicted_modified_data;
			svm_problem* prediction_data = NULL;
			SVMData training_samples;
			SVMData prediction_samples;
			UInt border_length = 0;
			UInt k_mer_length = 0;
			DoubleReal sigma = 0;
			DoubleReal sigma_0 = 0;
			DoubleReal sigma_max = 0;
			String temp_string = "";
			UInt maximum_length = 50;
			pair<DoubleReal, DoubleReal> temp_point;
			vector<Real> performance_retention_times;
			String inputfile_name = "";
			String outputfile_name = "";
			String outputfile_name_positive = "";
			String outputfile_name_negative = "";
			String svmfile_name = "";
			Real total_gradient_time = 1.f;
			bool separation_prediction = false;
			vector<PeptideIdentification> identifications_positive;
			vector<PeptideIdentification> identifications_negative;
			bool textfile_input = false;
			bool textfile_output = false;
			bool first_dim_rt = false;
			UInt number_of_peptides = 0;
			
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
			inputfile_name = getStringOption_("in");			
			outputfile_name_positive = getStringOption_("out_positive");
			if (outputfile_name_positive != "")
			{
				outputfile_name_negative = getStringOption_("out_negative");
				if (outputfile_name_negative != "")
				{
					separation_prediction = true;					
				}
				else
				{
					writeLog_("No file name given for negative output . Aborting!");
					return ILLEGAL_PARAMETERS;					
				}
			}
			else
			{
				outputfile_name = getStringOption_("out");	
			}
			textfile_output = getFlag_("textfile_output");			
			textfile_input = getFlag_("textfile_input");			
			if (textfile_input)
			{
				textfile_output = true;
			}

			svmfile_name = getStringOption_("svm_model");	
			total_gradient_time = getDoubleOption_("total_gradient_time");

			

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			svm.loadModel(svmfile_name);
						
			if ((svm.getIntParameter(SVM_TYPE) == C_SVC || svm.getIntParameter(SVM_TYPE) == NU_SVC) && !separation_prediction)
			{
					writeLog_("You cannot perform peptide separation prediction with a model trained for"
										+ String("\npeptide retention time prediction. Aborting!"));
					return ILLEGAL_PARAMETERS;					
			}
			if ((svm.getIntParameter(SVM_TYPE) != C_SVC && svm.getIntParameter(SVM_TYPE) != NU_SVC) && separation_prediction)
			{
					writeLog_("You cannot perform peptide retention time prediction with a model trained for\n"
										+ String("peptide separation prediction. Aborting!"));
					return ILLEGAL_PARAMETERS;					
			}

			// Since the POBK is not included in the libsvm we have to load
			// additional parameters from additional files.
			if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				inputFileReadable_(svmfile_name + "_additional_parameters");
	
				Param additional_parameters;
				
				additional_parameters.load(svmfile_name + "_additional_parameters");
				if (additional_parameters.exists("first_dim_rt") 
						&& additional_parameters.getValue("first_dim_rt") != DataValue::EMPTY)
				{
					first_dim_rt = additional_parameters.getValue("first_dim_rt").toBool();
				}
				if (additional_parameters.getValue("kernel_type") != DataValue::EMPTY)
				{
					svm.setParameter(KERNEL_TYPE, ((String) additional_parameters.getValue("kernel_type")).toInt());
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
					cout << "No sigma saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
				
				if (!separation_prediction && additional_parameters.getValue("sigma_0") == DataValue::EMPTY)
				{
					writeLog_("No sigma_0 saved in additional parameters file. Aborting!");
					cout << "No sigma_0 length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				if (!separation_prediction && additional_parameters.getValue("sigma_0") != DataValue::EMPTY)
				{
					sigma_0 = additional_parameters.getValue("sigma_0");
				}
				if (!separation_prediction && additional_parameters.getValue("sigma_max") == DataValue::EMPTY)
				{
					writeLog_("No sigma_max saved in additional parameters file. Aborting!");
					cout << "No sigma_max length saved in additional parameters file. Aborting!" << endl;
					return ILLEGAL_PARAMETERS;					
				}
				if (!separation_prediction && additional_parameters.getValue("sigma_max") != DataValue::EMPTY)
				{
					sigma_max = additional_parameters.getValue("sigma_max");
				}
			}				
			
			if (textfile_input)
			{
				loadStrings_(inputfile_name, peptides);
				if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
				{
					for(UInt i = 0; i < peptides.size(); ++i)
					{
						modified_peptides.push_back(AASequence(peptides[i]));
					}
					peptides.clear();
				}
			}
			else
			{
				IdXML_file.load(inputfile_name, protein_identifications, identifications);
			}
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			if (!textfile_input)
			{
				for(UInt i = 0; i < identifications.size(); i++)
				{
					temp_peptide_hits = identifications[i].getHits();
					for(UInt j = 0; j < temp_peptide_hits.size(); j++)
					{
						if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
						{
							modified_peptides.push_back(temp_peptide_hits[j].getSequence());
						}
						else
						{
							peptides.push_back(temp_peptide_hits[j].getSequence().toUnmodifiedString());
						}
					}
				}
			}
			if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				number_of_peptides = modified_peptides.size();
			}
			else
			{
				number_of_peptides = peptides.size();
			}			
										
			vector<DoubleReal> rts;
			rts.resize(number_of_peptides, 0);
			
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
				encoder.encodeProblemWithOligoBorderVectors(modified_peptides, 
																										k_mer_length, 
																										allowed_amino_acid_characters, 
																										border_length,
																										prediction_samples.sequences);
				prediction_samples.labels = rts;																																	
			}
					
			if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				inputFileReadable_((svmfile_name + "_samples").c_str());

				training_samples.load(svmfile_name + "_samples");
				svm.setTrainingSample(training_samples);				

				svm.setParameter(BORDER_LENGTH, (Int) border_length);
				svm.setParameter(SIGMA, sigma);
				svm.predict(prediction_samples, predicted_retention_times);
			}
			else
			{
				svm.predict(prediction_data, predicted_retention_times);
				LibSVMEncoder::destroyProblem(prediction_data);
			}

			for(UInt i = 0; i < number_of_peptides; i++)
			{
				if (svm.getIntParameter(KERNEL_TYPE) == OLIGO && !textfile_output)
				{
					predicted_modified_data.insert(make_pair(modified_peptides[i],
																					(predicted_retention_times[i] * total_gradient_time)));
				}
				else if (svm.getIntParameter(KERNEL_TYPE) != OLIGO)
				{
					predicted_data.insert(make_pair(peptides[i],
																					(predicted_retention_times[i] * total_gradient_time)));
				}
				else
				{
					predicted_data.insert(make_pair(modified_peptides[i].toString(),
																					(predicted_retention_times[i] * total_gradient_time)));
				}
			}
			if (!textfile_input)
			{
				if (!separation_prediction)
				{
					for(UInt i = 0; i < identifications.size(); i++)
					{
						temp_peptide_hits = identifications[i].getHits();
						for(UInt j = 0; j < temp_peptide_hits.size(); j++)
						{
							DoubleReal temp_rt = 0.;
							DoubleReal temp_p_value = 0.;
														
							if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
							{
								temp_rt = predicted_modified_data[temp_peptide_hits[j].getSequence()];
							}
							else
							{
								temp_rt = predicted_data[temp_peptide_hits[j].getSequence().toUnmodifiedString()];
							}
							
							if (first_dim_rt)
							{
								temp_point.first = identifications[i].getMetaValue("first_dim_rt");
							}
							else
							{
								temp_point.first = identifications[i].getMetaValue("RT");
							}
							if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
							{
								temp_point.second = temp_rt;
								temp_p_value = svm.getPValue(sigma_0, sigma_max, temp_point);
							}
							if (first_dim_rt)
							{
								if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
								{
									temp_peptide_hits[j].setMetaValue("predicted_RT_p_value_first_dim",temp_p_value);
								}
								temp_peptide_hits[j].setMetaValue("predicted_RT_first_dim",temp_rt);
							}
							else
							{
								if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
								{
									temp_peptide_hits[j].setMetaValue("predicted_RT_p_value",temp_p_value);
								}
								temp_peptide_hits[j].setMetaValue("predicted_RT",temp_rt);
							}
							performance_retention_times.push_back(identifications[i].getMetaValue("RT"));					
						}
						identifications[i].setHits(temp_peptide_hits);				
					}
				}
				else
				{
					vector<PeptideHit> hits_positive;
					vector<PeptideHit> hits_negative;
					
					PeptideIdentification temp_identification;
	
					for(UInt i = 0; i < identifications.size(); i++)
					{					
						hits_negative.clear();
						hits_positive.clear();
	
						temp_peptide_hits = identifications[i].getHits();
						for(vector<PeptideHit>::iterator it = temp_peptide_hits.begin();
								it != temp_peptide_hits.end();
								++it)
						{						
							if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
							{
								if (predicted_modified_data[it->getSequence()] > 0)
								{
									hits_positive.push_back(*it);
								}
								else
								{
									hits_negative.push_back(*it);
								}
							}
							else
							{
								if (predicted_data[it->getSequence().toUnmodifiedString()] > 0)
								{
									hits_positive.push_back(*it);
								}
								else
								{
									hits_negative.push_back(*it);
								}
							}
						}
						temp_identification.setMetaValue("MZ",identifications[i].getMetaValue("MZ"));
						temp_identification.setMetaValue("RT",identifications[i].getMetaValue("RT"));
											
						temp_identification = identifications[i];
						temp_identification.setHits(hits_positive);
						identifications_positive.push_back(temp_identification);																														
						temp_identification.setHits(hits_negative);
						identifications_negative.push_back(temp_identification);
					}																														
				}
			}
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			if (separation_prediction)
			{
				IdXML_file.store(outputfile_name_positive,
															 protein_identifications,
															 identifications_positive);
				IdXML_file.store(outputfile_name_negative,
															 protein_identifications,
															 identifications_negative);
			}
			else
			{
				if (textfile_output)
				{
					writeStringLabelLines_(outputfile_name, predicted_data);
				}
				else
				{
					IdXML_file.store(outputfile_name,
																 protein_identifications,
																 identifications);
					writeDebug_("Linear correlation between predicted and measured rt is: "
											+ String(Math::pearsonCorrelationCoefficient(predicted_retention_times.begin(), 
																				predicted_retention_times.end(), 
																				performance_retention_times.begin(), 
																				performance_retention_times.end())), 1);														 
					writeDebug_("MSE between predicted and measured rt is: "
											+ String(Math::meanSquareError(predicted_retention_times.begin(), 
																				predicted_retention_times.end(), 
																				performance_retention_times.begin(), 
																				performance_retention_times.end())), 1);														 
					
				}
			}
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPRTPredict tool;
	return tool.main(argc,argv);
}
  
/// @endcond





