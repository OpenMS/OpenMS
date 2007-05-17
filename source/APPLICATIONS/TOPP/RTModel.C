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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page RTModel RTModel
	
	@brief Used to train a prediction model for peptide retention 
				 time prediction or peptide separation prediction.
	
	For retention time prediction a support vector machine is 
	trained with peptide sequences and their measured retention 
	times.
	For peptide separation prediction two files have to be given.
	One file contains the positive examples (the peptides which
	are collected) and one file contains the negative peptides
	(the flowthrough peptides).
	
	There are a number of parameters which
	can be changed for the svm (specified in the ini file):
	<ul>
		<li>
			svm_type: the type of the svm (can be NU_SVR or 
			EPSILON_SVR for RT prediction and is C_SVC for separation
			prediction)
		</li>
		<li>
			kernel_type: the kernel function (can be POLY for the 
				polynomial kernel or LINEAR for the linear kernel, or 
				OLIGO for our POBK (recommended))
		</li>
		<li>
			border_length: border length for the POBK
		</li>
		<li>
			sigma: the amount of positional smoothing for the POBK
		</li>
		<li>
			k_mer_length: length of the signals considered in the 
			POBK
		</li>
		<li>
			degree: the degree parameter for the polynomial kernel
		</li>
		<li>
			c: the penalty parameter of the svm
		</li>
		<li>
			nu: the nu parameter for nu-SVR
		</li>
		<li>
			p: the epsilon parameter for epsilon-SVR
		</li>
	</ul>
	
	<br>
	
	The last 4 parameters (sigma, k_mer_length, degree, c, nu and p)
	can be used in a 
	cross validation (CV) to find the best parameters according to the 
	training set. Therefore you have to specify the start value of a
	parameter, the step size in which the parameters should be increased
	and a final value for the particular parameter such that the tested
	parameter is never bigger than the given final value. If you want
	to perform a cross validation for example for the parameter c, you
	have to specify <b>c_start</b>, <b>c_step_size</b> and <b>c_stop</b>
	in the ini file. So if you want to perform a CV for c from 0.1 to 2
	with step size 0.1 you include the following lines into your ini-file:
	<ul>
		<li>
			@code <ITEM name="c_start" value="0.1" type="float"/> @endcode
		</li>
		<li>
			@code <ITEM name="c_step_size" value="0.1" type="float"/> @endcode
		</li>
		<li>
			@code <ITEM name="c_stop" value="2" type="float"/> @endcode
		</li>
	</ul>
	If the CV should test additional parameters in a certain range 
	you just include them analogously to the example above.
	Furthermore you can specify the number of partitions for the CV with
	<b>number_of_partitions</b> in the ini file and the number of runs
	with <b>number_of_runs</b>.
	
	<br>
	Consequently you have two choices to use this application:
	
	<ol>
		<li> 
			Set the parameters of the svm: The RTModel application will train 
			the svm with the training data and store the svm model
		</li>
		<li>
			Give a range of parameters for which a CV should be performed:
			The RTModel application will perform a CV to find the best 
			parameter combination in the given range and afterwards train
			the svm with the best parameters and the whole training data.
			Then the model is stored.
		</li>
	</ol>
			
	<br>	
	The model can be used in RTPredict, to predict retention times 
	for peptides or peptide separation depending on how you trained 
	the model.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPRTModel
	: public TOPPBase
{
	public:
		TOPPRTModel()
			: TOPPBase("RTModel","Builds a model for retention time prediction of peptides from a training set."+
								String("\nFurthermore the tool can be used to build a model for peptide separation prediction.")
								+ "\nIn this case one file with positive examples and one file with negative examples have to be given.")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>","","input file in IdXML format (RT prediction)", false);
			registerStringOption_("in_positive","<file>","","input file in IdXML format with positive examples (peptide separation prediction)", false);
			registerStringOption_("in_negative","<file>","","input file in IdXML format with negative examples (peptide separation prediction)", false);
			registerStringOption_("out","<file>","","output file: the model in libsvm format");
			registerStringOption_("svm_type","<type>","NU_SVR","the type of the svm (NU_SVR or EPSILON_SVR for RT prediction, automatically set to C_SVC for separation prediction)",false);
			registerDoubleOption_("nu","<float>",0.5,"the nu parameter [0..1] of the svm (for nu-SVR)",false);
			registerDoubleOption_("p","<float>",0.1,"the epsilon parameter of the svm (for epsilon-SVR)",false);
			registerDoubleOption_("c","<float>",1,"the penalty parameter of the svm",false);
			registerStringOption_("kernel_type","<type>","OLIGO","the kernel type of the svm (LINEAR, RBF, POLY, SIGMOID or OLIGO)",false);
			registerIntOption_("degree","<int>",1,"the degree parameter of the kernel function of the svm (POLY kernel)",false);
			registerIntOption_("border_length","<int>",0,"length of the POBK",false);
			registerIntOption_("k_mer_length","<int>",0,"k_mer length of the POBK",false);
			registerDoubleOption_("sigma","<float>",-1.0,"sigma of the POBK",false);
			registerDoubleOption_("total_gradient_time","<time>",-1.0,"the time (in seconds) of the gradient (only for RT prediction)", false);
			registerFlag_("additive_cv","if the step sizes should be interpreted additively (otherwise the actual value is multiplied with the step size to get the new value");
			addEmptyLine_();
			addText_("Parameters for the grid search / cross validation:");
			registerIntOption_("number_of_runs","<int>",50,"number of runs for the CV",false);
			registerIntOption_("number_of_partitions","<int>",10,"number of CV partitions",false);
			registerIntOption_("degree_start","<int>",0,"starting point of degree",false);
			registerIntOption_("degree_step_size","<int>",0,"step size point of degree",false);
			registerIntOption_("degree_stop","<int>",0,"stopping point of degree",false);
			registerDoubleOption_("p_start","<float>",0.0,"starting point of p",false);
			registerDoubleOption_("p_step_size","<float>",0.0,"step size point of p",false);
			registerDoubleOption_("p_stop","<float>",0.0,"stopping point of p",false);
			registerDoubleOption_("c_start","<float>",0.0,"starting point of c",false);
			registerDoubleOption_("c_step_size","<float>",0.0,"step size of c",false);
			registerDoubleOption_("c_stop","<float>",0.0,"stopping point of c",false);
			registerDoubleOption_("nu_start","<float>",0.0,"starting point of nu",false);
			registerDoubleOption_("nu_step_size","<float>",0.0,"step size of nu",false);
			registerDoubleOption_("nu_stop","<float>",0.0,"stopping point of nu",false);
			registerDoubleOption_("sigma_start","<float>",0.0,"starting point of sigma",false);
			registerDoubleOption_("sigma_step_size","<float>",0.0,"step size of sigma",false);
			registerDoubleOption_("sigma_stop","<float>",0.0,"stopping point of sigma",false);
		}

		ExitCodes main_(Int , char**)
		{
			vector<Identification> protein_identifications;
		  vector<PeptideIdentification> identifications;
			vector<Identification> protein_identifications_negative;
		  vector<PeptideIdentification> identifications_negative;
		  vector< String > training_peptides;
		  vector< DoubleReal > training_retention_times;
		  PeptideHit temp_peptide_hit;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			svm_problem* encoded_training_sample = 0;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			map<SVM_parameter_type, DoubleReal> start_values;
			map<SVM_parameter_type, DoubleReal> step_sizes;
			map<SVM_parameter_type, DoubleReal> end_values;
			DoubleReal sigma_start = 0;
			DoubleReal sigma_step_size = 0;
			DoubleReal sigma_stop = 0;
			UInt number_of_partitions = 0;
			UInt number_of_runs = 0;
			DoubleReal cv_quality = 0;
			map<SVM_parameter_type, DoubleReal> optimized_parameters;
			map<SVM_parameter_type, DoubleReal>::iterator parameters_iterator;
			UInt maximum_sequence_length = 50;
			bool additive_cv = true;
			Param additional_parameters;
			pair<DoubleReal, DoubleReal> sigmas;
			Int temp_type = POLY;
			String debug_string = "";
			DoubleReal sigma = 0.1;
			UInt k_mer_length = 1;
			Int border_length = 0;
			bool separation_prediction = false;
			
	
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String inputfile_positives = getStringOption_("in_positive");
			String inputfile_negatives = "";
			String inputfile_name = "";
			if (inputfile_positives != "")
			{
				inputFileReadable_(inputfile_positives);
				inputfile_negatives = getStringOption_("in_negative");
				if (inputfile_negatives != "")
				{
					inputFileReadable_(inputfile_negatives);
					separation_prediction = true;					
				}
				else
				{
					writeLog_("Positive peptides for separation prediction set but no negative peptides. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;		
				}
			}
			else
			{										
				inputfile_name = getStringOption_("in");
				inputFileReadable_(inputfile_name);
			}
			String outputfile_name = getStringOption_("out");
			outputFileWritable_(outputfile_name);
			Real total_gradient_time = getDoubleOption_("total_gradient_time");
			if (!separation_prediction && total_gradient_time	< 0)
			{
					writeLog_("No total gradient time given for RT prediction. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;						
			}
 			//SVR type
 			String type = getStringOption_("svm_type");
			if (type == "NU_SVR" && !separation_prediction)
			{
				svm.setParameter(SVM_TYPE, NU_SVR);
			}
			else if (type == "EPSILON_SVR" && !separation_prediction)
			{
				svm.setParameter(SVM_TYPE, EPSILON_SVR);
			}
			else if ((separation_prediction && type == "C_SVC")
							 || separation_prediction)
			{
				svm.setParameter(SVM_TYPE, C_SVC);
			}
			else
			{
				writeLog_("Illegal svm type given. Svm type has to be either "
									+ String("NU_SVR or EPSILON_SVR for rt prediction and ")
									+ "C_SVC for separation prediction. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;		
			}
			//Kernel type
 			type = getStringOption_("kernel_type");
			if (type == "POLY")
			{
				svm.setParameter(KERNEL_TYPE, POLY);
				temp_type = POLY;
			}
			else if (type == "LINEAR")
			{
				svm.setParameter(KERNEL_TYPE, LINEAR);
				temp_type = LINEAR;
			}			
			else if (type == "RBF")
			{
				svm.setParameter(KERNEL_TYPE, RBF);
				temp_type = RBF;
			}
			else if (type == "OLIGO")
			{
				svm.setParameter(KERNEL_TYPE, OLIGO);
				temp_type = OLIGO;
			}
			else if (type == "SIGMOID")
			{
				svm.setParameter(KERNEL_TYPE, SIGMOID);
				temp_type = SIGMOID;
			}
			else
			{
				writeLog_("Unknown kernel type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;		
			}
			
			//parameters		
			svm.setParameter(C, getDoubleOption_("c"));
			svm.setParameter(DEGREE, getIntOption_("degree"));
 			if (svm.getIntParameter(SVM_TYPE) == NU_SVR)
 			{
				svm.setParameter(NU, getDoubleOption_("nu"));
			}
 			else if (svm.getIntParameter(SVM_TYPE) == EPSILON_SVR)
 			{
				svm.setParameter(P, getDoubleOption_("p"));
			}
			
			//grid search parameters
			UInt degree_start = getIntOption_("degree_start");
			UInt degree_step_size = getIntOption_("degree_step_size");
			UInt degree_stop = getIntOption_("degree_stop");
			if (degree_start != 0 && degree_step_size != 0 && degree_stop != 0)
			{
				start_values.insert(make_pair(DEGREE, degree_start));
				step_sizes.insert(make_pair(DEGREE, degree_step_size));
				end_values.insert(make_pair(DEGREE, degree_stop));	
			}
			
			DoubleReal p_start = getDoubleOption_("p_start");
			DoubleReal p_step_size = getDoubleOption_("p_step_size");
			DoubleReal p_stop = getDoubleOption_("p_stop");
			if (p_start != 0.0  && p_step_size != 0.0  && p_stop != 0.0  && svm.getIntParameter(SVM_TYPE) == EPSILON_SVR)
			{
				start_values.insert(make_pair(P, p_start));
				step_sizes.insert(make_pair(P, p_step_size));
				end_values.insert(make_pair(P, p_stop));	
			}
			
			DoubleReal c_start = getDoubleOption_("c_start");
			DoubleReal c_step_size = getDoubleOption_("c_step_size");
			DoubleReal c_stop = getDoubleOption_("c_stop");
			if (c_start != 0.0 && c_step_size != 0.0 && c_stop != 0.0)
			{
				start_values.insert(make_pair(C, c_start));
				step_sizes.insert(make_pair(C, c_step_size));
				end_values.insert(make_pair(C, c_stop));	
			}			

			DoubleReal nu_start = getDoubleOption_("nu_start");
			DoubleReal nu_step_size = getDoubleOption_("nu_step_size");
			DoubleReal nu_stop = getDoubleOption_("nu_stop");
			if (nu_start != 0.0 && nu_step_size != 0.0 && nu_stop != 0.0 && svm.getIntParameter(SVM_TYPE) == NU_SVR)
			{
				start_values.insert(make_pair(NU, nu_start));
				step_sizes.insert(make_pair(NU, nu_step_size));
				end_values.insert(make_pair(NU, nu_stop));	
			}			

 			border_length = getIntOption_("border_length");
 			if (border_length == 0 
 					&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
 			{
				writeLog_("No border length given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}
			svm.setParameter(BORDER_LENGTH, border_length);
 			sigma = getDoubleOption_("sigma");
 			if (sigma < 0 
 					&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
 			{
				writeLog_("No sigma given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}
			svm.setParameter(SIGMA, sigma);
 			k_mer_length = getIntOption_("k_mer_length");
 			if (k_mer_length == 0 
 					&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
 			{
				writeLog_("No k-mer length given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}

			sigma_start = getDoubleOption_("sigma_start");
			sigma_step_size = getDoubleOption_("sigma_step_size");
			sigma_stop = getDoubleOption_("sigma_stop");

			if (sigma_step_size != 0
					&& svm.getIntParameter(KERNEL_TYPE) == OLIGO)
			{
				start_values.insert(make_pair(SIGMA, sigma_start));
				step_sizes.insert(make_pair(SIGMA, sigma_step_size));
				end_values.insert(make_pair(SIGMA, sigma_stop));
				
				debug_string = "CV from sigma = " + String(sigma_start) +
					 " to sigma = " + String(sigma_stop) + " with step size " + 
					 String(sigma_step_size);
				writeDebug_(debug_string, 1);			
			}			

			if (start_values.size() > 0)
			{
 				number_of_runs = getIntOption_("number_of_runs");
				writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

 				number_of_partitions = getIntOption_("number_of_partitions");
				writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);
				
				additive_cv = getFlag_("additive_cv");
			}
			
			Int debug_level = getIntOption_("debug");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			if (!separation_prediction)
			{
				IdXMLFile().load(inputfile_name, protein_identifications, identifications);
			}
			else
			{
				IdXMLFile().load(inputfile_positives, protein_identifications, identifications);
				IdXMLFile().load(inputfile_negatives, protein_identifications_negative, identifications_negative);				
			}
		  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			for(UInt i = 0; i < identifications.size(); i++)
			{
				UInt temp_size = identifications[i].getHits().size();
				if (temp_size > 0)
				{
					if (temp_size == 1)
					{
						temp_peptide_hit = identifications[i].getHits()[0];
						training_peptides.push_back(temp_peptide_hit.getSequence());
						if (separation_prediction)
						{
							training_retention_times.push_back(1.0);
						}	
						else
						{
							training_retention_times.push_back((UInt)(identifications[i].getMetaValue("RT")));
						}
					}
					else
					{
						writeLog_("For one spectrum there should not be more than one peptide."
								      "Please use the IDFilter with the -best_hits option to achieve this. Aborting!");
						writeLog_("Hits: ");
						for(vector<PeptideHit>::const_iterator it = identifications[i].getHits().begin(); 
								it != identifications[i].getHits().end(); 
								it++)
						{
							writeLog_(String(it->getSequence()) + " score: " + String(it->getScore()));
						}
						return INPUT_FILE_CORRUPT;
					}
				}				
			}
			// For separation prediction there are two files needed
			if (separation_prediction)
			{
				for(UInt i = 0; i < identifications_negative.size(); i++)
				{
					UInt temp_size = identifications_negative[i].getHits().size();
					if (temp_size > 0)
					{
						if (temp_size == 1)
						{
							temp_peptide_hit = identifications_negative[i].getHits()[0];
							training_peptides.push_back(temp_peptide_hit.getSequence());

							training_retention_times.push_back(-1.0);
						}
						else
						{
							writeLog_("For one spectrum there should not be more than one peptide."
									      "Please use the IDFilter with the -best_hits option to achieve this. Aborting!");
							writeLog_("Hits: ");
							for(vector<PeptideHit>::const_iterator it = identifications_negative[i].getHits().begin(); 
									it != identifications_negative[i].getHits().end(); 
									it++)
							{
								writeLog_(String(it->getSequence()) + " score: " + String(it->getScore()));
							}
							return INPUT_FILE_CORRUPT;
						}
					}				
				}
			}

			if (!separation_prediction)
			{
				for(UInt i = 0; i < training_retention_times.size(); i++)
				{
					training_retention_times[i] = training_retention_times[i] / total_gradient_time;
				}
			}
			if (temp_type == LINEAR || temp_type == POLY || temp_type == RBF)
			{
				encoded_training_sample = 
					encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(training_peptides,
																																	training_retention_times,
																																	allowed_amino_acid_characters,
																																	maximum_sequence_length);
			}
			else if (temp_type == OLIGO)
			{
				encoded_training_sample = 
					encoder.encodeLibSVMProblemWithOligoBorderVectors(training_peptides,
																														training_retention_times,
																														k_mer_length,
																														allowed_amino_acid_characters,
																														svm.getIntParameter(BORDER_LENGTH));
			}			
																													
			if (start_values.size() > 0)
			{	
				String digest = "";
				bool output_flag = false;
				if (debug_level >= 1)
				{
					output_flag = true;
					vector<String> parts;
					inputfile_name.split('/', parts);
					if (parts.size() == 0)
					{
						digest = inputfile_name;
					}
					else
					{
						digest = parts[parts.size() - 1];
					}
				}				
				cv_quality = svm.performCrossValidation(encoded_training_sample,
																								start_values,
																	 							step_sizes,
																	 							end_values,
																	 							number_of_partitions,
																	 							number_of_runs,
																	 							optimized_parameters,
																	 							additive_cv,
																	 							output_flag,
																	 							"performances_" + digest + ".txt");
																	 												
				String debug_string = "Best parameters found in cross validation:";

				for(parameters_iterator = optimized_parameters.begin();
						parameters_iterator != optimized_parameters.end();
						parameters_iterator++)
				{
					svm.setParameter(parameters_iterator->first,
													 parameters_iterator->second);
					if (parameters_iterator->first == DEGREE)
					{
						debug_string += " degree: " + String(parameters_iterator->second);					
					}
					else if (parameters_iterator->first == C)
					{
						debug_string += " C: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == NU)
					{
						debug_string += " nu: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == P)
					{
						debug_string += " P: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == SIGMA)
					{
						debug_string += " sigma: " + String(parameters_iterator->second);
					}
				}
				debug_string += " with performance " + String(cv_quality);
				writeDebug_(debug_string, 1);
			}			

			svm.train(encoded_training_sample);
	
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			svm.saveModel(outputfile_name);

			// If the oligo-border kernel is used some additional information has to be stored
			if (temp_type == OLIGO)
			{
				encoder.storeLibSVMProblem(outputfile_name + "_samples", encoded_training_sample);
				additional_parameters.setValue((string)"kernel_type", (Int) temp_type);
				
				if (!separation_prediction)
				{	
					svm.getSignificanceBorders(encoded_training_sample, sigmas);
					additional_parameters.setValue((string)"sigma_0", sigmas.first); 
					additional_parameters.setValue((string)"sigma_max", sigmas.second);
				}
				if (temp_type == OLIGO)
				{
					additional_parameters.setValue((string)"border_length", svm.getIntParameter(BORDER_LENGTH));
					additional_parameters.setValue((string)"k_mer_length", (Int) k_mer_length);
					additional_parameters.setValue((string)"sigma", svm.getDoubleParameter(SIGMA));
				}
				
				additional_parameters.store(outputfile_name + "_additional_parameters");
			}
			
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPRTModel tool;
	return tool.main(argc,argv);
}

/// @endcond
