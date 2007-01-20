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

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page RTModel RTModel
	
	@brief Used to train a prediction model for peptide retention time prediction
	
	A support vector machine is trained with peptide sequences and their
	measured retention times. There are a number of parameters that
	can be changed for the svm (specified in the ini file):
	<ul>
		<li>
			svm_type: the type of the svm (can be NU_SVR or EPSILON_SVR)
		</li>
		<li>
			kernel_type: the kernel function (can be POLY for the 
				polynomial kernel or LINEAR for the linear kernel)
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
		<li>
			degree: the degree parameter for the polynomial kernel
		</li>
	</ul>
	
	<br>
	
	The last 4 parameters (c, nu, p and degree) can be used in a 
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
			Then the model is be stored.
		</li>
	</ol>
			
	<br>	
	The model can be used in RTPredict, to predict retention times 
	for peptides.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPRTModel
	: public TOPPBase
{
	public:
		TOPPRTModel()
			: TOPPBase("RTModel","Builds a model for retention time prediction of peptides from a training set")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerStringOption_("in","<file>","","input file in analysisXML format");
			registerStringOption_("out","<file>","","output file: the model in libsvm format");
			registerStringOption_("svm_type","<type>","NU_SVR","the type of the svm (NU_SVR or EPSILON_SVR)",false);
			registerDoubleOption_("nu","<float>",0.5,"the nu parameter [0..1] of the svm (for nu-SVR)",false);
			registerDoubleOption_("p","<float>",0.1,"the epsilon parameter of the svm (for epsilon-SVR)",false);
			registerDoubleOption_("c","<float>",1,"the penalty parameter of the svm",false);
			registerStringOption_("kernel_type","<type>","RBF","the kernel type of the svm (LINEAR, RBF, POLY or SIGMOID)",false);
			registerIntOption_("degree","<int>",1,"the degree parameter of the kernel function of the svm",false);
			registerDoubleOption_("total_gradient_time","<time>",0.0,"the time (in seconds) of the gradient");
			addEmptyLine_();
			addText_("Parameters for the grid search / cross validation:");
			registerIntOption_("number_of_runs","<n>",50,"number of runs for the CV",false);
			registerIntOption_("number_of_partitions","<n>",10,"number of CV partitions",false);
			registerIntOption_("degree_start","<int>",0,"starting point of degree",false);
			registerIntOption_("degree_step_size","<int>",0,"starting point of degree",false);
			registerIntOption_("degree_stop","<int>",0,"starting point of degree",false);
			registerDoubleOption_("p_start","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("p_step_size","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("p_stop","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("c_start","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("c_step_size","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("c_stop","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("nu_start","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("nu_step_size","<float>",0.0,"starting point of degree",false);
			registerDoubleOption_("nu_stop","<float>",0.0,"starting point of degree",false);
		}

		ExitCodes main_(int , char**)
		{
			vector<ProteinIdentification> protein_identifications;
		  vector<IdentificationData> identifications;
		  vector< String > training_peptides;
		  vector< DoubleReal > training_retention_times;
		  PeptideHit temp_peptide_hit;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			svm_problem* encoded_training_sample;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			map<SVM_parameter_type, DoubleReal> start_values;
			map<SVM_parameter_type, DoubleReal> step_sizes;
			map<SVM_parameter_type, DoubleReal> end_values;
			UnsignedInt number_of_partitions;
			UnsignedInt number_of_runs;
			DoubleReal cv_quality;
			map<SVM_parameter_type, DoubleReal>* optimized_parameters;
			map<SVM_parameter_type, DoubleReal>::iterator parameters_iterator;
			UnsignedInt maximum_sequence_length = 50;
	
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String inputfile_name = getStringOption_("in");
			inputFileReadable_(inputfile_name);
			String outputfile_name = getStringOption_("out");
			outputFileWritable_(outputfile_name);
			Real total_gradient_time = getDoubleOption_("total_gradient_time");		
 			//SVR type
 			String type = getStringOption_("svm_type");
			if (type == "NU_SVR")
			{
				svm.setParameter(SVM_TYPE, NU_SVR);
			}
			else if (type == "EPSILON_SVR")
			{
				svm.setParameter(SVM_TYPE, EPSILON_SVR);
			}
			else
			{
				writeLog_("Unknown svm type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;		
			}
			//Kernel type
 			type = getStringOption_("kernel_type");
			if (type == "POLY")
			{
				svm.setParameter(KERNEL_TYPE, POLY);
			}
			else if (type == "LINEAR")
			{
				svm.setParameter(KERNEL_TYPE, LINEAR);
			}			
			else if (type == "RBF")
			{
				svm.setParameter(KERNEL_TYPE, RBF);
			}
			else if (type == "SIGMOID")
			{
				svm.setParameter(KERNEL_TYPE, SIGMOID);
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
			UnsignedInt degree_start = getIntOption_("degree_start");
			UnsignedInt degree_step_size = getIntOption_("degree_step_size");
			UnsignedInt degree_stop = getIntOption_("degree_stop");
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

			number_of_runs = getIntOption_("number_of_runs");
			number_of_partitions = getIntOption_("number_of_partitions");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			AnalysisXMLFile().load(inputfile_name, protein_identifications, identifications);
		  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			for(UnsignedInt i = 0; i < identifications.size(); i++)
			{
				UnsignedInt temp_size = identifications[i].id.getPeptideHits().size();
				if (temp_size > 0)
				{
					if (temp_size == 1)
					{
						temp_peptide_hit = identifications[i].id.getPeptideHits()[0];
						training_peptides.push_back(temp_peptide_hit.getSequence());
						training_retention_times.push_back(identifications[i].rt);
					}
					else
					{
						writeLog_("For one spectrum there should not be more than one peptide."
								      "Please use the IDFilter with the -best_hits option to achieve this. Aborting!");
						writeLog_("Hits: ");
						for(vector<PeptideHit>::iterator it = identifications[i].id.getPeptideHits().begin(); 
								it != identifications[i].id.getPeptideHits().end(); 
								it++)
						{
							writeLog_(String(it->getSequence()) + " score: " + String(it->getScore()));
						}
						return INPUT_FILE_CORRUPT;
					}
				}				
			}

			for(UnsignedInt i = 0; i < training_retention_times.size(); i++)
			{
				training_retention_times[i] = training_retention_times[i] / total_gradient_time;
			}
			encoded_training_sample = 
				encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(training_peptides,
																																	&training_retention_times,
																																	allowed_amino_acid_characters,
																																	maximum_sequence_length);
																													
			if (start_values.size() > 0)
			{	
				optimized_parameters = svm.performCrossValidation(encoded_training_sample,
																												 	start_values,
																	 												step_sizes,
																	 												end_values,
																	 												&cv_quality,
																	 												number_of_partitions,
																	 												number_of_runs);
																	 												
				String debug_string = "Best parameters found in cross validation:";

				for(parameters_iterator = optimized_parameters->begin();
						parameters_iterator != optimized_parameters->end();
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
				}
				debug_string += " with performance " + String(cv_quality);
				writeDebug_(debug_string, 1);
			}			
			/// enabling probability estimates of the svm
			svm.setParameter(PROBABILITY, 1);
			
			svm.train(encoded_training_sample);
	
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			svm.saveModel(outputfile_name);
			
			return EXECUTION_OK;
		}
};


int main( int argc, char ** argv )
{
	TOPPRTModel tool;
	return tool.main(argc,argv);
}

/// @endcond
