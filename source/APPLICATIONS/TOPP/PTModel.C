// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_PTModel PTModel
	
	@brief Used to train a model for the prediction of proteotypic peptides.
	
	The input consists of two files: One file contains the positive examples (the peptides which
	are proteotypic) and the other contains the negative examples (the nonproteotypic peptides).
	
	Parts of this model has been described in the publication

	Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gr&ouml;pl, Oliver Kohlbacher and Knut Reinert
	LC-MSsim - a simulation software for Liquid Chromatography Mass Spectrometry data
	BMC Bioinformatics 2008, 9:423.

	There are a number of parameters which can be changed for the svm (specified in the ini file):
	<ul>
		<li>
			kernel_type: the kernel function (e.g., POLY for the 
				polynomial kernel, LINEAR for the linear kernel or RBF for the gaussian kernel); we recommend 
				SVMWrapper::OLIGO for our paired oligo-border kernel (POBK)
		</li>
		<li>
			border_length: border length for the POBK
		</li>
		<li>
			k_mer_length: length of the signals considered in the POBK
		</li>
		<li>
			sigma: the amount of positional smoothing for the POBK
		</li>
		<li>
			degree: the degree parameter for the polynomial kernel
		</li>
		<li>
			c: the penalty parameter of the svm
		</li>
		<li>
			nu: the nu parameter for nu-SVC
		</li>
	</ul>
	
	The last five parameters (sigma, degree, c, nu and p)
	are used in a cross validation (CV) to find the best parameters according to the 
	training set. Thus, you have to specify the start value of a
	parameter, the step size in which the parameters should be increased
	and a final value for the particular parameter such that the tested
	parameter is never bigger than the given final value. If you want
	to perform a cross validation, for example, for the parameter c, you
	have to specify <b>c_start</b>, <b>c_step_size</b> and <b>c_stop</b>
	in the ini file. Let's say you want to perform a CV for c from 0.1 to 2
	with step size 0.1. Open up your ini-file with INIFileEditor and modify the fields
	c_start, c_step_size, and c_stop accordingly.
	
	If the CV should test additional parameters in a certain range 
	you just include them analogously to the example above.
	Furthermore, you can specify the number of partitions for the CV with
	<b>number_of_partitions</b> in the ini file and the number of runs
	with <b>number_of_runs</b>.
	
	<br>
	Consequently you have two choices to use this application:
	
	<ol>
		<li> 
			Set the parameters of the svm: The PTModel application will train 
			the svm with the training data and store the svm model.
		</li>
		<li>
			Give a range of parameters for which a CV should be performed:
			The PTModel application will perform a CV to find the best 
			parameter combination in the given range and afterwards train
			the svm with the best parameters and the whole training data.
			Then the model is stored.
		</li>
	</ol>
			
	<br>	
	The model can be used in @ref TOPP_PTPredict, to predict the likelihood  
	for peptides to be proteotypic.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PTModel.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPTModel
	: public TOPPBase
{
	public:
		TOPPPTModel()
			: TOPPBase("PTModel","Trains a model for the prediction of proteotypic peptides from a training set.")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in_positive","<file>","","input file with positive examples\n");
			setValidFormats_("in_positive",StringList::create("idXML"));
			registerInputFile_("in_negative","<file>","","input file with negative examples\n");
			setValidFormats_("in_negative",StringList::create("idXML"));
			registerOutputFile_("out","<file>","","output file: the model in libsvm format");
			registerDoubleOption_("c","<float>",1,"the penalty parameter of the svm",false);
			registerStringOption_("svm_type","<type>","C_SVC","the type of the svm (NU_SVC or C_SVC)\n",false);
			setValidStrings_("svm_type",StringList::create("NU_SVC,C_SVC"));
			registerDoubleOption_("nu","<float>",0.5,"the nu parameter [0..1] of the svm (for nu-SVR)",false);
			setMinFloat_("nu", 0);
			setMaxFloat_("nu", 1);
			registerStringOption_("kernel_type","<type>","OLIGO","the kernel type of the svm",false);
			setValidStrings_("kernel_type",StringList::create("LINEAR,RBF,POLY,OLIGO"));
			registerIntOption_("degree","<int>",1,"the degree parameter of the kernel function of the svm (POLY kernel)\n",false);
			setMinInt_("degree", 1);
			registerIntOption_("border_length","<int>",22,"length of the POBK",false);
			setMinInt_("border_length", 1);
			registerIntOption_("k_mer_length","<int>",1,"k_mer length of the POBK",false);
			setMinInt_("k_mer_length", 1);
			registerDoubleOption_("sigma","<float>",5,"sigma of the POBK",false);
			registerIntOption_("max_positive_count","<int>",1000,"quantity of positive samples for training (randomly chosen if smaller than available quantity)",false);
			setMinInt_("max_positive_count", 1);
			registerIntOption_("max_negative_count","<int>",1000,"quantity of positive samples for training (randomly chosen if smaller than available quantity)",false);
			setMinInt_("max_negative_count", 1);
			registerFlag_("redundant","if the input sets are redundant and the redundant peptides should occur more than once in the training set, this flag has to be set");
			registerFlag_("additive_cv","if the step sizes should be interpreted additively (otherwise the actual value is multiplied\nwith the step size to get the new value");
			addEmptyLine_();
			addText_("Parameters for the grid search / cross validation:");
			registerIntOption_("number_of_runs","<int>",10,"number of runs for the CV",false);
			setMinInt_("number_of_runs", 1);
			registerIntOption_("number_of_partitions","<int>",10,"number of CV partitions",false);
			setMinInt_("number_of_partitions", 2);
			registerIntOption_("degree_start","<int>",1,"starting point of degree",false);
			setMinInt_("degree_start", 1);
			registerIntOption_("degree_step_size","<int>",2,"step size point of degree",false);
			registerIntOption_("degree_stop","<int>",4,"stopping point of degree",false);
			registerDoubleOption_("c_start","<float>",1,"starting point of c",false);
			registerDoubleOption_("c_step_size","<float>",100,"step size of c",false);
			registerDoubleOption_("c_stop","<float>",1000,"stopping point of c",false);
			registerDoubleOption_("nu_start","<float>",0.1,"starting point of nu",false);
			setMinFloat_("nu_start", 0);
			setMaxFloat_("nu_start", 1);
			registerDoubleOption_("nu_step_size","<float>",1.3,"step size of nu",false);
			registerDoubleOption_("nu_stop","<float>",0.9,"stopping point of nu",false);
			setMinFloat_("nu_stop", 0);
			setMaxFloat_("nu_stop", 1);
			registerDoubleOption_("sigma_start","<float>",1,"starting point of sigma",false);
			registerDoubleOption_("sigma_step_size","<float>",1.3,"step size of sigma",false);
			registerDoubleOption_("sigma_stop","<float>",15,"stopping point of sigma",false);
			registerFlag_("skip_cv", "Has to be set if the cv should be skipped and the model should just be trained with the specified parameters.");			
		}

		ExitCodes main_(Int , const char**)
		{
			vector<ProteinIdentification> protein_identifications;
		  vector<PeptideIdentification> identifications;
			vector<ProteinIdentification> protein_identifications_negative;
		  vector<PeptideIdentification> identifications_negative;
		  vector<String> training_peptides;
		  vector< DoubleReal > training_labels;
		  PeptideHit temp_peptide_hit;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			svm_problem* encoded_training_sample = 0;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			map<SVMWrapper::SVM_parameter_type, DoubleReal> start_values;
			map<SVMWrapper::SVM_parameter_type, DoubleReal> step_sizes;
			map<SVMWrapper::SVM_parameter_type, DoubleReal> end_values;
			DoubleReal sigma_start = 0;
			DoubleReal sigma_step_size = 0;
			DoubleReal sigma_stop = 0;
			UInt number_of_partitions = 0;
			UInt number_of_runs = 0;
			map<SVMWrapper::SVM_parameter_type, DoubleReal> optimized_parameters;
			map<SVMWrapper::SVM_parameter_type, DoubleReal>::iterator parameters_iterator;
			bool additive_cv = true;
			Param additional_parameters;
			Int temp_type = POLY;
			String debug_string = "";
			DoubleReal sigma = 0.1;
			UInt k_mer_length = 1;
			Int border_length = 0;
			bool non_redundant = false;
			bool skip_cv = getFlag_("skip_cv");
			
			svm.setParameter(SVMWrapper::PROBABILITY, 1);
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String inputfile_positives = getStringOption_("in_positive");
			String inputfile_negatives = getStringOption_("in_negative");;
			String temp_string = "";

			String outputfile_name = getStringOption_("out");

			UInt max_positive_count = getIntOption_("max_positive_count");
			UInt max_negative_count = getIntOption_("max_negative_count");

 			//SVM type
 			String type = getStringOption_("svm_type");
			if (type == "NU_SVC")
			{
				svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVC);
			}
			else if (type == "C_SVC")
			{
				svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
			}
			else
			{
				writeLog_("Illegal svm type given. Svm type has to be either "
									+ String("NU_SVC or C_SVC. Aborting!"));
				printUsage_();
				return ILLEGAL_PARAMETERS;		
			}
			//Kernel type
 			type = getStringOption_("kernel_type");
			if (type == "POLY")
			{
				svm.setParameter(SVMWrapper::KERNEL_TYPE, POLY);
				temp_type = POLY;
			}
			else if (type == "LINEAR")
			{
				svm.setParameter(SVMWrapper::KERNEL_TYPE, LINEAR);
				temp_type = LINEAR;
			}			
			else if (type == "RBF")
			{
				svm.setParameter(SVMWrapper::KERNEL_TYPE, RBF);
				temp_type = RBF;
			}
			else if (type == "OLIGO")
			{
				svm.setParameter(SVMWrapper::KERNEL_TYPE, SVMWrapper::OLIGO);
				temp_type = SVMWrapper::OLIGO;
			}
			else if (type == "SIGMOID")
			{
				svm.setParameter(SVMWrapper::KERNEL_TYPE, SIGMOID);
				temp_type = SIGMOID;
			}
			else
			{
				writeLog_("Unknown kernel type given. Aborting!");
				printUsage_();
				return ILLEGAL_PARAMETERS;		
			}
			
			//parameters		
			svm.setParameter(SVMWrapper::C, getDoubleOption_("c"));
			svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));
 			if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC)
 			{
				svm.setParameter(SVMWrapper::NU, getDoubleOption_("nu"));
			}

			//grid search parameters
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == POLY)
			{
				svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));
				if(!skip_cv)
					{
						DoubleReal degree_start = getIntOption_("degree_start");
						DoubleReal degree_step_size = getIntOption_("degree_step_size");
						if (!additive_cv && degree_step_size <= 1)
							{
								writeLog_("Step size of degree <= 1 and additive_cv is false. Aborting!");
								return ILLEGAL_PARAMETERS;
							}
						DoubleReal degree_stop = getIntOption_("degree_stop");
				
						start_values.insert(make_pair(SVMWrapper::DEGREE, degree_start));
						step_sizes.insert(make_pair(SVMWrapper::DEGREE, degree_step_size));
						end_values.insert(make_pair(SVMWrapper::DEGREE, degree_stop));
					}
			}			
			
      if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == C_SVC && !skip_cv)
			{			
        DoubleReal c_start = getDoubleOption_("c_start");
        DoubleReal c_step_size = getDoubleOption_("c_step_size");
        if (!additive_cv && c_step_size <= 1)
        {
          writeLog_("Step size of c <= 1 and additive_cv is false. Aborting!");
          return ILLEGAL_PARAMETERS;
        }
        DoubleReal c_stop = getDoubleOption_("c_stop");

        start_values.insert(make_pair(SVMWrapper::C, c_start));
        step_sizes.insert(make_pair(SVMWrapper::C, c_step_size));
        end_values.insert(make_pair(SVMWrapper::C, c_stop));
      }

      if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC && !skip_cv)
			{			
        DoubleReal nu_start = getDoubleOption_("nu_start");
        DoubleReal nu_step_size = getDoubleOption_("nu_step_size");
        if (!additive_cv && nu_step_size <= 1)
        {
          writeLog_("Step size of nu <= 1 and additive_cv is false. Aborting!");
          return ILLEGAL_PARAMETERS;
        }
        DoubleReal nu_stop = getDoubleOption_("nu_stop");

        start_values.insert(make_pair(SVMWrapper::NU, nu_start));
        step_sizes.insert(make_pair(SVMWrapper::NU, nu_step_size));
        end_values.insert(make_pair(SVMWrapper::NU, nu_stop));
			}			

			border_length = getIntOption_("border_length");
			svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);

			sigma = getDoubleOption_("sigma");
	 		svm.setParameter(SVMWrapper::SIGMA, sigma);

			k_mer_length = getIntOption_("k_mer_length");
			
			sigma_start = 0.;
			sigma_step_size = 0.;
			sigma_stop = 0.;
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO 
					&& !skip_cv)
			{
				sigma_start = getDoubleOption_("sigma_start");
				sigma_step_size = getDoubleOption_("sigma_step_size");
				if (!additive_cv && sigma_step_size <= 1)
					{
						writeLog_("Step size of sigma <= 1 and additive_cv is false. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
				sigma_stop = getDoubleOption_("sigma_stop");
				
				start_values.insert(make_pair(SVMWrapper::SIGMA, sigma_start));
				step_sizes.insert(make_pair(SVMWrapper::SIGMA, sigma_step_size));
				end_values.insert(make_pair(SVMWrapper::SIGMA, sigma_stop));
				
				debug_string = "CV from sigma = " + String(sigma_start) +
					" to sigma = " + String(sigma_stop) + " with step size " + 
					String(sigma_step_size);
				writeDebug_(debug_string, 1);			
			}

      if (!skip_cv && !start_values.empty())
			{
 				number_of_runs = getIntOption_("number_of_runs");
				writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

 				number_of_partitions = getIntOption_("number_of_partitions");
				writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);
				
				additive_cv = getFlag_("additive_cv");
			}
			
			Int debug_level = getIntOption_("debug");
			non_redundant = !(getFlag_("redundant"));
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			String document_id;
			IdXMLFile().load(inputfile_positives, protein_identifications, identifications, document_id);
			IdXMLFile().load(inputfile_negatives, protein_identifications_negative, identifications_negative, document_id);				
		  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			for (Size i = 0; i < identifications.size(); i++)
			{
				const vector<PeptideHit>& temp_peptide_hits = identifications[i].getHits();
				Size temp_size = temp_peptide_hits.size();
				if (temp_size > 0)
				{
					for (Size j = 0; j < temp_size; ++j)
					{
						temp_peptide_hit = temp_peptide_hits[j];
						temp_string = temp_peptide_hit.getSequence().toUnmodifiedString();
						if (!non_redundant 
								|| find(training_peptides.begin(), training_peptides.end(), temp_string) == training_peptides.end())
						{
							training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
						}
					}
				}
			}
			training_labels.resize(training_peptides.size(), 1.0);
			debug_string = String(training_labels.size()) + " positive sequences read";
			writeDebug_(debug_string, 1);							

			if (training_peptides.size() > max_positive_count)
			{
				random_shuffle(training_peptides.begin(), training_peptides.end());
				training_peptides.resize(max_positive_count, "");
				training_labels.resize(max_positive_count, 1.);
			}	
			debug_string = String(training_peptides.size()) + " positive sequences for training";				
			writeDebug_(debug_string, 1);							

			UInt counter = 0;
			
			vector<String> temp_training_peptides;
			for (Size i = 0; i < identifications_negative.size(); i++)
			{
				const vector<PeptideHit>& temp_peptide_hits = identifications_negative[i].getHits();
				Size temp_size = temp_peptide_hits.size();
				if (temp_size > 0)
				{
					for (Size j = 0; j < temp_size; ++j)
					{
						temp_peptide_hit = temp_peptide_hits[j];
						temp_string = temp_peptide_hit.getSequence().toUnmodifiedString();
						if (find(training_peptides.begin(), training_peptides.end(), temp_string) != training_peptides.end())
						{
							writeLog_("Peptides are not allowed to occur in the positive and the negative set. Example: '" + temp_string + "'");
							return ILLEGAL_PARAMETERS;
						}		

						if (!non_redundant 
								|| find(training_peptides.begin(), training_peptides.end(), temp_string) == training_peptides.end())
						{
							temp_training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
							training_labels.push_back(-1.0);
							++counter;
						}	
					}
				}
			}
			if (non_redundant)
			{
				debug_string = String(counter) + " non redundant negative sequences read";				
			}
			else
			{
				debug_string = String(counter) + " negative sequences read";				
			}
			writeDebug_(debug_string, 1);							
			if (temp_training_peptides.size() > max_negative_count)
			{
				random_shuffle(temp_training_peptides.begin(), temp_training_peptides.end());
				temp_training_peptides.resize(max_negative_count, "");
				training_labels.resize(training_peptides.size() + max_negative_count, -1.);
			}	
			training_peptides.insert(training_peptides.end(), 
															 temp_training_peptides.begin(), 
															 temp_training_peptides.end());
															 			
			debug_string = String(temp_training_peptides.size()) + " negative sequences for training";				
			writeDebug_(debug_string, 1);							
			temp_training_peptides.clear();

      if (temp_type == LINEAR || temp_type == POLY || temp_type == RBF)
      {
        UInt maximum_sequence_length = 50;
        encoded_training_sample =
            encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(training_peptides,
                                                                       training_labels,
                                                                       allowed_amino_acid_characters,
                                                                       maximum_sequence_length);
      }
      else if (temp_type == SVMWrapper::OLIGO)
      {
        encoded_training_sample =
            encoder.encodeLibSVMProblemWithOligoBorderVectors(training_peptides,
                                                              training_labels,
                                                              k_mer_length,
                                                              allowed_amino_acid_characters,
                                                              svm.getIntParameter(SVMWrapper::BORDER_LENGTH));
      }

      if ( !start_values.empty() )
			{	
				String digest = "";
				bool output_flag = false;
				if (debug_level >= 1)
				{
					output_flag = true;
					vector<String> parts;
					outputfile_name.split('/', parts);
					if (parts.empty())
					{
						digest = outputfile_name;
					}
					else
					{
						digest = parts[parts.size() - 1];
					}
				}				
        SVMData dummy;
        DoubleReal cv_quality = svm.performCrossValidation(encoded_training_sample,
                                                dummy,
                                                false,
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

        for ( parameters_iterator = optimized_parameters.begin();
						parameters_iterator != optimized_parameters.end();
            ++parameters_iterator )
				{
					svm.setParameter(parameters_iterator->first,
													 parameters_iterator->second);
					if (parameters_iterator->first == SVMWrapper::DEGREE)
					{
						debug_string += " degree: " + String(parameters_iterator->second);					
					}
					else if (parameters_iterator->first == SVMWrapper::C)
					{
						debug_string += " C: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == SVMWrapper::NU)
					{
						debug_string += " nu: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == SVMWrapper::SIGMA)
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
			if (temp_type == SVMWrapper::OLIGO)
			{
				encoder.storeLibSVMProblem(outputfile_name + "_samples", encoded_training_sample);
				additional_parameters.setValue("kernel_type", temp_type);
				
				if (temp_type == SVMWrapper::OLIGO)
				{
					additional_parameters.setValue("border_length", svm.getIntParameter(SVMWrapper::BORDER_LENGTH));
					additional_parameters.setValue("k_mer_length", k_mer_length);
					additional_parameters.setValue("sigma", svm.getDoubleParameter(SVMWrapper::SIGMA));
				}
				
				additional_parameters.store(outputfile_name + "_additional_parameters");
			}
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPPTModel tool;
	return tool.main(argc,argv);
}

/// @endcond
