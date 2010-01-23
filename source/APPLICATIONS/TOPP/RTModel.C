// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_RTModel RTModel
	
	@brief Used to train a prediction model for peptide retention 
				 time prediction or peptide separation prediction.

	For retention time prediction, a support vector machine is 
	trained with peptide sequences and their measured retention 
	times.
	For peptide separation prediction two files have to be given.
	One file contains the positive examples (the peptides which
	are collected) and one file contains the negative peptides
	(the flowthrough peptides).
	
	This methods and applications of this model are described 
	in several publications:

	Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
	Statistical learning of peptide retention behavior in chromatographic separations: A new kernel-based approach for computational proteomics.
	BMC Bioinformatics 2007, 8:468 

	Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
	Improving Peptide Identification in Proteome Analysis by a Two-Dimensional Retention Time Filtering Approach
	J. Proteome Res. 2009, 8(8):4109-15


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
				SVMWrapper::OLIGO for our POBK (recommended))
		</li>
		<li>
			border_length: border length for the POBK
		</li>
		<li>
			k_mer_length: length of the signals considered in the 
			POBK
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
			nu: the nu parameter for nu-SVR
		</li>
		<li>
			p: the epsilon parameter for epsilon-SVR
		</li>
	</ul>
	
	<br>
	
	The last five parameters (sigma, degree, c, nu and p)
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
	Furthermore, you can specify the number of partitions for the CV with
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
	The model can be used in @ref TOPP_RTPredict, to predict retention times 
	for peptides or peptide separation depending on how you trained 
	the model.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_RTModel.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPRTModel
	: public TOPPBase
{
	public:
		TOPPRTModel()
			: TOPPBase("RTModel","Trains a model for the retention time prediction of peptides from a training set.")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","This is the name of the input file (RT prediction). It is assumed that the file type is IdXML. If it is just a textfile having a sequence and the corresponding rt per line, the 'textfile_input' flag has to be set.\n", false);
			setValidFormats_("in",StringList::create("idXML"));
			registerFlag_("textfile_input", "Has to be set if the input file is a textfile contatining a sequence with corresponding rt per line (separated by space).");
			registerInputFile_("in_positive","<file>","","input file with positive examples (peptide separation prediction)\n", false);
			setValidFormats_("in_positive",StringList::create("idXML"));
			registerInputFile_("in_negative","<file>","","input file with negative examples (peptide separation prediction)\n", false);
			setValidFormats_("in_negative",StringList::create("idXML"));
			registerOutputFile_("out","<file>","","output file: the model in libsvm format");
			registerStringOption_("svm_type","<type>","NU_SVR","the type of the svm (NU_SVR or EPSILON_SVR for RT prediction, automatically set\nto C_SVC for separation prediction)\n",false);
			setValidStrings_("svm_type",StringList::create("NU_SVR,NU_SVC,EPSILON_SVR,C_SVC"));
			registerDoubleOption_("nu","<float>",0.5,"the nu parameter [0..1] of the svm (for nu-SVR)",false);
			setMinFloat_("nu", 0);
			setMaxFloat_("nu", 1);
			registerDoubleOption_("p","<float>",0.1,"the epsilon parameter of the svm (for epsilon-SVR)",false);
			registerDoubleOption_("c","<float>",1,"the penalty parameter of the svm",false);
			registerStringOption_("kernel_type","<type>","OLIGO","the kernel type of the svm",false);
			setValidStrings_("kernel_type",StringList::create("LINEAR,RBF,POLY,OLIGO"));
			registerIntOption_("degree","<int>",1,"the degree parameter of the kernel function of the svm (POLY kernel)\n",false);
			setMinInt_("degree", 1);
			registerIntOption_("border_length","<int>",22,"length of the POBK",false);
			setMinInt_("border_length", 1);
			registerDoubleOption_("max_std","<float>",10,"max standard deviation for a peptide to be included (if there are several ones for one peptide string)(median is taken)",false);
			setMinFloat_("max_std", 0.);
			registerIntOption_("k_mer_length","<int>",1,"k_mer length of the POBK",false);
			setMinInt_("k_mer_length", 1);
			registerDoubleOption_("sigma","<float>",5,"sigma of the POBK",false);
			registerDoubleOption_("total_gradient_time","<time>",1,"the time (in seconds) of the gradient (only for RT prediction)", false);
			setMinFloat_("total_gradient_time", 0.00001);
			registerFlag_("first_dim_rt","if set the model will be built for first_dim_rt");
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
			registerDoubleOption_("p_start","<float>",1,"starting point of p",false);
			registerDoubleOption_("p_step_size","<float>",10,"step size point of p",false);
			registerDoubleOption_("p_stop","<float>",1000,"stopping point of p",false);
			registerDoubleOption_("c_start","<float>",1,"starting point of c",false);
			registerDoubleOption_("c_step_size","<float>",10,"step size of c",false);
			registerDoubleOption_("c_stop","<float>",1000,"stopping point of c",false);
			registerDoubleOption_("nu_start","<float>",0.3,"starting point of nu",false);
			setMinFloat_("nu_start", 0);
			setMaxFloat_("nu_start", 1);
			registerDoubleOption_("nu_step_size","<float>",1.2,"step size of nu",false);
			registerDoubleOption_("nu_stop","<float>",0.7,"stopping point of nu",false);
			setMinFloat_("nu_stop", 0);
			setMaxFloat_("nu_stop", 1);
			registerDoubleOption_("sigma_start","<float>",1,"starting point of sigma",false);
			registerDoubleOption_("sigma_step_size","<float>",1.3,"step size of sigma",false);
			registerDoubleOption_("sigma_stop","<float>",15,"stopping point of sigma",false);
			registerFlag_("skip_cv", "Has to be set if the cv should be skipped and the model should just be trained with the specified parameters.");
		}
		
	  void loadStringLabelLines_(String                		filename, 
														  std::vector<String>&  		sequences, 
														  std::vector<DoubleReal>&  labels)
	  {
	      TextFile text_file(filename.c_str(), true);
	      TextFile::iterator it;
	      std::vector<String> parts;
	      labels.clear();	
	      
	      it = text_file.begin();	
	      while(it != text_file.end())
	      {
				  it->split(' ', parts);
				  if (parts.size() == 2)
				  {
				      sequences.push_back(parts[0].trim());
				      labels.push_back(parts[1].trim().toDouble());
				      ++it;
				  }
				  else
				  {
				    it->split('\v', parts);	    
				    if (parts.size() == 2)
				    {
				      sequences.push_back(parts[0].trim());
				      labels.push_back(parts[1].trim().toDouble());
				      ++it;
				    }
				    else
				    {
				      it->split('\t', parts);	    
				      if (parts.size() == 2)
				      {
								sequences.push_back(parts[0].trim());
								labels.push_back(parts[1].trim().toDouble());
								++it;
				      }
				      else
				      {
				      	String debug_string = "found line '" + *it + "' in file which is not of the form <string> <label>\n";
								writeDebug_(debug_string, 1);
								++it;
				      }
				    }
				  }
	      }		
	  }
		

		ExitCodes main_(Int , const char**)
		{
			vector<ProteinIdentification> protein_identifications;
		  vector<PeptideIdentification> identifications;
			vector<ProteinIdentification> protein_identifications_negative;
		  vector<PeptideIdentification> identifications_negative;
		  vector<String> training_peptides;
		  vector<AASequence> training_modified_peptides;
		  vector< DoubleReal > training_retention_times;
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
			DoubleReal cv_quality = 0;
			map<SVMWrapper::SVM_parameter_type, DoubleReal> optimized_parameters;
			map<SVMWrapper::SVM_parameter_type, DoubleReal>::iterator parameters_iterator;
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
			map<String, DoubleReal> redundant_peptides;
			map<AASequence, DoubleReal> redundant_modified_peptides;
			DoubleReal max_std = 0.;
			bool textfile_input = false;
			SVMData training_sample;
			bool first_dim_rt = false;
			bool skip_cv = false;
			
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String inputfile_positives = getStringOption_("in_positive");
			String inputfile_negatives = "";
			String inputfile_name = "";
			if (inputfile_positives != "")
			{
				inputfile_negatives = getStringOption_("in_negative");
				if (inputfile_negatives != "")
				{
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
			}			
			String outputfile_name = getStringOption_("out");
			textfile_input = getFlag_("textfile_input");
			additive_cv = getFlag_("additive_cv");
			skip_cv = getFlag_("skip_cv");

			Real total_gradient_time = getDoubleOption_("total_gradient_time");
			max_std = getDoubleOption_("max_std");
			if (!separation_prediction && total_gradient_time	< 0)
			{
					writeLog_("No total gradient time given for RT prediction. Aborting!");
					printUsage_();
					return ILLEGAL_PARAMETERS;						
			}
 			//SVM type
 			String type = getStringOption_("svm_type");
			if (type == "NU_SVR" && !separation_prediction)
			{
				svm.setParameter(SVMWrapper::SVM_TYPE, NU_SVR);
			}
			else if (type == "EPSILON_SVR" && !separation_prediction)
			{
				svm.setParameter(SVMWrapper::SVM_TYPE, EPSILON_SVR);
			}
			else if ((separation_prediction && type == "C_SVC")
							 || separation_prediction)
			{
				svm.setParameter(SVMWrapper::SVM_TYPE, C_SVC);
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
 			if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVR || svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC)
 			{
				svm.setParameter(SVMWrapper::NU, getDoubleOption_("nu"));
			}
 			else if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == EPSILON_SVR)
 			{
				svm.setParameter(SVMWrapper::P, getDoubleOption_("p"));
			}
			
			//grid search parameters
			UInt degree_start = 0;
			UInt degree_step_size = 0;
			UInt degree_stop = 0;
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == POLY)
			{
				svm.setParameter(SVMWrapper::DEGREE, getIntOption_("degree"));

				if (setByUser_("degree_start") 
						&& setByUser_("degree_step_size") 
						&& setByUser_("degree_stop")
						&& !skip_cv)
				{
					degree_start = getIntOption_("degree_start");
					degree_step_size = getIntOption_("degree_step_size");
					if (!additive_cv && degree_step_size <= 1)
					{
						writeLog_("Step size of degree <= 1 and additive_cv is false. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					degree_stop = getIntOption_("degree_stop");

					start_values.insert(make_pair(SVMWrapper::DEGREE, degree_start));
					step_sizes.insert(make_pair(SVMWrapper::DEGREE, degree_step_size));
					end_values.insert(make_pair(SVMWrapper::DEGREE, degree_stop));	
				}
			}			
			DoubleReal p_start = 0.;
			DoubleReal p_step_size = 0.;
			DoubleReal p_stop = 0.;
			if (svm.getIntParameter(SVMWrapper::SVM_TYPE) == EPSILON_SVR)
			{							
				if (setByUser_("p_start") 
						&& setByUser_("p_step_size") 
						&& setByUser_("p_stop")
						&& !skip_cv)
				{
					p_start = getDoubleOption_("p_start");
					p_step_size = getDoubleOption_("p_step_size");
					if (!additive_cv && p_step_size <= 1)
					{
						writeLog_("Step size of p <= 1 and additive_cv is false. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					p_stop = getDoubleOption_("p_stop");

					start_values.insert(make_pair(SVMWrapper::P, p_start));
					step_sizes.insert(make_pair(SVMWrapper::P, p_step_size));
					end_values.insert(make_pair(SVMWrapper::P, p_stop));	
				}
			}
			DoubleReal c_start = 0.;
			DoubleReal c_step_size = 0.;
			DoubleReal c_stop = 0.;

			if (setByUser_("c_start") 
					&& setByUser_("c_step_size") 
					&& setByUser_("c_stop") 
					&& !skip_cv)
			{
				c_start = getDoubleOption_("c_start");
				c_step_size = getDoubleOption_("c_step_size");
				if (!additive_cv && c_step_size <= 1)
				{
					writeLog_("Step size of c <= 1 and additive_cv is false. Aborting!");
					return ILLEGAL_PARAMETERS;
				}
				c_stop = getDoubleOption_("c_stop");

				start_values.insert(make_pair(SVMWrapper::C, c_start));
				step_sizes.insert(make_pair(SVMWrapper::C, c_step_size));
				end_values.insert(make_pair(SVMWrapper::C, c_stop));	
			}			

			DoubleReal nu_start = 0.;
			DoubleReal nu_step_size = 0.;
			DoubleReal nu_stop = 0.;
			if (((svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVR 
					|| svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC))
					&& !skip_cv)
			{
				if (setByUser_("nu_start") && setByUser_("nu_step_size") && setByUser_("nu_stop"))
				{
					nu_start = getDoubleOption_("nu_start");
					nu_step_size = getDoubleOption_("nu_step_size");
					if (!additive_cv && nu_step_size <= 1)
					{
						writeLog_("Step size of nu <= 1 and additive_cv is false. Aborting!");
						return ILLEGAL_PARAMETERS;
					}
					nu_stop = getDoubleOption_("nu_stop");

					start_values.insert(make_pair(SVMWrapper::NU, nu_start));
					step_sizes.insert(make_pair(SVMWrapper::NU, nu_step_size));
					end_values.insert(make_pair(SVMWrapper::NU, nu_stop));	
				}			
			}
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO && setByUser_("border_length"))
			{
 				border_length = getIntOption_("border_length");
 			}
 			if (!setByUser_("border_length")
 					&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
 			{
				writeLog_("No border length given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}
			svm.setParameter(SVMWrapper::BORDER_LENGTH, border_length);
			if (setByUser_("sigma"))
			{
	 			sigma = getDoubleOption_("sigma");
	 		}
 			if ((!setByUser_("sigma") && !setByUser_("sigma_start"))
 					&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
 			{
				writeLog_("No sigma given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}
 			if (setByUser_("sigma"))
			{
				svm.setParameter(SVMWrapper::SIGMA, sigma);
			}
			else if (setByUser_("sigma_start"))
			{
				sigma = getDoubleOption_("sigma_start");
				svm.setParameter(SVMWrapper::SIGMA, sigma);
			}

			if (setByUser_("k_mer_length"))
			{
	 			k_mer_length = getIntOption_("k_mer_length");
			}
 			if (!setByUser_("k_mer_length")
 					&& svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
 			{
				writeLog_("No k-mer length given for POBK. Aborting!");
				return ILLEGAL_PARAMETERS;		
 			}

			sigma_start = 0.;
			sigma_step_size = 0.;
			sigma_stop = 0.;
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO 
					&& !skip_cv)
			{
				if (setByUser_("sigma_start") 
						&& setByUser_("sigma_step_size") 
						&& setByUser_("sigma_stop"))
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
			}
			if (start_values.size() > 0)
			{
 				number_of_runs = getIntOption_("number_of_runs");
				writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

 				number_of_partitions = getIntOption_("number_of_partitions");
				writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);				
			}
			
			first_dim_rt = getFlag_("first_dim_rt");
			
			Int debug_level = getIntOption_("debug");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			if (!separation_prediction)
			{
				if (textfile_input)
				{
					loadStringLabelLines_(inputfile_name, training_peptides, training_retention_times);
					for (Size i = 0; i < training_peptides.size(); ++i)
					{
						if (temp_type == SVMWrapper::OLIGO)
						{
							redundant_modified_peptides.insert(make_pair(AASequence(training_peptides[i]),
																									training_retention_times[i]));
						}
						else
						{
							redundant_peptides.insert(make_pair(training_peptides[i], training_retention_times[i]));
						}
					}
					training_peptides.clear();
					training_retention_times.clear();
				}
				else
				{
					String document_id;
					IdXMLFile().load(inputfile_name, protein_identifications, identifications, document_id);
				}
			}
			else
			{
				String document_id;
				IdXMLFile().load(inputfile_positives, protein_identifications, identifications, document_id);
				IdXMLFile().load(inputfile_negatives, protein_identifications_negative, identifications_negative, document_id);				
			}
		  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			if (!textfile_input)
			{
				for (Size i = 0; i < identifications.size(); i++)
				{
					Size temp_size = identifications[i].getHits().size();
					if (temp_size > 0)
					{
						if (temp_size == 1)
						{
							temp_peptide_hit = identifications[i].getHits()[0];
							if (separation_prediction)
							{
								training_retention_times.push_back(1.0);
								if (temp_type == SVMWrapper::OLIGO)
								{
									training_modified_peptides.push_back(temp_peptide_hit.getSequence());
								}
								else
								{
									training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
								}
							}	
							else
							{
								if (first_dim_rt)
								{
									if (temp_type != SVMWrapper::OLIGO)
									{
										redundant_peptides.insert(make_pair(temp_peptide_hit.getSequence().toUnmodifiedString(),
																												(DoubleReal)(identifications[i].getMetaValue("first_dim_rt"))));
									}
									else
									{
										redundant_modified_peptides.insert(make_pair(temp_peptide_hit.getSequence(),
																												(DoubleReal)(identifications[i].getMetaValue("first_dim_rt"))));
									}
								}
								else
								{
									if (temp_type != SVMWrapper::OLIGO)
									{
										redundant_peptides.insert(make_pair(temp_peptide_hit.getSequence().toUnmodifiedString(),
																												(DoubleReal)(identifications[i].getMetaValue("RT"))));
									}
									else
									{
										redundant_modified_peptides.insert(make_pair(temp_peptide_hit.getSequence(),
																												(DoubleReal)(identifications[i].getMetaValue("RT"))));
									}
								}
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
								writeLog_(String(it->getSequence().toUnmodifiedString()) + " score: " + String(it->getScore()));
							}
							return INPUT_FILE_CORRUPT;
						}
					}				
				}
			}
			
			// Getting a non redundant training set. If there are several copies of one peptide,
			// the standard deviation is calculated. If this std is less or equal to the
			// maximal allowed std 'max_std', the peptide is added to the training set
			// with the median as retention time. Unique peptides are added immediately.
			if (!separation_prediction && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
			{
				map<AASequence, DoubleReal>::iterator it = redundant_modified_peptides.begin();
				DoubleReal temp_variance = 0;
				DoubleReal temp_median = 0;
				DoubleReal temp_mean = 0;
				vector<DoubleReal> temp_values;
				pair< map<AASequence, DoubleReal>::iterator, map<AASequence, DoubleReal>::iterator > it_pair;
					
				while(it != redundant_modified_peptides.end())
			  {
			    temp_values.clear();
			    temp_variance = 0;
			 
			    it_pair = redundant_modified_peptides.equal_range(it->first);
			    for(it = it_pair.first; it != it_pair.second; ++it)
			    {
			      temp_values.push_back(it->second);
			    }
			    if (temp_values.size() == 1)
			    {
			      temp_median = temp_values[0];
			      temp_mean = temp_values[0];
			    }
			    else
			    {
			      sort(temp_values.begin(), temp_values.end());
			      if (temp_values.size() % 2 == 1)
			      {
							temp_median = temp_values[temp_values.size() / 2];
			      }
			      else
			      {
							temp_median = ((DoubleReal) temp_values[temp_values.size() / 2] 
															+ temp_values[temp_values.size() / 2 - 1]) / 2;
			      }
			
			      temp_mean = accumulate(temp_values.begin(), temp_values.end(), 0.) / temp_values.size();
			
			      for (Size j =0; j < temp_values.size(); ++j)
			      {
							temp_variance += (temp_values[j] - temp_mean) * (temp_values[j] - temp_mean);
			      }
			      temp_variance /= temp_values.size();
			    }
			    if (sqrt(temp_variance) <= max_std)
			    {
			    	training_modified_peptides.push_back(it_pair.first->first);
			    	training_retention_times.push_back(temp_median);
			    }
			    else
			    {
			    	debug_string = "Did not take peptide " + it_pair.first->first.toString() + " for training because"
			    		+ " there were several copies and std was " + String(temp_median) 
			    		+ " while " + String(max_std) + " was allowed.";
			    	writeDebug_(debug_string, 1);
			    }
			  }
			}
			
			if (!separation_prediction && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
			{
				map<String, DoubleReal>::iterator it = redundant_peptides.begin();
				DoubleReal temp_variance = 0;
				DoubleReal temp_median = 0;
				DoubleReal temp_mean = 0;
				vector<DoubleReal> temp_values;
				pair< map<String, DoubleReal>::iterator, map<String, DoubleReal>::iterator > it_pair;
					
				while(it != redundant_peptides.end())
			  {
			    temp_values.clear();
			    temp_variance = 0;
			 
			    it_pair = redundant_peptides.equal_range(it->first);
			    for(it = it_pair.first; it != it_pair.second; ++it)
			    {
			      temp_values.push_back(it->second);
			    }
			    if (temp_values.size() == 1)
			    {
			      temp_median = temp_values[0];
			      temp_mean = temp_values[0];
			    }
			    else
			    {
			      sort(temp_values.begin(), temp_values.end());
			      if (temp_values.size() % 2 == 1)
			      {
							temp_median = temp_values[temp_values.size() / 2];
			      }
			      else
			      {
							temp_median = ((DoubleReal) temp_values[temp_values.size() / 2] 
															+ temp_values[temp_values.size() / 2 - 1]) / 2;
			      }
			
			      temp_mean = accumulate(temp_values.begin(), temp_values.end(), 0.) / temp_values.size();
			
			      for (Size j =0; j < temp_values.size(); ++j)
			      {
							temp_variance += (temp_values[j] - temp_mean) * (temp_values[j] - temp_mean);
			      }
			      temp_variance /= temp_values.size();
			    }
			    if (sqrt(temp_variance) <= max_std)
			    {
			    	training_peptides.push_back(it_pair.first->first);
			    	training_retention_times.push_back(temp_median);
			    }
			    else
			    {
			    	debug_string = "Did not take peptide " + it_pair.first->first + " for training because"
			    		+ " there were several copies and std was " + String(temp_median) 
			    		+ " while " + String(max_std) + " was allowed.";
			    	writeDebug_(debug_string, 1);
			    }
			  }
			}
			
			// For separation prediction there are two files needed
			if (separation_prediction)
			{
				for (Size i = 0; i < identifications_negative.size(); i++)
				{
					Size temp_size = identifications_negative[i].getHits().size();
					if (temp_size > 0)
					{
						if (temp_size == 1)
						{
							temp_peptide_hit = identifications_negative[i].getHits()[0];
							if (temp_type == SVMWrapper::OLIGO)
							{
								training_modified_peptides.push_back(temp_peptide_hit.getSequence());
							}
							else
							{
								training_peptides.push_back(temp_peptide_hit.getSequence().toUnmodifiedString());
							}

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
								writeLog_(String(it->getSequence().toUnmodifiedString()) + " score: " + String(it->getScore()));
							}
							return INPUT_FILE_CORRUPT;
						}
					}				
				}
			}

			if (!separation_prediction)
			{
				for (Size i = 0; i < training_retention_times.size(); i++)
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
			else if (temp_type == SVMWrapper::OLIGO)
			{
				encoder.encodeProblemWithOligoBorderVectors(training_modified_peptides,
																										k_mer_length,
																										allowed_amino_acid_characters,
																										svm.getIntParameter(SVMWrapper::BORDER_LENGTH),
																										training_sample.sequences);
			}			
																													
			if (temp_type == SVMWrapper::OLIGO)
			{
				training_sample.labels = training_retention_times;
			}
			
			if (!skip_cv && start_values.size() > 0)
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
				if (temp_type == SVMWrapper::OLIGO)
				{
					debug_string = String(training_sample.sequences.size()) + " sequences for training, "
						+ training_sample.labels.size() + " labels for training";
					writeDebug_(debug_string, 1);
					
					cv_quality = svm.performCrossValidation(training_sample,
																								start_values,
																	 							step_sizes,
																	 							end_values,
																	 							number_of_partitions,
																	 							number_of_runs,
																	 							optimized_parameters,
																	 							additive_cv,
																	 							output_flag,
																	 							"performances_" + digest + ".txt");
				}
				else
				{				
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
				}												 												
				String debug_string = "Best parameters found in cross validation:";

				for(parameters_iterator = optimized_parameters.begin();
						parameters_iterator != optimized_parameters.end();
						parameters_iterator++)
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
					else if (parameters_iterator->first == SVMWrapper::P)
					{
						debug_string += " P: " + String(parameters_iterator->second);
					}
					else if (parameters_iterator->first == SVMWrapper::SIGMA)
					{
						debug_string += " sigma: " + String(parameters_iterator->second);
					}
				}
				debug_string += " with performance " + String(cv_quality);
				writeDebug_(debug_string, 1);
			}			

			if (temp_type == SVMWrapper::OLIGO)
			{
				svm.train(training_sample);
			}
			else
			{
				svm.train(encoded_training_sample);
			}
	
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			svm.saveModel(outputfile_name);

			// If the oligo-border kernel is used some additional information has to be stored
			if (temp_type == SVMWrapper::OLIGO)
			{
				training_sample.store(outputfile_name + "_samples");
				additional_parameters.setValue("kernel_type", temp_type);
				
				if (!separation_prediction)
				{
					svm.getSignificanceBorders(training_sample, sigmas);

					additional_parameters.setValue("sigma_0", sigmas.first); 
					additional_parameters.setValue("sigma_max", sigmas.second);
					if (first_dim_rt)
					{
						additional_parameters.setValue("first_dim_rt", "true");
					}
				}
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
	TOPPRTModel tool;
	return tool.main(argc,argv);
}

/// @endcond
