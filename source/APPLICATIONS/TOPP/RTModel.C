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
			parameter combination in the given range and afterwars train
			the svm with the best parameters and the whole training data.
			Afterwards the model will be stored.
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
			: TOPPBase("RTModel")
		{
			
		}
	
	protected:
		void printToolUsage_() const
		{
			cerr << endl
		       << getToolName() << " -- Builds a model for retention time" 
		       << " prediction of peptides. Peptides with the associated"
		       << " retention times are used to train the model."
		       << "Version: " << VersionInfo::getVersion() << endl
		       << endl
		       << "Usage:" << endl
					 << " " << getToolName() << " [options]" << endl
					 << endl
					 << "Options are:" << endl
					 << "  -in <file>              input file in analysisXML format (default read from INI file)" << endl
					 << "  -out <file>             output file: the model in libsvm format (default read from INI file)" << endl
					 << "  -total_gradient_time    the time (in seconds) of the gradient (default read from INI file)" << endl
					 << "  -c                      the penalty parameter of the svm (default read from INI file)" << endl
					 << "  -nu                     the nu parameter of the svm (for nu-SVR) (default read from INI file)" << endl
					 << "  -degree                 the degree parameter of the kernel function of the svm (default read from INI file)" << endl
					 << "  -p                      the epsilon parameter of the svm (for epsilon-SVR) (default read from INI file)" << endl
					 << "  -kernel_type            the kernel type of the svm (LINEAR, RBF, POLY or SIGMOID) (default read from INI file)" << endl
					 << "  -svm_type               the type of the svm (nu-SVR or epsilon-SVR) (default read from INI file)" << endl
					 << endl ;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-out"] = "out";
			options_["-in"] = "in";
			options_["-total_gradient_time"] = "total_gradient_time";
			options_["-c"] = "c";
			options_["-nu"] = "nu";
			options_["-degree"] = "degree";
			options_["-p"] = "p";
			options_["-kernel_type"] = "kernel_type";
			options_["-svm_type"] = "svm_type";
			options_["--help"] = "help";
		}
	
		void printToolHelpOpt_() const
		{
			cerr << endl
		       << getToolName() << endl
		       << endl
		       << "INI options:" << endl
					 << "  in                        input file" << endl
					 << "  out                       output file" << endl
					 << "  total_gradient_time       the time (in seconds) of the gradient" << endl
					 << "  c                         the penalty parameter of the svm" << endl
					 << "  nu                        the nu parameter of the svm (for nu-SVR)" << endl
					 << "  degree                    the degree parameter of the kernel function of the svm" << endl
					 << "  p                         the epsilon parameter of the svm (for epsilon-SVR)" << endl
					 << "  kernel_type               the kernel type of the svm (LINEAR, RBF, POLY or SIGMOID)" << endl
					 << "  svm_type                  the type of the svm (nu-SVR or epsilon-SVR)" << endl
					 << endl << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"in\" value=\"input.analysisXML\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"svm.model\" type=\"string\"/>" << endl
					 << "  <!-- The penalty parameter for generalisation. -->" << endl
  		     << "  <ITEM name=\"c\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"c_start\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"c_step_size\" value=\"0.3\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"c_stop\" value=\"2\" type=\"float\"/>" << endl
  		     << "  <!-- The nu parameter in NU_SVR. -->" << endl
  		     << "  <ITEM name=\"nu\" value=\"0.5\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"nu_start\" value=\"0.4\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"nu_step_size\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"nu_stop\" value=\"0.6\" type=\"float\"/>" << endl
  		     << "  <!-- The degree of the polynomial kernel. -->" << endl
  		     << "  <ITEM name=\"degree\" value=\"1\" type=\"int\"/>" << endl
  		     << "  <ITEM name=\"degree_start\" value=\"1\" type=\"int\"/>" << endl
  		     << "  <ITEM name=\"degree_step_size\" value=\"1\" type=\"int\"/>" << endl
  		     << "  <ITEM name=\"degree_stop\" value=\"3\" type=\"int\"/>" << endl
  		     << "  <!-- The epsilon parameter in EPSILON_SVR (not used in this example)-->" << endl
  		     << "  <ITEM name=\"p\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"p_start\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"p_step_size\" value=\"0.1\" type=\"float\"/>" << endl
  		     << "  <ITEM name=\"p_stop\" value=\"0.2\" type=\"float\"/>" << endl;
 
 					 
		}

		ExitCodes main_(int , char**)
		{
			// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
			String ini_location;
			// path to the log file
			String logfile = "";
			// log filestream (as long as the real logfile is not setermined yet)
			ofstream log;
			String inputfile_name;
			String outputfile_name;
		  vector<ProteinIdentification> protein_identifications;
		  vector<IdentificationData> identifications;
		  vector<DoubleReal> training_retention_times_double;
		  vector< String > training_peptides;
		  vector< DoubleReal > training_retention_times;
		  UnsignedInt temp_size = 0;
		  PeptideHit temp_peptide_hit;
			SVMWrapper svm;
			LibSVMEncoder encoder;
			svm_problem* encoded_training_sample;
			String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
			String type = "";
			String parameter = "";
			Real total_gradient_time = 1.f;
			map<SVM_parameter_type, DoubleReal> start_values;
			map<SVM_parameter_type, DoubleReal> step_sizes;
			map<SVM_parameter_type, DoubleReal> end_values;
			UnsignedInt number_of_partitions = 5;
			UnsignedInt number_of_runs = 20;
			DoubleReal cv_quality;
			UnsignedInt degree_start = 0;
			UnsignedInt degree_step_size = 0;
			UnsignedInt degree_stop = 0;
			DoubleReal c_start = 0;
			DoubleReal c_step_size = 0;
			DoubleReal c_stop = 0;
			DoubleReal nu_start = 0;
			DoubleReal nu_step_size = 0;
			DoubleReal nu_stop = 0;
			DoubleReal p_start = 0;
			DoubleReal p_step_size = 0;
			DoubleReal p_stop = 0;
			map<SVM_parameter_type, DoubleReal>* optimized_parameters;
			map<SVM_parameter_type, DoubleReal>::iterator parameters_iterator;
			String start;
			String stop;
			String step_size;
			String debug_string;
			UnsignedInt maximum_sequence_length = 50;
	
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

 			type = getParamAsString_("svm_type", "NU_SVR");
			writeDebug_(String("Svm type: ") + type, 1);
			
			if (type == "NU_SVR")
			{
				svm.setParameter(SVM_TYPE, NU_SVR);
			}
			else if (type == "EPSILON_SVR")
			{
				svm.setParameter(SVM_TYPE, EPSILON_SVR);
			}			

 			type = getParamAsString_("kernel_type", "POLY");
			writeDebug_(String("Kernel type: ") + type, 1);
			if (type == "POLY")
			{
				svm.setParameter(KERNEL_TYPE, POLY);
			}
			else if (type == "LINEAR")
			{
				svm.setParameter(KERNEL_TYPE, LINEAR);
			}			

 			parameter = getParamAsString_("c", "1");
			writeDebug_(String("c: ") + parameter, 1);
			svm.setParameter(C, parameter.toDouble());

 			parameter = getParamAsString_("nu", "0.5");
 			if (svm.getIntParameter(SVM_TYPE) == NU_SVR)
 			{
				writeDebug_(String("nu: ") + parameter, 1);
				svm.setParameter(NU, parameter.toDouble());
			}

 			parameter = getParamAsString_("degree", "1");
			writeDebug_(String("degree: ") + parameter, 1);
			svm.setParameter(DEGREE, parameter.toInt());

 			parameter = getParamAsString_("p", "0.1");
 			if (svm.getIntParameter(SVM_TYPE) == EPSILON_SVR)
 			{
				writeDebug_(String("p (epsilon in Epsilon SVR): ") + parameter, 1);
				svm.setParameter(P, parameter.toDouble());
			}

			start = getParamAsString_("degree_start");
			step_size = getParamAsString_("degree_step_size");
			stop = getParamAsString_("degree_stop");

			if (start != "" && step_size != "" && stop != "")
			{
				degree_start = start.toInt();
				degree_step_size = step_size.toInt();
				degree_stop = stop.toInt();
				start_values.insert(make_pair(DEGREE, degree_start));
				step_sizes.insert(make_pair(DEGREE, degree_step_size));
				end_values.insert(make_pair(DEGREE, degree_stop));
				
				debug_string = "CV from degree = " + String(degree_start) +
					 " to degree = " + String(degree_stop) + " with step size " + 
					 String(degree_step_size);
				writeDebug_(debug_string, 1);			
			}			

			start = getParamAsString_("p_start");
			step_size = getParamAsString_("p_step_size");
			stop = getParamAsString_("p_stop");

			if (start != "" 
					&& step_size != "" 
					&& stop != "" 
					&& svm.getIntParameter(SVM_TYPE) == EPSILON_SVR)
			{
				p_start = start.toFloat();
				p_step_size = step_size.toFloat();
				p_stop = stop.toFloat();
				start_values.insert(make_pair(P, p_start));
				step_sizes.insert(make_pair(P, p_step_size));
				end_values.insert(make_pair(P, p_stop));
				
				debug_string = "CV from p = " + String(p_start) +
					 " to p = " + String(p_stop) + " with step size " + 
					 String(p_step_size);
				writeDebug_(debug_string, 1);			
			}			

			start = getParamAsString_("c_start");
			step_size = getParamAsString_("c_step_size");
			stop = getParamAsString_("c_stop");

			if (start != "" && step_size != "" && stop != "")
			{
				c_start = start.toFloat();
				c_step_size = step_size.toFloat();
				c_stop = stop.toFloat();
				start_values.insert(make_pair(C, c_start));
				step_sizes.insert(make_pair(C, c_step_size));
				end_values.insert(make_pair(C, c_stop));
				
				debug_string = "CV from c = " + String(c_start) +
					 " to c = " + String(c_stop) + " with step size " + 
					 String(c_step_size);
				writeDebug_(debug_string, 1);			
			}			

			start = getParamAsString_("nu_start");
			step_size = getParamAsString_("nu_step_size");
			stop = getParamAsString_("nu_stop");

			if (start != "" 
					&& step_size != "" 
					&& stop != "" 
					&& svm.getIntParameter(SVM_TYPE) == NU_SVR)
			{
				nu_start = start.toFloat();
				nu_step_size = step_size.toFloat();
				nu_stop = stop.toFloat();
				start_values.insert(make_pair(NU, nu_start));
				step_sizes.insert(make_pair(NU, nu_step_size));
				end_values.insert(make_pair(NU, nu_stop));
				
				debug_string = "CV from nu = " + String(nu_start) +
					 " to nu = " + String(nu_stop) + " with step size " + 
					 String(nu_step_size);
				writeDebug_(debug_string, 1);			
			}			

			if (start_values.size() > 0)
			{
 				number_of_runs = getParamAsString_("number_of_runs", "50").toInt();
				writeDebug_(String("Number of CV runs: ") + String(number_of_runs), 1);

 				number_of_partitions = getParamAsString_("number_of_partitions", "10").toInt();
				writeDebug_(String("Number of CV partitions: ") + String(number_of_partitions), 1);
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
			
			if (!File::writable(outputfile_name))
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, outputfile_name);
			}
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			
			
			AnalysisXMLFile().load(inputfile_name,
														protein_identifications,
														identifications);
		  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			for(UnsignedInt i = 0; i < identifications.size(); i++)
			{
				if ((temp_size = identifications[i].id.getPeptideHits().size()) > 0)
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
								      "Please use the IDFilter with the -strict option to achieve this. Aborting!");
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
				encoder.encodeLIBSVMProblemWithCompositionAndLengthVectors(training_peptides,
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
																	 												
				debug_string = "Best parameters found in cross validation:";

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
