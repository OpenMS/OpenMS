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
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_RTPredict RTPredict
	
	@brief This application is used to predict retention times 
				 for peptides or peptide separation.

  This methods and applications of this model are described 
  in several publications:

  Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
  Statistical learning of peptide retention behavior in chromatographic separations: A new kernel-based approach for computational proteomics.
  BMC Bioinformatics 2007, 8:468 

  Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
  Improving Peptide Identification in Proteome Analysis by a Two-Dimensional Retention Time Filtering Approach
  J. Proteome Res. 2009, 8(8):4109-15


	The input of this application 
	is an svm model and an IdXML
	file with peptide identifications. The svm model file is specified
	by the <b>svm_model</b> parameter in the command line or the ini file. 
	This file should have been produced by the @ref TOPP_RTModel application. 
	<br>
	For retention time prediction the peptide sequences are extracted 
	from the IdXML inputfile 
	and passed to the svm. The svm then predicts retention times
	according to the trained model. The predicted retention times
	are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" /> 
	@endcode inside the peptide entities in the IdXML output file.

	For separation prediction you have to specify two output file names.
	'out_id:positive' is the filename of the peptides which are predicted
	to be collected by the column and 'out_id:negative' is the file
	of the predicted flowthrough peptides.

  Retention time prediction and separation prediction cannot be combined!

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_RTPredict.cli

	@todo This need serious clean up! Combining certain input and output options will
          result in strange behaviour, especially when using text output/input.
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
      //input
			registerInputFile_("in_id","<file>","","peptides with precursor information", false);
			setValidFormats_("in_id",StringList::create("idXML"));
      registerInputFile_("in_text","<file>", "","peptides as text-based file", false);
      
      //output
      registerTOPPSubsection_("out_id", "Output files in idXML format");
      registerOutputFile_("out_id:file","<file>","","Output file with peptide RT prediction", false);
      setValidFormats_("out_id:file",StringList::create("idXML"));
      registerOutputFile_("out_id:positive","<file>","","Output file in IdXML format containing positive predictions (peptide separation prediction - requires negative file to be present as well)\n", false);
      setValidFormats_("out_id:positive",StringList::create("idXML"));
      registerOutputFile_("out_id:negative","<file>","","Output file in IdXML format containing negative predictions (peptide separation prediction - requires positive file to be present as well)\n", false);
      setValidFormats_("out_id:negative",StringList::create("idXML"));
      
      registerTOPPSubsection_("out_text","Output files in text format");
      registerOutputFile_("out_text:file","<file>","","Output file with predicted RT values", false);

			registerInputFile_("svm_model","<file>","","svm model in libsvm format (can be produced by RTModel)");
			registerDoubleOption_("total_gradient_time","<time>",1.0,"the time (in seconds) of the gradient (peptide RT prediction)", false);
			setMinFloat_("total_gradient_time", 0.00001);
			registerIntOption_("max_number_of_peptides", "<int>",100000,"the maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).\n", false, true);
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
					os << it->first << " " << it->second << "\n";
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
			vector<DoubleReal> all_predicted_retention_times;
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
			String svmfile_name = "";
			Real total_gradient_time = 1.f;
			bool separation_prediction = false;
			vector<PeptideIdentification> identifications_positive;
			vector<PeptideIdentification> identifications_negative;
			bool first_dim_rt = false;
			Size number_of_peptides = 0;
			Size max_number_of_peptides = getIntOption_("max_number_of_peptides");
			
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			
      String outputfile_name_positive = getStringOption_("out_id:positive");
      String outputfile_name_negative = getStringOption_("out_id:negative");
      // for separation prediction, we require both files to be present!
      if (outputfile_name_positive != "" || outputfile_name_negative != "")
			{
        if (outputfile_name_positive != "" && outputfile_name_negative != "")
				{
					separation_prediction = true;					
				}
				else
				{
					writeLog_("Both files for separation prediction required. Please specify the other one as well. Aborting!");
					return ILLEGAL_PARAMETERS;					
				}
			}

      // either or
      String input_id = getStringOption_("in_id");			
      String input_text = getStringOption_("in_text");			
      if (input_text != "" && input_id != "")
			{
        writeLog_("Two input parameter files given, only one allowed! Use either -in_id:file or -in_text:file!");
        return ILLEGAL_PARAMETERS;
			}
      else if (input_text == "" && input_id == "")
      {
        writeLog_("No input file given. Aborting...");
        return ILLEGAL_PARAMETERS;
      }

      // OUTPUT
      // (can use both)
      String output_id = getStringOption_("out_id:file");	
      String output_text = getStringOption_("out_text:file");			
      if (output_text == "" && output_id == "" && !separation_prediction)
      {
        writeLog_("No output files given. Aborting...");
        return ILLEGAL_PARAMETERS;
      }

			svmfile_name = getStringOption_("svm_model");	
			total_gradient_time = getDoubleOption_("total_gradient_time");

			

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			svm.loadModel(svmfile_name);
						
			if ((svm.getIntParameter(SVMWrapper::SVM_TYPE) == C_SVC || svm.getIntParameter(SVMWrapper::SVM_TYPE) == NU_SVC) && !separation_prediction)
			{
					writeLog_("You cannot perform peptide separation prediction with a model trained for"
										+ String("\npeptide retention time prediction. Aborting!"));
					return ILLEGAL_PARAMETERS;					
			}
			if ((svm.getIntParameter(SVMWrapper::SVM_TYPE) != C_SVC && svm.getIntParameter(SVMWrapper::SVM_TYPE) != NU_SVC) && separation_prediction)
			{
					writeLog_("You cannot perform peptide retention time prediction with a model trained for\n"
										+ String("peptide separation prediction. Aborting!"));
					return ILLEGAL_PARAMETERS;					
			}

			// Since the POBK is not included in the libsvm we have to load
			// additional parameters from additional files.
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
			{
				inputFileReadable_(svmfile_name + "_additional_parameters", "svm_model (derived)");
	
				Param additional_parameters;
				
				additional_parameters.load(svmfile_name + "_additional_parameters");
				if (additional_parameters.exists("first_dim_rt") 
						&& additional_parameters.getValue("first_dim_rt") != DataValue::EMPTY)
				{
					first_dim_rt = additional_parameters.getValue("first_dim_rt").toBool();
				}
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
			
			if ( input_text != "" )
			{
				loadStrings_(input_text, peptides);
				if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					for (Size i = 0; i < peptides.size(); ++i)
					{
						modified_peptides.push_back(AASequence(peptides[i]));
					}
					peptides.clear();
				}
			}
			else
			{
				String document_id;
				IdXML_file.load(input_id, protein_identifications, identifications, document_id);
			}
	  													
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			if (input_id != "")
			{
				for (Size i = 0; i < identifications.size(); i++)
				{
					temp_peptide_hits = identifications[i].getHits();
					for (Size j = 0; j < temp_peptide_hits.size(); j++)
					{
						if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
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
			if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
			{
				number_of_peptides = modified_peptides.size();
			}
			else
			{
				number_of_peptides = peptides.size();
			}			
										
			vector<DoubleReal> rts;
			rts.resize(number_of_peptides, 0);
			
			vector<String>::iterator it_from = peptides.begin();
			vector<String>::iterator it_to = peptides.begin();
			vector<AASequence>::iterator it_from_mod = modified_peptides.begin();
			vector<AASequence>::iterator it_to_mod = modified_peptides.begin();
			Size counter = 0;
			while(counter < number_of_peptides)
			{
				vector<String> temp_peptides;
				vector<AASequence> temp_modified_peptides;
				vector<DoubleReal> temp_rts;

				Size temp_counter = 0;
				if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
				{
					while(temp_counter <= max_number_of_peptides && it_to != peptides.end())
					{
						++it_to;
						++temp_counter;
					}
					temp_peptides.insert(temp_peptides.end(), it_from, it_to);
					//temp_peptides.insert(temp_peptides.end(), peptides.begin(), peptides.end());
					temp_rts.resize(temp_peptides.size(), 0);

					prediction_data = 
						encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(temp_peptides,
																																				temp_rts,
																															 					allowed_amino_acid_characters,
																															 					maximum_length);
					it_from = it_to;
				}
				else if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					while(temp_counter < max_number_of_peptides && it_to_mod != modified_peptides.end())
					{
						++it_to_mod;
						++temp_counter;
					}
					temp_modified_peptides.insert(temp_modified_peptides.end(), it_from_mod, it_to_mod);
					//					temp_modified_peptides.insert(temp_modified_peptides.end(), modified_peptides.begin(), modified_peptides.end());
					temp_rts.resize(temp_modified_peptides.size(), 0);

					encoder.encodeProblemWithOligoBorderVectors(temp_modified_peptides, 
																											k_mer_length, 
																											allowed_amino_acid_characters, 
																											border_length,
																											prediction_samples.sequences);
					prediction_samples.labels = temp_rts;
					it_from_mod = it_to_mod;
				}
				counter += temp_counter;

				if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
				{
					inputFileReadable_((svmfile_name + "_samples").c_str(), "svm_model (derived)");
	
					training_samples.load(svmfile_name + "_samples");
					svm.setTrainingSample(training_samples);				
	
					svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
					svm.setParameter(SVMWrapper::SIGMA, sigma);
					svm.predict(prediction_samples, predicted_retention_times);
					prediction_samples.labels.clear();
					prediction_samples.sequences.clear();
				}
				else
				{
					svm.predict(prediction_data, predicted_retention_times);
					LibSVMEncoder::destroyProblem(prediction_data);
				}
				for (Size i = 0; i < temp_counter; ++i)
				{
					if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO && output_text=="")
					{
						predicted_modified_data.insert(make_pair(temp_modified_peptides[i],
																						(predicted_retention_times[i] * total_gradient_time)));
					}
					else if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
					{
						predicted_data.insert(make_pair(temp_peptides[i],
																						(predicted_retention_times[i] * total_gradient_time)));
					}
					else
					{
						predicted_data.insert(make_pair(temp_modified_peptides[i].toString(),
																						(predicted_retention_times[i] * total_gradient_time)));
					}												
				}
				all_predicted_retention_times.insert(all_predicted_retention_times.end(), predicted_retention_times.begin(), predicted_retention_times.end());
				predicted_retention_times.clear();
			}						

			if (input_id != "")
			{
				if (!separation_prediction)
				{
					for (Size i = 0; i < identifications.size(); i++)
					{
						temp_peptide_hits = identifications[i].getHits();
						for (Size j = 0; j < temp_peptide_hits.size(); j++)
						{
							DoubleReal temp_rt = 0.;
							DoubleReal temp_p_value = 0.;
														
							if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
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
								temp_point.first = 0;
								if (identifications[i].metaValueExists("RT"))
								{
									temp_point.first = identifications[i].getMetaValue("RT");
								}
							}
							if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
							{
								temp_point.second = temp_rt;
								temp_p_value = svm.getPValue(sigma_0, sigma_max, temp_point);
							}
							if (first_dim_rt)
							{
								if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
								{
									temp_peptide_hits[j].setMetaValue("predicted_RT_p_value_first_dim",temp_p_value);
								}
								temp_peptide_hits[j].setMetaValue("predicted_RT_first_dim",temp_rt);
								performance_retention_times.push_back(identifications[i].getMetaValue("first_dim_rt"));					
							}
							else
							{
								if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
								{
									temp_peptide_hits[j].setMetaValue("predicted_RT_p_value",temp_p_value);
								}
								temp_peptide_hits[j].setMetaValue("predicted_RT",temp_rt);

								if (identifications[i].metaValueExists("RT"))
								{
									performance_retention_times.push_back(identifications[i].getMetaValue("RT"));
								}
								else
								{
									performance_retention_times.push_back(0);
								}
							}
						}
						identifications[i].setHits(temp_peptide_hits);				
					}
				}
				else // separation prediction
				{
					vector<PeptideHit> hits_positive;
					vector<PeptideHit> hits_negative;
					
					PeptideIdentification temp_identification;
	
					for (Size i = 0; i < identifications.size(); i++)
					{					
						hits_negative.clear();
						hits_positive.clear();
	
						temp_peptide_hits = identifications[i].getHits();
						for(vector<PeptideHit>::iterator it = temp_peptide_hits.begin();
								it != temp_peptide_hits.end();
								++it)
						{						
							if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
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

						if (identifications[i].metaValueExists("MZ"))
						{
							temp_identification.setMetaValue("MZ",identifications[i].getMetaValue("MZ"));
						}
						if (identifications[i].metaValueExists("RT"))
						{
							temp_identification.setMetaValue("RT",identifications[i].getMetaValue("RT"));
						}
											
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
				if (output_text != "") // text
				{
					writeStringLabelLines_(output_text, predicted_data);
				}
				if (output_id != "") // idXML
				{
					IdXML_file.store(output_id,
													 protein_identifications,
													 identifications);
					writeDebug_("Linear correlation between predicted and measured rt is: "
											+ String(Math::pearsonCorrelationCoefficient(all_predicted_retention_times.begin(), 
																				all_predicted_retention_times.end(), 
																				performance_retention_times.begin(), 
																				performance_retention_times.end())), 1);														 
					writeDebug_("MSE between predicted and measured rt is: "
											+ String(Math::meanSquareError(all_predicted_retention_times.begin(), 
																				all_predicted_retention_times.end(), 
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





