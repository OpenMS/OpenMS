#// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Identification.h>

#include "TOPPBase.h"

#include <qfileinfo.h>
#include <qfile.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IDFilter IDFilter
	
	@brief Filters identifications depending on certain criteria.
	
	This component is used to filter the identifications found by
	an identification tool like Mascot&copy;. The identifications 
	can be filtered by different criteria.
	
	<ul>
		<li> 
			peptide significance threshold fraction: This parameter 
			specifies which amount of the significance threshold should 
			be reached by a peptide to be kept. If for example a peptide
			has score 30 and the significance threshold is 40, the 
			peptide will only be kept by the filter if the significance 
			threshold fraction is set to 0.75 or lower. The value for this
			parameter can be set in the ini file with the 
			<b>peptide_significance_threshold_fraction</b> parameter.
		</li> 
		<li> 
			protein significance threshold fraction: This parameter 
			behaves in the same way as the peptide significance threshold 
			fraction parameter. The only difference is that it is used
			to filter protein hits. The value for this
			parameter can be set in the ini file with the 
			<b>protein_significance_threshold_fraction</b> parameter.
		</li>
		<li>
			sequences file (in FASTA format): If you know which proteins
			are in the measured sample you can specify a FASTA file 
			which contains the protein sequences of the sample. All 
			peptides which are not a substring of a protein contained
			in the sequences file will be filtered out. The name of the
			sequences file can be specified by the <b>sequences_file</b>
			parameter in the ini file.
		</li>
		<li>
			retention time: To filter identifications according to their 
			predicted retention times you have to set two additional parameters
			in the ini file <b>total_gradient_time</b> which is the total 
			number of seconds that the gradient took and <b>allowed_deviation</b>
			which defines a factor for the amount of deviation from the original
			retention time. To use this filter mode you have to supply an
			analysisXML file that is produced by the RTPredict component.
		</li>
	</ul>
	
	<br>
		The first three filtering possibilities have in common that for 
		one spectrum only the hits with the maximal score are kept. 
		If there is more than one peptide hit with maximal score for 
		one spectrum you can specify by the <b>strict</b> option in the 
		commandline or the ini file that you want to drop all of them. 
		If you do not specify this option these identifications will be kept.						 

	@ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFilter
	: public TOPPBase
{
 public:
	TOPPIDFilter()
		: TOPPBase("IDFilter")
	{
			
	}
	
 protected:
	void printToolUsage_()
	{
		cerr << endl
				 << tool_name_ << " -- annotates MS/MS spectra using Mascot" << endl
				 << endl
				 << "Usage:" << endl
				 << " " << tool_name_ << " [options]" << endl
				 << endl
				 << "Options are:" << endl
				 << "  -in <file>   		          input file in mzData " 
				 << "(default read from INI file)" << endl
				 << "  -out <file>  		          output file in analysisXML/Mascot generic format "
				 << "(default read from INI file)" << endl
				 << "  -strict      		          flag indicating strict filtering (default read from INI file)" << endl
				 << "  -sequences_file	          Filename of a fasta file containing protein sequences. "
				 << "All peptides that are not a substring of" 
				 << " a sequence in this file are "
				 << "filtered out (default read from INI file)." << endl
				 << "  -pepfr                          the fraction of the peptide significance threshold that should be"
				 << " reached by a peptide hit" << endl
				 << "  -protfr                         the fraction of the protein significance threshold that should be"
				 << " reached by a protein hit" << endl
				 << "  -exclusion_peptides_file        if this AnalysisXML file is given all peptides having the same"
				 << " sequence as any in the AnalysisXML file will be dropped" << endl
				 << endl ;
	}
		
	void printToolHelpOpt_()
	{
		cerr << endl
				 << tool_name_ << endl
				 << endl
				 << "INI options:" << endl
				 << "  in                                         input file" << endl
				 << "  out                                        output file" << endl
				 << "  strict                                     flag indicating strict filtering" << endl
				 << "  sequences_file                             Filename of a fasta file containing protein sequences. "
				 << "All peptides that are not a substring of" 
				 << " a sequence in this file are filtered out." << endl
				 << "  peptide_significance_threshold_fraction    the fraction of the peptide significance threshold that should be"
				 << " reached by a peptide hit" << endl
				 << "  protein_significance_threshold_fraction    the fraction of the protein significance threshold that should be"
				 << " reached by a protein hit" << endl
				 << "  exclusion_peptides_file                    if this AnalysisXML file is given all peptides having the same"
				 << " sequence as any in the AnalysisXML file will be dropped" << endl
				 << endl
				 << "INI File example section:" << endl
				 << "  <ITEM name=\"in\" value=\"input.analysisXML\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"out\" value=\"output.analysisXML\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"strict\" value=\"false\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"protein_significance_threshold_fraction\" value=\"1\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"peptide_significance_threshold_fraction\" value=\"1\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"sequences_file\" value=\"sequences.fasta\" type=\"string\"/>" << endl
				 << "  <ITEM name=\"exclusion_peptides_file\" value=\"peptides.analysisXML\" type=\"string\"/>" << endl; 					 
	}		

	void setOptionsAndFlags_()
	{
		options_["-out"] = "out";
		options_["-in"] = "in";
		options_["-pepfr"] = "peptide_significance_threshold_fraction";
		options_["-protfr"] = "protein_significance_threshold_fraction";
		options_["-strict"] = "strict";
		options_["-sequences_file"] = "sequences_file";
		options_["-exclusion_peptides_file"] = "exclusion_peptides_file";
	}

	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name;
		String outputfile_name;
		IDFilter filter;
		AnalysisXMLFile analysisXML_file;
		ContactPerson contact_person;
		ContactPerson contact_person_exclusion;
		vector<ProteinIdentification> protein_identifications;
		vector<Identification> identifications;
		vector<Identification> identifications_exclusion;
		vector<Real> precursor_retention_times;
		vector<Real> precursor_retention_times_exclusion;
		vector<Real> precursor_mz_values;
		vector<Real> precursor_mz_values_exclusion;
		vector<Identification> filtered_identifications;
		vector<Real> filtered_precursor_retention_times;
		vector<Real> filtered_precursor_mz_values;
		Identification filtered_identification;
		vector<UnsignedInt> charges;
		Real protein_significance_threshold_fraction = 1;
		Real peptide_significance_threshold_fraction = 1;
		bool strict = false;
		QFileInfo file_info;
		QFile file;
		String sequences_file_name = "";
		String exclusion_peptides_file_name = "";
		FASTAFile fasta_file;
		vector< pair< String, String > > sequences;
		map<String, double> predicted_retention_times;
		DoubleReal predicted_sigma;
		Real total_gradient_time = 0.f;
		DoubleReal allowed_deviation = 10;
		vector<String> exclusion_peptides;
							
		//input file names and types
		inputfile_name = getParamAsString_("in");			
		writeDebug_(String("Input file: ") + inputfile_name, 1);
		if (inputfile_name == "")
		{
			writeLog_("No input file specified. Aborting!");
			cout << "No input file specified. Aborting!" << endl;
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
	
		//output file names and types
		outputfile_name = getParamAsString_("out");
		writeDebug_(String("Output file: ") + outputfile_name, 1);
		if (outputfile_name == "")
		{
			writeLog_("No output file specified. Aborting!");
			cout << "No output file specified. Aborting!" << endl;
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}				

		contact_person.setName(getParamAsString_("contactName", "unknown"));
		writeDebug_(String("Contact name: ") + contact_person.getName(), 1);

		contact_person.setInstitution(getParamAsString_("contactInstitution", "unknown"));
		writeDebug_(String("Contact institution: ") + contact_person.getInstitution(), 1);
			
		contact_person.setContactInfo(getParamAsString_("contactInfo"));
		if (contact_person.getContactInfo() != "")
		{
			writeDebug_(String("Contact info: ") + contact_person.getContactInfo(), 1);
		}

		peptide_significance_threshold_fraction = getParamAsString_("peptide_significance_threshold_fraction", String("0.f")).toFloat();
		writeDebug_(String("Peptide significance threshold fraction: ") + 
								String(peptide_significance_threshold_fraction), 1);

		protein_significance_threshold_fraction = getParamAsString_("protein_significance_threshold_fraction", String("0.f")).toFloat();
		writeDebug_(String("Protein significance threshold fraction: ") + 
								String(protein_significance_threshold_fraction), 1);

		if (getParamAsString_("strict", "false") != "false")
		{				
			writeDebug_(String("strict filtering (if there is more than one best hit")
									+ String(" for one spectrum, discard all of them)"), 1);
			strict = true;
		}
		else
		{
			writeDebug_(String("no strict filtering (if there is more than one best hit")
									+ String(" for one spectrum, take all of them)"), 1);				
		}


		sequences_file_name = getParamAsString_("sequences_file");
		if (sequences_file_name != "")
		{
			writeDebug_(String("Filter sequences in file: ") + sequences_file_name, 1);
		}

		exclusion_peptides_file_name = getParamAsString_("exclusion_peptides_file");
		if (exclusion_peptides_file_name != "")
		{
			writeDebug_(String("Exclusion peptides File: ") + exclusion_peptides_file_name, 1);
		}
			
		if ((total_gradient_time 
				 = getParamAsString_("total_gradient_time", String("0.f")).toFloat())
				!= 0.f)
		{
			writeDebug_(String("Total gradient time: ") + String(total_gradient_time) +
									String(" used for filtering"), 1);
		}
			
		if ((allowed_deviation 
				 = getParamAsString_("allowed_deviation", "0.f").toFloat()) 
				!= 0.f)
		{
			writeDebug_(String("allowed deviation: ") + String(allowed_deviation), 1);
		}
			
		//-------------------------------------------------------------
		// testing whether input and output files are accessible
		//-------------------------------------------------------------
	
		file_info.setFile(inputfile_name.c_str());
		if (!file_info.exists())
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
		}
		if (!file_info.isReadable())
		{
			throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);			
		}
		if (file_info.size() == 0)
		{
			throw Exception::FileEmpty(__FILE__, __LINE__, __PRETTY_FUNCTION__, inputfile_name);
		}		
		if (sequences_file_name != "")
		{
			file_info.setFile(sequences_file_name.c_str());
			if (!file_info.exists())
			{
				sequences_file_name = "";
			}
			if (!file_info.isReadable())
			{
				sequences_file_name = "";
			}
			if (sequences_file_name != "")
			{ 
				writeDebug_("Sequences file <" + sequences_file_name + "> used for"
										+ " filtering", 1);
			}
		}
		file.setName(outputfile_name.c_str());
		file.open( IO_WriteOnly );
		if (!file.isWritable())
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, outputfile_name);
		}
		file.close();				
	
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
	
		if (total_gradient_time != 0.f)
		{
			analysisXML_file.load(inputfile_name,
														&protein_identifications, 
														&identifications, 
														&precursor_retention_times, 
														&precursor_mz_values,
														&contact_person,
														&predicted_retention_times,
														&predicted_sigma);
		}
		else
		{
			analysisXML_file.load(inputfile_name, 
														&protein_identifications, 
														&identifications, 
														&precursor_retention_times, 
														&precursor_mz_values,
														&contact_person);				
		}
		if (sequences_file_name != "")
		{
			fasta_file.load(sequences_file_name,sequences);				
		}
			
		if (exclusion_peptides_file_name  != "")
		{
			analysisXML_file.load(exclusion_peptides_file_name, 
													  &protein_identifications, 
													 	&identifications_exclusion, 
													 	&precursor_retention_times_exclusion, 
													 	&precursor_mz_values_exclusion,
													 	&contact_person_exclusion);
			for(UnsignedInt i = 0; i < identifications_exclusion.size(); i++)
			{
				for(vector<PeptideHit>::const_iterator it = identifications_exclusion[i].getPeptideHits().begin();
						it != identifications_exclusion[i].getPeptideHits().end();
						it++)
				{
					exclusion_peptides.push_back(it->getSequence());
				}
			} 
		}												 
			
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
						
		/// Filtering identifications	by thresholds
		for(UnsignedInt i = 0; i < identifications.size(); i++)
		{	
																 
			filter.filterIdentificationsByThresholds(identifications[i], 
																							 peptide_significance_threshold_fraction, 
																							 protein_significance_threshold_fraction,
																							 filtered_identification, 
																							 strict);
			if (sequences_file_name != "")
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByProteins(temp_identification, 
																							 sequences,
																							 filtered_identification);
			}
				
			if (total_gradient_time != 0.f)
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByRetentionTimes(temp_identification, 
																										 predicted_retention_times,
																										 precursor_retention_times[i],
																										 predicted_sigma,
																										 allowed_deviation,
																										 total_gradient_time,
																										 filtered_identification);
			}
				
			if (exclusion_peptides_file_name != "")
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByExclusionPeptides(temp_identification,
																												exclusion_peptides,
																												filtered_identification); 				
			}
	
			if(!filtered_identification.empty())
			{
				filtered_identifications.push_back(filtered_identification);
				filtered_precursor_retention_times.push_back(precursor_retention_times[i]);
				filtered_precursor_mz_values.push_back(precursor_mz_values[i]);
			}
		}
						
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		if (total_gradient_time != 0.f)
		{
			analysisXML_file.store(outputfile_name,
														 protein_identifications, 
												 		 filtered_identifications, 
												 		 filtered_precursor_retention_times, 
												 		 filtered_precursor_mz_values,
												 		 contact_person,
												 		 predicted_retention_times,
												 		 predicted_sigma);
		}
		else
		{
			analysisXML_file.store(outputfile_name,
														 protein_identifications, 
												 		 filtered_identifications, 
												 		 filtered_precursor_retention_times, 
												 		 filtered_precursor_mz_values,
												 		 contact_person);
		}
		return OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPIDFilter tool;

	return tool.main(argc,argv);
}

/// @endcond
