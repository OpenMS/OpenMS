#// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IDFilter IDFilter
	
	@brief Filters Filters identification engine results by different criteria.
	
	This tool is used to filter the identifications found by
	an identification tool like Mascot&copy;. The identifications 
	can be filtered by different criteria.
	
	<ul>
		<li> 
			<b>peptide significance threshold</b>:<br> This parameter 
			specifies which amount of the significance threshold should 
			be reached by a peptide to be kept. If for example a peptide
			has score 30 and the significance threshold is 40, the 
			peptide will only be kept by the filter if the significance 
			threshold fraction is set to 0.75 or lower.
		</li> 
		<li> 
			<b>protein significance threshold</b>:<br> This parameter 
			behaves in the same way as the peptide significance threshold 
			fraction parameter. The only difference is that it is used
			to filter protein hits.
		</li>
		<li>
			<b>peptide seqences</b>:<br> If you know which proteins
			are in the measured sample you can specify a FASTA file 
			which contains the protein sequences of those proteins. All 
			peptides which are not a substring of a protein contained
			in the sequences file will be filtered out. 
		</li>
		<li>
			<b>predicted retention time</b>:<br> To filter identifications according to their 
			predicted retention times you have to set two parameters:<br>  
			The total number of seconds that the gradient ran. (The model is learnt for normalized retention times and the 
			sigma that is calculated and stored in the analysisXMLFile corresponds to these normalized retention times.) 
			The maximum allowed deviation from the original
			retention time using the laplace error model that is learnt for confidently assigned peptides in RTModel. 
			It serves as a scaling of standard 
			deviations that are allowed for the predicted retention times. If set to 1 this means that one standard 
			deviation unit is allowed. 
			This filter can only be applied to AnalysisXML files produced by RTPredict.
		</li>
		<li>
			<b>exclusion peptides</b>:<br> For this option you specify an AnalysisXML file.
			All peptides that are present in both files (in-file and exclusion peptides
			file) will be dropped.
		</li>
		<li>
			<b>best hits only</b>:<br> Only the best hit of a spectrum is kept.
			If there is more than one hit for a spectrum with the maximal score then
			none of the hits will be kept.
		</li>
	</ul>
	
	@todo write test for filtering by retention time (Nico)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFilter
	: public TOPPBase
{
 public:
	TOPPIDFilter()
		: TOPPBase("IDFilter","filters identification engine results by different criteria")
	{
		
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in AnalysisXML format");
		registerStringOption_("out","<file>","","output file in AnalysisXML format");	
		registerStringOption_("sequences_file","<file>","","filename of a fasta file containing protein sequences.\n"
																											 "All peptides that are not a substring of a sequence in this file are filtered out",false);
		registerStringOption_("exclusion_peptides_file","<file>","","An AnalysisXML file. Peptides having the same sequence as any peptide in this file will be filtered out",false);
		registerDoubleOption_("pep_fraction","<fraction>",0.0,"the fraction of the peptide significance threshold that should be reached by a peptide hit",false);	
		registerDoubleOption_("prot_fraction","<fraction>",0.0,"the fraction of the protein significance threshold that should be reached by a protein hit",false);
		registerDoubleOption_("total_gradient_time","<time>",0.0,"the total time the HPLC gradient ran",false);	
		registerDoubleOption_("allowed_deviation","<dev>",0.0,"standard deviation allowed for the predicted retention times",false);	
		registerFlag_("best_hits","If this flag is set only the highest scoring hit is kept.\n"
															"If there is are two or more highest scoring hits, none are kept.");
	}

	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// varaibles
		//-------------------------------------------------------------
			
		IDFilter filter;
		AnalysisXMLFile analysisXML_file;
		vector<ProteinIdentification> protein_identifications;
		vector<IdentificationData> identifications;
		vector<IdentificationData> identifications_exclusion;
		vector<IdentificationData> filtered_identifications;
		Identification filtered_identification;
		vector<UnsignedInt> charges;
		vector< pair< String, String > > sequences;
		map<String, double> predicted_retention_times;
		DoubleReal predicted_sigma;
		vector<String> exclusion_peptides;
		
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name = getStringOption_("in");			
		inputFileReadable_(inputfile_name);
		
		String outputfile_name = getStringOption_("out");
		outputFileWritable_(outputfile_name);
		
		double peptide_significance_threshold_fraction = getDoubleOption_("pep_fraction");
		double protein_significance_threshold_fraction = getDoubleOption_("prot_fraction");
		
		String sequences_file_name = getStringOption_("sequences_file");
		if (sequences_file_name!="")
		{
			inputFileReadable_(sequences_file_name);
		}
		
		String exclusion_peptides_file_name = getStringOption_("exclusion_peptides_file");
		if (exclusion_peptides_file_name!="")
		{
			inputFileReadable_(exclusion_peptides_file_name);
		}
		
		double total_gradient_time = getDoubleOption_("total_gradient_time");
		double allowed_deviation = getDoubleOption_("allowed_deviation");

		bool strict = getFlag_("best_hits");
	
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
	
		if (total_gradient_time != 0.0)
		{
			analysisXML_file.load(inputfile_name,
														protein_identifications, 
														identifications,
														predicted_retention_times,
														predicted_sigma);
		}
		else
		{
			analysisXML_file.load(inputfile_name, 
														protein_identifications, 
														identifications);				
		}
		if (sequences_file_name != "")
		{
			FASTAFile().load(sequences_file_name,sequences);				
		}
			
		if (exclusion_peptides_file_name  != "")
		{
			analysisXML_file.load(exclusion_peptides_file_name, 
													  protein_identifications, 
													 	identifications_exclusion);
			for(UnsignedInt i = 0; i < identifications_exclusion.size(); i++)
			{
				for(vector<PeptideHit>::const_iterator it = identifications_exclusion[i].id.getPeptideHits().begin();
						it != identifications_exclusion[i].id.getPeptideHits().end();
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
			filter.filterIdentificationsByThresholds(identifications[i].id, 
																							 peptide_significance_threshold_fraction, 
																							 protein_significance_threshold_fraction,
																							 filtered_identification);
			if (sequences_file_name != "")
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByProteins(temp_identification, 
																							 sequences,
																							 filtered_identification);
			}
			
			if (total_gradient_time != 0.0)
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByRetentionTimes(temp_identification, 
																										 predicted_retention_times,
																										 identifications[i].rt,
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
	
			if (strict)
			{
				Identification temp_identification = filtered_identification;
				filter.filterIdentificationsByBestHits(temp_identification,
																							 filtered_identification,
																							 strict); 				
			}

			if(!filtered_identification.empty())
			{
			  IdentificationData tmp;
			  tmp.id = filtered_identification;
			  tmp.rt = identifications[i].rt;
			  tmp.mz = identifications[i].mz;
				filtered_identifications.push_back(tmp);
			}
		}
						
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		if (total_gradient_time != 0.0)
		{
			analysisXML_file.store(outputfile_name,
														 protein_identifications, 
												 		 filtered_identifications,
												 		 predicted_retention_times,
												 		 predicted_sigma);
		}
		else
		{
			analysisXML_file.store(outputfile_name,
														 protein_identifications, 
												 		 filtered_identifications);
		}
		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPIDFilter tool;

	return tool.main(argc,argv);
}

/// @endcond
