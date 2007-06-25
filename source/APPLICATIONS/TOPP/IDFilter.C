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

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <limits>
#include <cmath>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page IDFilter_TOPP IDFilter
	
	@brief Filters Filters ProteinIdentification engine results by different criteria.
	
	This tool is used to filter the identifications found by
	an ProteinIdentification tool like Mascot&copy;. The identifications 
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
			sigma that is calculated and stored in the IdXMLFile corresponds to these normalized retention times.) 
			The maximum allowed deviation from the original
			retention time using the laplace error model that is learnt for confidently assigned peptides in RTModel. 
			It serves as a scaling of standard 
			deviations that are allowed for the predicted retention times. If set to 1 this means that one standard 
			deviation unit is allowed. 
			This filter can only be applied to IdXML files produced by RTPredict.
		</li>
		<li>
			<b>exclusion peptides</b>:<br> For this option you specify an IdXML file.
			All peptides that are present in both files (in-file and exclusion peptides
			file) will be dropped.
		</li>
		<li>
			<b>best hits only</b>:<br> Only the best hit of a spectrum is kept.
			If there is more than one hit for a spectrum with the maximal score then
			none of the hits will be kept.
		</li>
	</ul>
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFilter
	: public TOPPBase
{
 public:
	TOPPIDFilter()
		: TOPPBase("IDFilter","filters ProteinIdentification engine results by different criteria")
	{
		
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in IdXML format");
		registerStringOption_("out","<file>","","output file in IdXML format");	
		registerStringOption_("sequences_file","<file>","","filename of a fasta file containing protein sequences.\n"
																											 "All peptides that are not a substring of a sequence in this file are filtered out",false);
		registerStringOption_("exclusion_peptides_file","<file>","","An IdXML file. Peptides having the same sequence as any peptide in this file will be filtered out",false);
		registerDoubleOption_("pep_fraction","<fraction>",0.0,"the fraction of the peptide significance threshold that should be reached by a peptide hit",false);	
		registerDoubleOption_("prot_fraction","<fraction>",0.0,"the fraction of the protein significance threshold that should be reached by a protein hit",false);
		registerDoubleOption_("pep_score","<score>", 0,"the score which should be reached by a peptide hit to be kept",false);	
		registerDoubleOption_("prot_score","<score>", 0,"the score which should be reached by a protein hit to be kept",false);
		registerDoubleOption_("p_value","<significance>",0.05,"The probability of a correct ProteinIdentification having a deviation between observed and predicted rt equal or bigger than allowed",false);	
		registerIntOption_("best_n_peptide_hits","<score>", 0, "If this value is set only the n highest scoring peptide hits are kept.", false);
		registerIntOption_("best_n_protein_hits","<score>", 0, "If this value is set only the n highest scoring protein hits are kept.", false);
		registerFlag_("best_hits", "If this flag is set only the highest scoring hit is kept.\n"
															"If there are two or more highest scoring hits, none are kept.");
		registerFlag_("rt_filtering","If this flag is set rt filtering will be pursued.");
	}

	ExitCodes main_(int , char**)
	{
	
		//-------------------------------------------------------------
		// varaibles
		//-------------------------------------------------------------
			
		IDFilter filter;
		IdXMLFile IdXML_file;
		vector<ProteinIdentification> protein_identifications;
		vector<PeptideIdentification> identifications;
		vector<PeptideIdentification> identifications_exclusion;
		vector<PeptideIdentification> filtered_peptide_identifications;
		vector<ProteinIdentification> filtered_protein_identifications;
		PeptideIdentification filtered_identification;
		ProteinIdentification filtered_protein_identification;
		vector<UInt> charges;
		vector< pair< String, String > > sequences;
		vector<String> exclusion_peptides;
		bool rt_filtering = false;
		DoubleReal p_value = 0.05;
		
		
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name = getStringOption_("in");			
		inputFileReadable_(inputfile_name);
		
		String outputfile_name = getStringOption_("out");
		outputFileWritable_(outputfile_name);
		
		DoubleReal peptide_significance_threshold_fraction = getDoubleOption_("pep_fraction");
		DoubleReal protein_significance_threshold_fraction = getDoubleOption_("prot_fraction");
		DoubleReal peptide_threshold_score = getDoubleOption_("pep_score");
		DoubleReal protein_threshold_score = getDoubleOption_("prot_score");
		
		Int best_n_peptide_hits = getIntOption_("best_n_peptide_hits");
		Int best_n_protein_hits = getIntOption_("best_n_protein_hits");
		
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
		
		p_value = getDoubleOption_("p_value");
		rt_filtering = getFlag_("rt_filtering");

		bool strict = getFlag_("best_hits");
	
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
	
		if (rt_filtering)
		{
			IdXML_file.load(inputfile_name, protein_identifications, identifications);
		}
		else
		{
			IdXML_file.load(inputfile_name, protein_identifications, identifications);				
		}
		if (sequences_file_name != "")
		{
			FASTAFile().load(sequences_file_name,sequences);				
		}
			
		if (exclusion_peptides_file_name  != "")
		{
			IdXML_file.load(exclusion_peptides_file_name, protein_identifications, identifications_exclusion);
			for(UInt i = 0; i < identifications_exclusion.size(); i++)
			{
				for(vector<PeptideHit>::const_iterator it = identifications_exclusion[i].getHits().begin();
						it != identifications_exclusion[i].getHits().end();
						it++)
				{
					exclusion_peptides.push_back(it->getSequence());
				}
			} 
		}												 
			
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
						
		// Filtering peptide identifications	according to set criteria
		for(UInt i = 0; i < identifications.size(); i++)
		{
			if (fabs(peptide_significance_threshold_fraction - 0) < 0.00001)
			{
				filtered_identification = identifications[i];
			}
			else
			{
				filter.filterIdentificationsByThreshold(identifications[i], peptide_significance_threshold_fraction, filtered_identification);
			}
			if (sequences_file_name != "")
			{
				PeptideIdentification temp_identification = filtered_identification;				
				filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification);
			}

			if (rt_filtering)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByRTPValues(temp_identification, filtered_identification, p_value);																																
			}

			if (exclusion_peptides_file_name != "")
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByExclusionPeptides(temp_identification, exclusion_peptides, filtered_identification); 				
			}
	
			if (strict)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByBestHits(temp_identification, filtered_identification, strict); 				
			}
			
			if (setByUser_("pep_score"))
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByScore(temp_identification, peptide_threshold_score, filtered_identification); 				
			}

			if (setByUser_("best_n_peptide_hits"))
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByBestNHits(temp_identification, best_n_peptide_hits, filtered_identification); 				
			}

			if(!filtered_identification.getHits().empty())
			{
			  PeptideIdentification tmp;
			  tmp = filtered_identification;
			  tmp.setMetaValue("RT", identifications[i].getMetaValue("RT"));
			  tmp.setMetaValue("MZ", identifications[i].getMetaValue("MZ"));
				filtered_peptide_identifications.push_back(tmp);
			}
		}
						
		// Filtering protein identifications	according to set criteria
		for(UInt i = 0; i < protein_identifications.size(); i++)
		{
			if (fabs(protein_significance_threshold_fraction - 0) < 0.00001)
			{       
				filtered_protein_identification = protein_identifications[i];
			}
			else
			{
				filter.filterIdentificationsByThreshold(protein_identifications[i], protein_significance_threshold_fraction, filtered_protein_identification);
			}
			
			if (sequences_file_name != "")
			{
				ProteinIdentification temp_identification = filtered_protein_identification;				
				filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_protein_identification);
			}

			if (setByUser_("prot_score"))
			{
				ProteinIdentification temp_identification = filtered_protein_identification;
				filter.filterIdentificationsByScore(temp_identification, protein_threshold_score, filtered_protein_identification); 				
			}

			if (setByUser_("best_n_protein_hits"))
			{
				ProteinIdentification temp_identification = filtered_protein_identification;
				filter.filterIdentificationsByBestNHits(temp_identification, best_n_protein_hits, filtered_protein_identification); 				
			}

			if(!(filtered_protein_identification.getHits().empty()))
			{
				filtered_protein_identifications.push_back(filtered_protein_identification);
			}
		}
						
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		IdXML_file.store(outputfile_name, filtered_protein_identifications, filtered_peptide_identifications);

		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPIDFilter tool;

	return tool.main(argc,argv);
}

/// @endcond
