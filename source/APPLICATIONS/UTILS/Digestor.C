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

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_Digestor Digestor
	
	@brief Digests a protein database in-silico.
	
	This application is used to digest a protein database to get all
	peptides given a cleavage enzyme. At the moment only trypsin is supported.
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_Digestor.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDigestor
	: public TOPPBase
{
	public:
		TOPPDigestor()
			: TOPPBase("Digestor","Digests a protein database in-silico.",false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","input file");
			registerOutputFile_("out","<file>","","output file (peptides)\n");
			registerIntOption_("missed_cleavages","<number>",1,"the number of allowed missed cleavages", false);
			registerIntOption_("min_length","<number>",6,"minimum length of peptide", false);
			registerStringOption_("enzyme","<string>","Trypsin","the digestion enzyme", false);
			setMinInt_("missed_cleavages", 0);
		}

		ExitCodes main_(int , const char**)
		{
			IdXMLFile IdXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			vector< String > peptides;
			vector<PeptideHit> temp_peptide_hits;
			std::vector<FASTAFile::FASTAEntry> protein_data;
			FASTAFile file;
			EnzymaticDigestion digestor;
			vector<AASequence> temp_peptides;
			PeptideIdentification peptide_identification;
			ProteinIdentification protein_identification;
			PeptideHit temp_peptide_hit;
			ProteinHit temp_protein_hit;
			vector<String> protein_accessions;
			vector<String> parts;
			String inputfile_name = "";
			String outputfile_name = "";
			UInt min_size = 0;
			UInt missed_cleavages = 0;
			ProteinIdentification::SearchParameters search_parameters;
			
			protein_identifications.push_back(ProteinIdentification());
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			inputfile_name = getStringOption_("in");			
			outputfile_name = getStringOption_("out");	
			min_size = getIntOption_("min_length");
			missed_cleavages = getIntOption_("missed_cleavages");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			file.load(inputfile_name, protein_data);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			// This should be updated if more cleavage enzymes are available
			digestor.setEnzyme(EnzymaticDigestion::TRYPSIN);
			search_parameters.enzyme = ProteinIdentification::TRYPSIN;
			digestor.setMissedCleavages(missed_cleavages);
			
			protein_accessions.resize(1, String(""));
			for (Size i = 0; i < protein_data.size(); ++i)
			{
				protein_accessions[0] = protein_data[i].identifier;
				temp_protein_hit.setSequence(protein_data[i].sequence);
				temp_protein_hit.setAccession(protein_accessions[0]);
				
				digestor.digest(AASequence(protein_data[i].sequence), temp_peptides);
				temp_peptide_hit.setProteinAccessions(protein_accessions);
				for (Size j = 0; j < temp_peptides.size(); ++j)
				{
					if (temp_peptides[j].size() >= min_size)
					{
						temp_peptide_hit.setSequence(temp_peptides[j]);
						peptide_identification.insertHit(temp_peptide_hit);
					}
				}
				protein_identifications[0].insertHit(temp_protein_hit);
			}
			DateTime date_time = DateTime::now();
			String date_time_string = "";
			
			date_time_string = date_time.get();
			protein_identifications[0].setSearchParameters(search_parameters);
			protein_identifications[0].setDateTime(date_time);
			protein_identifications[0].setSearchEngine("In-silico digestion");
			protein_identifications[0].setIdentifier("In-silico_digestion" + date_time_string);
			peptide_identification.setIdentifier("In-silico_digestion" + date_time_string);
			identifications.push_back(peptide_identification);
			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			
			IdXML_file.store(outputfile_name,
											 protein_identifications,
											 identifications);

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPDigestor tool;
	return tool.main(argc,argv);
}
  
/// @endcond





