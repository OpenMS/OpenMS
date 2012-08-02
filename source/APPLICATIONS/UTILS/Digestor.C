// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeiffer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
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
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Digestor \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter (peptide blacklist)</td>
		</tr>
	</table>
</CENTER>

	This application is used to digest a protein database to get all
	peptides given a cleavage enzyme. At the moment only trypsin is supported.
	
  The output can be used as a blacklist filter input to @ref TOPP_IDFilter, to remove certain peptides.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_Digestor.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude UTILS_Digestor.html
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
      setValidFormats_("in",StringList::create("FASTA"));
			registerOutputFile_("out","<file>","","Output file (peptides)");
      setValidFormats_("out", StringList::create("idXML,FASTA"));
      registerStringOption_("out_type", "<type>","", "Set this if you cannot control the filename of 'out', e.g., in TOPPAS.", false);
      setValidStrings_("out_type", StringList::create("idXML,FASTA"));

      registerIntOption_("missed_cleavages","<number>",1,"The number of allowed missed cleavages", false);
			setMinInt_("missed_cleavages", 0);
			registerIntOption_("min_length","<number>",6,"Minimum length of peptide", false);
			registerIntOption_("max_length","<number>",40,"Maximum length of peptide", false);
			registerStringOption_("enzyme","<string>","Trypsin","The type of digestion enzyme", false);
		}

		ExitCodes main_(int , const char**)
		{
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			PeptideIdentification peptide_identification;
			ProteinIdentification protein_identification;
			
			protein_identifications.push_back(ProteinIdentification());
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String inputfile_name = getStringOption_("in");			
			String outputfile_name = getStringOption_("out");	

		  //input file type
		  FileHandler fh;
		  FileTypes::Type out_type = fh.nameToType(getStringOption_("out_type"));

		  if (out_type==FileTypes::UNKNOWN)
		  {
			  out_type = fh.getTypeByFileName(outputfile_name);
			  writeDebug_(String("Output file type: ") + fh.typeToName(out_type), 2);
		  }

		  if (out_type==FileTypes::UNKNOWN)
		  {
			  LOG_ERROR << ("Error: Could not determine output file type!") << std::endl;
			  return PARSE_ERROR;
		  }

			Size min_size = getIntOption_("min_length");
      Size max_size = getIntOption_("max_length");
			Size missed_cleavages = getIntOption_("missed_cleavages");
			

      bool has_FASTA_output = (out_type==FileTypes::FASTA);

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			std::vector<FASTAFile::FASTAEntry> protein_data;
			FASTAFile().load(inputfile_name, protein_data);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
		
			// This should be updated if more cleavage enzymes are available
			EnzymaticDigestion digestor;
			digestor.setEnzyme(EnzymaticDigestion::TRYPSIN);
			ProteinIdentification::SearchParameters search_parameters;
			search_parameters.enzyme = ProteinIdentification::TRYPSIN;
			digestor.setMissedCleavages(missed_cleavages);
			
      vector<String> protein_accessions(1);
			PeptideHit temp_peptide_hit;

			std::vector<FASTAFile::FASTAEntry> all_peptides;

      for (Size i = 0; i < protein_data.size(); ++i)
			{
        if (!has_FASTA_output)
        {
				  protein_accessions[0] = protein_data[i].identifier;
  			  ProteinHit temp_protein_hit;
          temp_protein_hit.setSequence(protein_data[i].sequence);
				  temp_protein_hit.setAccession(protein_accessions[0]);
				  protein_identifications[0].insertHit(temp_protein_hit);
  				temp_peptide_hit.setProteinAccessions(protein_accessions);
        }
				
   			vector<AASequence> temp_peptides;
				digestor.digest(AASequence(protein_data[i].sequence), temp_peptides);
	
        for (Size j = 0; j < temp_peptides.size(); ++j)
				{
					if ((temp_peptides[j].size() >= min_size) &&
              (temp_peptides[j].size() <= max_size))
					{
            if (!has_FASTA_output)
            {
						  temp_peptide_hit.setSequence(temp_peptides[j]);
						  peptide_identification.insertHit(temp_peptide_hit);
            }
            else // for FASTA file output
            {             
              FASTAFile::FASTAEntry pep(protein_data[i].identifier, protein_data[i].description, temp_peptides[j].toString());
              all_peptides.push_back(pep);
            }
					}
				}
			}

      if (!has_FASTA_output)
      {
			  DateTime date_time = DateTime::now();
			  String date_time_string = "";
  			
			  date_time_string = date_time.get();
			  protein_identifications[0].setSearchParameters(search_parameters);
			  protein_identifications[0].setDateTime(date_time);
			  protein_identifications[0].setSearchEngine("In-silico digestion");
			  protein_identifications[0].setIdentifier("In-silico_digestion" + date_time_string);
			  peptide_identification.setIdentifier("In-silico_digestion" + date_time_string);
			  identifications.push_back(peptide_identification);
      }
  			
			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
      
      if (has_FASTA_output)
      {
        FASTAFile().store(outputfile_name, all_peptides);
      }
      else
      {
			  IdXMLFile().store(outputfile_name,
										      protein_identifications,
											    identifications);
      }

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPDigestor tool;
	return tool.main(argc,argv);
}
  
/// @endcond





