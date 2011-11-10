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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_SequenceCoverageCalculator SequenceCoverageCalculator
	
	@brief Prints information about IdXML files.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SequenceCoverageCalculator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSequenceCoverageCalculator
	: public TOPPBase
{
	public:
		TOPPSequenceCoverageCalculator()
			: TOPPBase("SequenceCoverageCalculator","Prints information about IdXML files.",false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in_database","<file>","","input file containing the database in fasta format");
			registerInputFile_("in_peptides","<file>","","input file containing the identified peptides");
		}

		void getStartAndEndIndex(const String& sequence, const String& substring, pair<Size, Size>& indices)
		{
			indices.first = 0;
			indices.second = 0;
			Size temp_index = 0;
			Size temp_count = 0;
			
			if (sequence.hasSubstring(substring))
			{
				for (Size i = 0; i <= sequence.size() - substring.size(); ++i)
				{
					temp_index = i;
					temp_count = 0;
					while(temp_index < sequence.size()								 
								&& temp_count < substring.size()
								&& sequence.at(temp_index) == substring.at(temp_index - i))
					{
						++temp_index;
						++temp_count;
					}
					if (temp_count == substring.size())
					{
						indices.first = i;
						indices.second = temp_index;
						i = sequence.size();
					}
				}
			}
		}

		ExitCodes main_(int , const char**)
		{
			IdXMLFile idXML_file;
			vector<ProteinIdentification> protein_identifications;
			vector<PeptideIdentification> identifications;
			String inputfile_name = "";
			String database_name = "";
			vector< FASTAFile::FASTAEntry > proteins;
			vector<DoubleReal> statistics;
			vector<Size> counts;
			vector<Size> mod_counts;
			vector<PeptideHit> temp_hits;
			vector<Size> coverage;
			Size spectrum_count = 0;
			map<String, Size> unique_peptides;
			map<String, Size> temp_unique_peptides;
			map<String, Size> temp_modified_unique_peptides;

			protein_identifications.push_back(ProteinIdentification());
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			inputfile_name = getStringOption_("in_peptides");			
			database_name = getStringOption_("in_database");			
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			String document_id;
			idXML_file.load(inputfile_name, protein_identifications, identifications, document_id);
			FASTAFile().load(database_name, proteins);				

			statistics.resize(proteins.size(), 0.);
			counts.resize(proteins.size(), 0);
			mod_counts.resize(proteins.size(), 0);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------					

			for (Size j = 0; j < proteins.size(); ++j)
			{
				coverage.clear();
				coverage.resize(proteins[j].sequence.size(), 0);
				temp_unique_peptides.clear();
				temp_modified_unique_peptides.clear();
				
				for (Size i = 0; i < identifications.size(); ++i)
				{
					if (!identifications[i].empty())
					{
						if (identifications[i].getHits().size() > 1)
						{
							cout << "Spectrum with more than one identification found, which is not allowed"
									 << endl << "use the IDFilter with the -best_hits option to filter for best hits." << endl;
							return ILLEGAL_PARAMETERS;
            }
            temp_hits.clear();
						identifications[i].getReferencingHits(proteins[j].identifier, temp_hits);

						if (temp_hits.size() == 1)
						{
							pair<Size, Size> indices;
							getStartAndEndIndex(proteins[j].sequence, temp_hits[0].getSequence().toUnmodifiedString(), indices);
							for (Size k = indices.first; k < indices.second; ++k)
							{
								coverage[k] = 1;
							}
							if (indices.first != indices.second)
							{
//								cout <<  temp_hits[0].getSequence().toUnmodifiedString() << endl;
							}
							++spectrum_count;
							if (unique_peptides.find(temp_hits[0].getSequence().toString()) == unique_peptides.end())
							{
								unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
							}
							if (temp_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_unique_peptides.end())
							{
								temp_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toUnmodifiedString(), 0));
							}
							if (temp_modified_unique_peptides.find(temp_hits[0].getSequence().toUnmodifiedString()) == temp_modified_unique_peptides.end())
							{
								temp_modified_unique_peptides.insert(make_pair(temp_hits[0].getSequence().toString(), 0));
							}
						}
					}					
				}
/*				cout << proteins[j].sequence << endl;
				for (Size k = 0; k < coverage.size(); ++k)
				{
					cout << coverage[k];
				}
				cout << endl;
*/
				//				statistics[j] = make_pair(, 
				//													accumulate(coverage.begin(), coverage.end(), 0) / proteins[j].sequence.size());
				statistics[j] = ((DoubleReal) accumulate(coverage.begin(), coverage.end(), Size(0))) / proteins[j].sequence.size();
				counts[j] = temp_unique_peptides.size();
				mod_counts[j] = temp_modified_unique_peptides.size();
//				cout << statistics[j] << endl;
			}

//			cout << "Sum of coverage is " << accumulate(statistics.begin(), statistics.end(), 0.) << endl;
			cout << "Average coverage per protein is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
			cout << "Average number of peptides per protein is " << (((DoubleReal) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
			cout << "Average number of un/modified peptides per protein is " << (((DoubleReal) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;
			cout << "Number of identified spectra: " << spectrum_count << endl;
			cout << "Number of unique identified peptides: " << unique_peptides.size() << endl;
			
			vector<DoubleReal>::iterator it = statistics.begin(); 
			vector<Size>::iterator it2 = counts.begin(); 
			vector<Size>::iterator it3 = mod_counts.begin(); 
			while(it != statistics.end())
			{
				if (*it == 0.)
				{
					it = statistics.erase(it);
					it2 = counts.erase(it2);
					it3 = mod_counts.erase(it3);
				}
				else
				{				
					++it;
					++it2;
					++it3;
				}
			}
			cout << "Average coverage per found protein (" << statistics.size() << ") is " << (accumulate(statistics.begin(), statistics.end(), 0.) / statistics.size()) << endl;
			cout << "Average number of peptides per found protein is " << (((DoubleReal) accumulate(counts.begin(), counts.end(), 0.)) / counts.size()) << endl;
			cout << "Average number of un/modified peptides per protein is " << (((DoubleReal) accumulate(mod_counts.begin(), mod_counts.end(), 0.)) / mod_counts.size()) << endl;
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPSequenceCoverageCalculator tool;
	return tool.main(argc,argv);
}
  
/// @endcond





