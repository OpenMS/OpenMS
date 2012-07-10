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
// $Maintainer: Sven Nahnsen $
// $Authors: Sven Nahnsen, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_DecoyDatabase DecoyDatabase
	
	@brief Create decoy peptide databases from normal ones.
	
  Decoy databases are useful to control false discovery rates and thus estimate score cutoffs for identified spectra.
  
  The decoy can either be generated from reversed or shuffled sequences.
  
  To get a 'contaminants' database have a look at http://www.thegpm.org/crap/index.html or find/create your own contaminant database.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_DecoyDatabase.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude UTILS_DecoyDatabase.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDecoyDatabase
	: public TOPPBase
{
	public:
		TOPPDecoyDatabase()
			: TOPPBase("DecoyDatabase","Create decoy peptide databases from normal ones.", false)
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input FASTA file containing the database.");
      setValidFormats_("in", StringList::create("FASTA"));
			registerOutputFile_("out","<file>","","Output FASTA file where the decoy database will be written to.");
      setValidFormats_("out", StringList::create("FASTA"));
			registerStringOption_("decoy_string", "<string>", "_rev", "String that is appended to the accession of the protein database to indicate a decoy protein.", false);
      registerStringOption_("decoy_string_position", "<enum>", "suffix", "Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?", false);
      setValidStrings_("decoy_string_position", StringList::create("prefix,suffix"));
      registerFlag_("append", "If this flag is used, the decoy database is appended to the target database, allowing combined target decoy searches.");
			registerFlag_("shuffle","If 'true' then the decoy hit are shuffled from the target sequences, otherwise they are reversed");
			registerInputFile_("contaminants","<file>","","Input a FASTA file containing contaminants - if given they are included in the database (recommended)", false);
		}

    String getIdentifier_(const String& identifier, const String& decoy_string, const bool as_prefix)
    {
      if (as_prefix) return decoy_string + identifier;
      else return identifier + decoy_string;
    }

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String cont(getStringOption_("contaminants"));
			String out(getStringOption_("out"));
			bool append = getFlag_("append");
			bool shuffle = getFlag_("shuffle");
			
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			vector<FASTAFile::FASTAEntry> proteins;
			FASTAFile().load(in, proteins);

      if (!cont.empty())
			{
        vector<FASTAFile::FASTAEntry> contaminants;
				FASTAFile().load(cont, contaminants);
				Size num_contaminants = contaminants.size();
				for (Size k = 0; k < num_contaminants; ++k)
				{
					proteins.push_back(contaminants[k]);
				}
			}
			else
			{
				LOG_WARN << "Warning, no contaminant sequences have been included!"<<endl;
			}
			
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------

			String decoy_string(getStringOption_("decoy_string"));
      bool decoy_string_position_prefix =	(String(getStringOption_("decoy_string_position")) == "prefix" ? true : false);
			Size num_proteins = proteins.size();
			set<String> identifiers;
			if (shuffle)
			{
				for (Size i = 0; i < num_proteins; ++i)
				{
					if (identifiers.find(proteins[i].identifier) != identifiers.end())
					{
						LOG_WARN << "DecoyDatabase: Warning, identifier is not unique to sequence file: '" << proteins[i].identifier << "'!" << endl;
					}
					identifiers.insert(proteins[i].identifier);

					FASTAFile::FASTAEntry entry = proteins[i];

					String pro_seq, temp;
					pro_seq = entry.sequence;
					Size x = pro_seq.size();
					srand(time(0));
					while(x != 0)
					{
						Size y = rand()%x;
						temp +=pro_seq[y];
						pro_seq[y] = pro_seq[x-1];
						--x;
					}
					entry.sequence = temp;

					if (append)
					{
						entry.identifier = getIdentifier_(entry.identifier, decoy_string, decoy_string_position_prefix);
						proteins.push_back(entry);
					}
					else
					{
						proteins[i].sequence = entry.sequence;
						proteins[i].identifier = getIdentifier_(proteins[i].identifier, decoy_string, decoy_string_position_prefix);
					}
				}
			}
			else // !shuffle
			{
				for (Size i = 0; i < num_proteins; ++i)
				{
					if (identifiers.find(proteins[i].identifier) != identifiers.end())
					{
						LOG_WARN << "DecoyDatabase: Warning, identifier is not unique to sequence file: '" << proteins[i].identifier << "'!" << endl;
					}
					identifiers.insert(proteins[i].identifier);

					if (append)
					{
						FASTAFile::FASTAEntry entry = proteins[i];
						entry.sequence.reverse();
						entry.identifier = getIdentifier_(entry.identifier, decoy_string, decoy_string_position_prefix);
						proteins.push_back(entry);
					}
					else
					{
						proteins[i].sequence.reverse();
						proteins[i].identifier = getIdentifier_(proteins[i].identifier, decoy_string, decoy_string_position_prefix);
					}
				}
			}
			
			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
		
			FASTAFile().store(out, proteins);
	
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPDecoyDatabase tool;
	return tool.main(argc,argv);
}
  
/// @endcond





