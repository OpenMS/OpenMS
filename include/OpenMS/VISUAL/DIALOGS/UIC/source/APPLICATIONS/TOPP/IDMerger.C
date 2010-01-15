// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDMerger IDMerger
	
	@brief Merges several IdXML files into one IdXML file.
	
	You can merge an unlimited number of files into one IdXML file.
	
	This tool is typically applied before @ref TOPP_ConsensusID or @ref TOPP_IDMapper.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDMerger.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMerger
	: public TOPPBase
{
 public:
	TOPPIDMerger()
		: TOPPBase("IDMerger","Merges several protein/peptide identification files into one file.")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>",StringList(),"two or more input files separated by blank");
		setValidFormats_("in",StringList::create("idXML"));
		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("idXML"));
	}
	
	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		StringList file_names = getStringList_("in");
		String out = getStringOption_("out");
		
		if (file_names.size() < 2)
		{
			writeLog_("Less than two filenames given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
				
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		vector<ProteinIdentification> protein_identifications;
		vector<PeptideIdentification> identifications;
		String document_id;
		IdXMLFile().load(file_names[0], protein_identifications, identifications, document_id);

		vector<String> used_ids;
		for (Size i=1; i<file_names.size(); ++i)
		{
			vector<ProteinIdentification> additional_protein_identifications;
			vector<PeptideIdentification> additional_identifications;
			IdXMLFile().load(file_names[i], additional_protein_identifications, additional_identifications, document_id);
			
			for (Size i=0; i<additional_protein_identifications.size();++i)
			{
				if (find(used_ids.begin(), used_ids.end(), additional_protein_identifications[i].getIdentifier())!=used_ids.end())
				{
					writeLog_(String("Error: The idenitifier '") + additional_protein_identifications[i].getIdentifier() + "' was used before!");
					return INCOMPATIBLE_INPUT_DATA;
				}
				used_ids.push_back(additional_protein_identifications[i].getIdentifier());
			}
			
			protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
			identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
		}
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		IdXMLFile().store(out, protein_identifications, identifications);
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
