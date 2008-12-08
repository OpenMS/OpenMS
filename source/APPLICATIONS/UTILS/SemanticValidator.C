// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/CVMappings.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page SemanticValidator SemanticValidator
		
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SemanticValidator.cli
	
	@todo Docu (Andreas)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSemanticValidator
	: public TOPPBase
{
 public:
	TOPPSemanticValidator()
		: TOPPBase("SemanticValidator","TODO",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "");
		registerInputFile_("mapping_file", "<file>", "", "", false);
	}	
	
	ExitCodes main_(int , const char**)
	{
		String in_file = getStringOption_("in");
		String mapping_file = getStringOption_("mapping_file");

		CVMappings mappings;
		CVMappingFile().load(mapping_file, mappings, true);

		ControlledVocabulary cv;
		cv.loadFromOBO("PSI-PI", "psi-pi.obo");
		cv.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));
	  cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("brenda",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		cv.loadFromOBO("UNIMOD",File::find("/CV/unimod.obo"));
		//cv.loadFromOBO("NCBITaxon", "ncbi_taxonomy.obo");

		// check cv params
		Internal::SemanticValidator semantic_validator(mappings, cv);
		semantic_validator.setCheckTermValueTypes(true);
		semantic_validator.setCheckUnits(true);
		StringList errors, warnings;
		/*bool valid =*/ semantic_validator.validate(in_file, errors, warnings);
    for (UInt i=0; i<warnings.size(); ++i)
    {
    	cout << "Warning: " << warnings[i] << endl;
    }
    for (UInt i=0; i<errors.size(); ++i)
    {
      cout << "Error: " << errors[i] << endl;
    }

		if (warnings.size() == 0 && errors.size() == 0)
		{
			cout << "Congratulations, the file is valid!" << endl;
		}

						

		/*
		// check units cv
		Internal::SemanticValidator semantic_validator_u(mappings, cv);
		semantic_validator_u.setAccessionAttribute("unitAccession");
		semantic_validator_u.setNameAttribute("unitName");
		semantic_validator_u.setValueAttribute("value");
		semantic_validator_u.setCheckTermValueTypes(true);
		StringList unit_tags(StringList::create("MinusValue,PlusValue,MassDelta"));

		for (StringList::const_iterator it = unit_tags.begin(); it != unit_tags.end(); ++it)
		{
			semantic_validator_u.setTag(*it);
    	StringList errors_u, warnings_u;
    	bool valid_u = semantic_validator_u.validate(in_file, errors_u, warnings_u);
    	for (UInt i=0; i<warnings_u.size(); ++i)
    	{
      	cout << "Warning: " << warnings_u[i] << endl;
    	}
    	for (UInt i=0; i<errors_u.size(); ++i)
    	{
      	cout << "Error: " << errors_u[i] << endl;
    	}
		}*/
	
		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPSemanticValidator tool;
	return tool.main(argc,argv);
}

/// @endcond



