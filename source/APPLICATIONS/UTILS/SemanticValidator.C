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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_SemanticValidator SemanticValidator
	
	@brief SemanticValidator for analysisXML and mzML files.
	
	This util is able to validate analysisXML and mzML files
	using an instance document and a mapping file.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SemanticValidator.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude UTILS_SemanticValidator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSemanticValidator
	: public TOPPBase
{
 public:
	TOPPSemanticValidator()
		: TOPPBase("SemanticValidator","SemanticValidator for analysisXML and mzML files.",false)
	{
	}
	
 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file, either analysisXML or mzML.");
		registerInputFile_("mapping_file", "<file>", "", "Mapping file which is used to semantically validate the given XML file against this mapping file (see 'share/OpenMS/MAPPING' for templates).");
	}	
	
	ExitCodes main_(int , const char**)
	{
		String in_file = getStringOption_("in");
		String mapping_file = getStringOption_("mapping_file");

		CVMappings mappings;
		CVMappingFile().load(mapping_file, mappings, false);

		ControlledVocabulary cv;
		cv.loadFromOBO("PSI-MOD", File::find("/CHEMISTRY/PSI-MOD.obo"));
	  cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("brenda",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		cv.loadFromOBO("UNIMOD",File::find("/CV/unimod.obo"));
		cv.loadFromOBO("PSI-MS",File::find("/CV/psi-ms.obo"));

		// check cv params
		Internal::SemanticValidator semantic_validator(mappings, cv);
		semantic_validator.setCheckTermValueTypes(true);
		semantic_validator.setCheckUnits(true);

		StringList errors, warnings;

		/*bool valid =*/ semantic_validator.validate(in_file, errors, warnings);
    for (Size i=0; i<warnings.size(); ++i)
    {
    	cout << "Warning: " << warnings[i] << endl;
    }
    for (Size i=0; i<errors.size(); ++i)
    {
      cout << "Error: " << errors[i] << endl;
    }

		if (warnings.empty() && errors.empty())
		{
			cout << "Congratulations, the file is valid!" << endl;
		}

		return EXECUTION_OK;
	}
};

int main( int argc, const char** argv )
{
	TOPPSemanticValidator tool;
	return tool.main(argc,argv);
}

/// @endcond



