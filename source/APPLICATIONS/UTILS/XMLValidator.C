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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_XMLValidator XMLValidator

	@brief Validates XML files against an XSD schema.

	When a schema file is given, the input file is simply validated against the schema.

	When no schema file is given, the tool tries to determine the file type and validates the file against
	the latest schema version.

	@note XML schema files for the %OpenMS XML formats and several other XML formats can be found in the folder
	      OpenMS/share/OpenMS/SCHEMAS/

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_XMLValidator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPXMLValidator
	: public TOPPBase
{
 public:
	TOPPXMLValidator()
		: TOPPBase("XMLValidator","Validates XML files against an XSD schema.",false)
	{
	}

 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","file to validate");
		registerInputFile_("schema","<file>","","schema to validate against.\nIf no schema is given, the file is validated against the latest schema of the file type.", false);
	}

	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in");
		String schema = getStringOption_("schema");
		bool valid = true;

		if (schema!="") //schema explicitly given
		{
			XMLValidator xmlv;
			valid = xmlv.isValid(in,schema);
		}
		else //no schema given
		{
			//determine input type
			FileTypes::Type in_type = FileHandler::getType(in);
			if (in_type==FileTypes::UNKNOWN)
			{
				writeLog_("Error: Could not determine input file type!");
				return PARSE_ERROR;
			}

			cout << endl << "Validating " << FileHandler::typeToName(in_type) << " file";
			switch(in_type)
			{
				case FileTypes::MZDATA :
					cout << " against schema version " << MzDataFile().getVersion() << endl;
					valid = MzDataFile().isValid(in);
					break;
				case FileTypes::FEATUREXML :
					cout << " against schema version " << FeatureXMLFile().getVersion() << endl;
					valid = FeatureXMLFile().isValid(in);
					break;
				case FileTypes::IDXML :
					cout << " against schema version " << IdXMLFile().getVersion() << endl;
					valid = IdXMLFile().isValid(in);
					break;
				case FileTypes::CONSENSUSXML :
					cout << " against schema version " << ConsensusXMLFile().getVersion() << endl;
					valid = ConsensusXMLFile().isValid(in);
					break;
				case FileTypes::MZXML :
					cout << " against schema version " << MzXMLFile().getVersion() << endl;
					valid = MzXMLFile().isValid(in);
					break;
				case FileTypes::INI :
					cout << " against schema version " << Param().getVersion() << endl;
					valid = Param().isValid(in);
					break;
			  case FileTypes::PEPXML :
					cout << " against schema version " << PepXMLFile().getVersion() << endl;
					valid = PepXMLFile().isValid(in);
					break;
				default:
					cout << endl << "Aborted: Validation of this file type is not supported!" << endl;
					return EXECUTION_OK;
			};
		}

		//Result
		if (valid)
		{
			cout << "Success: the file is valid!" << endl;
		}
		else
		{
			cout << "Failed: errors are listed above!" << endl;
		}

		return EXECUTION_OK;
		}
};

int main( int argc, const char** argv )
{
	TOPPXMLValidator tool;
	return tool.main(argc,argv);
}

/// @endcond



