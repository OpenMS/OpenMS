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
// $Maintainer: Andreas Bertsch$
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_GenericWrapper GenericWrapper
	
	@brief Allows generically the wrapping of external tools.
	
	@todo Write docu

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_GenericWrapper.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPGenericWrapper
	: public TOPPBase
{
	public:
		TOPPGenericWrapper()
			: TOPPBase("GenericWrapper", "Allows the generic wrapping of external tools.")
		{
		}
	
	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file ");
			//setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out", "<file>", "", "output file ");
	  	//setValidFormats_("out",StringList::create("mzML"));
			registerStringOption_("call", "<call>", "", "Command line which calls the external tool, e.g. 'ProteinProphet $ini $out'");

			addEmptyLine_();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String call(getStringOption_("call"));
			String logfile(getStringOption_("log"));
		
      //-------------------------------------------------------------
      // call external program
      //-------------------------------------------------------------

			writeDebug_("Original call: '" + call, 1);
			call.substitute("$in", in);
			call.substitute("$out", out);
			writeDebug_("Final call: '" + call , 1);

			Int status(system(call.c_str()));
      if (status != 0)
      {
        writeLog_("Error: External program problem! (Details can be seen in the logfile: \"" + logfile + "\")");
				writeLog_("Call was '" + call + "'");
        return EXTERNAL_PROGRAM_ERROR;
      }
			
			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPGenericWrapper tool;
	return tool.main(argc,argv);
}

