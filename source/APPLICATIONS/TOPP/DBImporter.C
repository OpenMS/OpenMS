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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_DBImporter DBImporter
	
	@brief Imports an mzML file to an %OpenMS database.
	
	Besides the file to import, only the connection data has to be given.
	The data can then be retrieved by the @ref TOPP_DBExporter.
	
	The @em init flag can be used to create a new %OpenMS database.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_DBImporter.cli
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES 

class TOPPDBImporter
	: public TOPPBase
{
	public:
		TOPPDBImporter()
			: TOPPBase("DBImporter","Imports data to an OpenMS database.")
		{
			
		}
	
	protected:
		void registerOptionsAndFlags_()
		{			
			registerStringOption_("user", "<user>", "", "user/login of the DB");
			registerStringOption_("host", "<host>", "localhost", "host name of the DB server", false);
			registerStringOption_("password", "<password>", "", "password for the user");
			registerIntOption_("port", "<port>", 3306, "port the DB server is running on", false);
			registerStringOption_("db", "<name>", "", "DB name");
			registerInputFile_("in", "<file>", "", "input file ", false);
			setValidFormats_("in",StringList::create("mzML"));
			registerFlag_("init", "Deletes all tables and sets up a new OpenMS database.\n"
														"The data of 'in' is not imported!");
		}
	
		ExitCodes main_(int , const char**)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			
			//varaibles
			String db,user,password,host,in;
			Int port;
			
			bool init = getFlag_("init");
			if (!init)
			{
				in = getStringOption_("in");
			}
			
			db = getStringOption_("db");
			user = getStringOption_("user");
			password = getStringOption_("password");
			host = getStringOption_("host");
			port = getIntOption_("port");
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			DBConnection con;
			con.connect(db, user, password, host, port);
			DBAdapter a(con);
			
			if (init)
			{
				a.createDB();
			}
			else
			{
				//load input file data
				MSExperiment<Peak1D> exp;
				MzMLFile f;
				f.setLogType(log_type_);
				f.load(in,exp);			
				
				//annotate output with data processing info
				addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));
				
				//store data
				a.storeExperiment(exp);
			
				writeLog_( String(" written file to DB (id: ") + (double)(exp.getPersistenceId()) + ")");	
			}
			
			return EXECUTION_OK;
		}
};

/// @endcond

int main( int argc, const char** argv )
{
	TOPPDBImporter tool;
	return tool.main(argc,argv);
}


