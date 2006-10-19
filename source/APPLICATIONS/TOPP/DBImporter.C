// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DBAdapter.h>

#include "TOPPBase.h"

#include <qapplication.h>

using namespace OpenMS;
using namespace std;

/**
	@page DBImporter DBImporter
	
	@brief Imports an mzData file to an OpenMS database.
	
	Besides the file to import only the connection data has to be given.
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES 

class TOPPDBImporter
	: public TOPPBase
{
	public:
		TOPPDBImporter()
			: TOPPBase("DBImporter")
		{
			
		}
	
	protected:
		void printToolUsage_() const
		{
			cerr << endl
		       << getToolName() << " -- Imports an mzData file to an OpenMS database." << endl
		       << "Version: " << VersionInfo::getVersion() << endl
		       << endl
		       << "Usage:" << endl
					 << " " << getToolName() << " [options]" << endl
					 << endl
					 << "Options are:" << endl
					 << "  -u <DB user>      user/login of the DB" << endl
					 << "  -h <DB host>      host name of the DB server (default: localhost)" << endl
					 << "  -p <DB password>  password on the DB" << endl
					 << "  -P <DB port>      port the DB server is running on (default: 3306)" << endl			 
					 << "  -db <DB name>     DB name" << endl
					 << "  -in <file>        input file in mzData format" << endl;
		}
	
		void printToolHelpOpt_() const
		{
			cerr << endl
		       << getToolName() << endl
		       << endl
		       << "INI options:" << endl
					 << "  user      user/login of the DB" << endl
					 << "  host      host name of the DB server (default: localhost)" << endl
					 << "  password  password on the DB" << endl
					 << "  port      port the DB server is running on (default: 3306)" << endl			 
					 << "  db        DB name" << endl
					 << "  in        input file in mzData format" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"user\" value=\"user\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"host\" value=\"192.168.0.16\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"password\" value=\"password\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"port\" value=\"3307\" type=\"int\"/>" << endl
					 << "  <ITEM name=\"db\" value=\"OpenMS_DB\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-p"] = "password";
			options_["-u"] = "user";
			options_["-h"] = "host";
			options_["-P"] = "port";
			options_["-db"] = "db";
			options_["-in"] = "in";
		}
	
		ExitCodes main_(int argc , char** argv)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			
			//varaibles
			String db,user,password,host,in;
			SignedInt port;
			
			//input file names and types
			in = getParamAsString_("in");			
			writeDebug_(String("Input file: ") + in, 1);
			
			//db
			db = getParamAsString_("db");
			writeDebug_(String("db: ") + db,1);	
	
			//user
			user = getParamAsString_("user");
			writeDebug_(String("user: ") + user,1);	
	
			//password
			password = getParamAsString_("password");
			writeDebug_(String("password: ") + password,5);	
	
			//host
			host = getParamAsString_("host","localhost");
			writeDebug_(String("host: ") + host,1);	
	
			//port
			port = getParamAsInt_("port",3306);
			writeDebug_(String("port: ") + String(port),1);	
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			//load input file data
			MSExperiment< > exp;
			MzDataFile f;
			
			f.load(in,exp);			
			
			QApplication app(argc,argv,false);
			
			DBConnection con;
			con.connect(db, user, password, host, port);
			DBAdapter a(con);
			
			a.storeExperiment(exp);
			
			writeLog_( String(" written file to DB (id: ") + String((double)(exp.getPersistenceId())) + ")");	
			
			return OK;
		}
};

/// @endcond

int main( int argc, char ** argv )
{
	TOPPDBImporter tool;
	return tool.main(argc,argv);
}


