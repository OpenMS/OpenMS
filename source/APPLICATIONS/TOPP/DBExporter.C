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
	@page DBExporter DBExporter
	
	@brief Extracts MS data from a OpenMS database
	
	Extracts arbitrary MS data (MS, LC-MS, MS/MS) from a OpenMS database.
	A single dataset can be exported by giving one id contained in the 'MSExperiment' table.
	A query that returns several ids of the 'MSExperiment' table can be used to export several datasets at a time.
	
	If only one dataset is exported, it is stored with the given name.
	If several datasets are exported, the given name is prefixed with the DB id and a underscore.
	
	@ingroup TOPP
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES 

class TOPPDBExporter
	: public TOPPBase
{
	public:
		TOPPDBExporter()
			: TOPPBase("DBExporter")
		{
			
		}
	
	protected:
		void printToolUsage_() const
		{
			cerr << endl
       << getToolName() << " -- Extracts MS data from a OpenMS database." << endl
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
			 << "  -id <DB id>       id of the the map to export" << endl
			 << "  -query <query>    a SQL query that returns one or several DB ids of the MSExperiment table" << endl
			 << "  -out <file>       output file in mzData format (prefixed with DB id and '_' if several files are exported)" << endl;
		}
	
		void printToolHelpOpt_() const
		{
			cerr << endl
		       << getToolName() << endl
		       << endl
		       << "INI options:" << endl
					 << "  user     user/login of the DB" << endl
					 << "  host     host name of the DB server (default: localhost)" << endl
					 << "  password password on the DB" << endl
					 << "  port     port the DB server is running on (default: 3306)" << endl			 
					 << "  db       DB name" << endl
					 << "  id       DB id of the data in the table MSExperiment" << endl
					 << "  query    a SQL query that returns one or several DB ids of the MSExperiment table" << endl
					 << "  out      output file name (stored in mzData format)" << endl
					 << endl
					 << "INI File example section:" << endl
					 << "  <ITEM name=\"user\" value=\"user\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"host\" value=\"192.168.0.16\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"password\" value=\"password\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"port\" value=\"3307\" type=\"int\"/>" << endl
					 << "  <ITEM name=\"db\" value=\"OpenMS_DB\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"id\" value=\"1234567\" type=\"int\"/>" << endl
					 << "  <ITEM name=\"query\" value=\"SELECT id from MSExperiment WHERE 1\" type=\"string\"/>" << endl
					 << "  <ITEM name=\"out\" value=\"ouput.mzData\" type=\"string\"/>" << endl;
		}
	
		void setOptionsAndFlags_()
		{
			options_["-p"] = "password";
			options_["-u"] = "user";
			options_["-h"] = "host";
			options_["-P"] = "port";
			options_["-db"] = "db";
			options_["-id"] = "id";
			options_["-query"] = "query";
			options_["-out"] = "out";
		}
	
		ExitCodes main_(int argc , char** argv)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			
			//varaibles
			String db,user,password,host,query,out;
			UnsignedInt port;
			UID id;
			
			//input file names and types
			out = getParamAsString_("out");			
			writeDebug_(String("Output file: ") + out, 1);
			
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

			//query
			query = getParamAsString_("query");
			writeDebug_(String("query: ") + query,1);	
	
			//id
			id = getParamAsInt_("id",0);
			writeDebug_(String("id: ") + String(double(id)),1);	

	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			
			QApplication app(argc,argv,false);
			
			vector<UID> ids;
			
			if (id!=0)
			{
				ids.push_back(id);
			}
			
			if (query!="")
			{
				DBConnection con;
				con.connect(db, user, password, host, port);
				QSqlQuery result;
				con.executeQuery(query,result);
				while (result.isValid())
				{
					ids.push_back(result.value(0).toInt());
					result.next();
				}
			}
			
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			
			if (ids.size()>=1)
			{
				writeDebug_("Opening DB connection ...",1);	
				DBConnection con;
				con.connect(db, user, password, host, port);
				DBAdapter a(con);
				
				MzDataFile f;
				MSExperiment<> exp;
				
				if (ids.size()==1)
				{
					writeDebug_("Writing single file...",1);	
					//load from db
					a.loadExperiment(*(ids.begin()),exp);
					
					//write to file
					f.store(out, exp);
				}
				else
				{
					writeDebug_("Writing multiple files...",1);	
					for (vector<UID>::iterator it = ids.begin(); it!=ids.end(); ++it)
					{
						//load from db
						a.loadExperiment(*it,exp);
					
						//write to file
						stringstream ss;
						ss << *it;
						ss << "_" << out;
						f.store(ss.str(), exp);
					}
				}		
			}		
			
			return OK;
		}
};

/// @endcond

int main( int argc, char ** argv )
{
	TOPPDBExporter tool;
	return tool.main(argc,argv);
}



