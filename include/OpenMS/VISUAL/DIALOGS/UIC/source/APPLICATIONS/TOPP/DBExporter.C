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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_DBExporter DBExporter
	
	@brief Extracts MS data from a %OpenMS database.
	
	Extracts arbitrary MS data (MS, LC-MS, MS/MS) from a %OpenMS database.
	A single dataset can be exported by giving one id contained in the 'MSExperiment' table.
	A query that returns several ids of the 'MSExperiment' table can be used to export several datasets at a time.
	
	If only one dataset is exported, it is stored with the given name.
	If several datasets are exported, the given name is prefixed with the DB id and an underscore.
	
	In order to create a new %OpenMS database, use the @ref TOPP_DBImporter.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_DBExporter.cli
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES 

class TOPPDBExporter
	: public TOPPBase
{
	public:
		TOPPDBExporter()
			: TOPPBase("DBExporter","Exports data from an OpenMS database to a file.")
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
			registerIntOption_("id", "<DB id>", 0, "id of the the map to export", false);
			registerStringOption_("query", "<query>", "", "a SQL query that returns one or several DB ids of the MSExperiment table", false);
			registerStringOption_("out", "<file>", "", "output file in mzML format (prefixed with DB id and '_' if several files are exported)");
			addEmptyLine_();
			addText_("In order to create a new OpenMS database, please use the 'init' flag of the DBImporter.");
		}
	
		ExitCodes main_(int , const char**)
		{
	
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
			
			//varaibles
			String db,user,password,host,query,out;
			UInt port;
			UID id;
			
			out = getStringOption_("out");			
			db = getStringOption_("db");
			user = getStringOption_("user");
			password = getStringOption_("password");
			host = getStringOption_("host");
			port = getIntOption_("port");
			query = getStringOption_("query");
			id = getIntOption_("id");
			
			if (id==0 && query == "")
			{
				writeLog_("Error: You have to give weither the 'id' option or the 'query' option! Aborting.");
				return ILLEGAL_PARAMETERS;
			}
	
			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------
			vector<UID> ids;
			
			if (id!=0)
			{
				ids.push_back(id);
			}
			
			if (query!="")
			{
				DBConnection con;
				con.connect(db, user, password, host, port);
				QSqlQuery result = con.executeQuery(query);
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
				
				MzMLFile f;
				f.setLogType(log_type_);
				
				MSExperiment<Peak1D> exp;
				
				if (ids.size()==1)
				{
					writeDebug_("Writing single file...",1);	
					//load from db
					a.loadExperiment(*(ids.begin()),exp);

					//annotate output with data processing info
					addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

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
						
						//annotate output with data processing info
						addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));
						
						//write to file
						stringstream ss;
						ss << *it;
						ss << "_" << out;
						f.store(ss.str(), exp);
					}
				}		
			}		
			
			return EXECUTION_OK;
		}
};

/// @endcond

int main( int argc, const char** argv )
{
	TOPPDBExporter tool;
	return tool.main(argc,argv);
}



