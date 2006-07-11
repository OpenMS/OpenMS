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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DBConnection.h>
#include <qapplication.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DBConnection, "$Id: DBConnection_test.C,v 1.4 2006/03/09 13:01:19 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// start QApplication. Otherwise the DB classes are not available
QApplication qapp(argc,argv,false);

DBConnection* ptr = 0;
CHECK(DBConnection())
	ptr = new DBConnection();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DBConnection())
	delete ptr;
RESULT

// check for credentials
// if they are not present, abort the test (successfully)
bool do_tests=true;
TextFile credentials;
try
{
	credentials.load("DB_credentials.txt",true);
}
catch(...)
{
	do_tests=false;
}

String db,host,user,password,port;

//read out connection data
for (TextFile::iterator it = credentials.begin(); it!= credentials.end(); ++it)
{
	//comments and empty lines
	if (it->hasPrefix('#') || *it == "")
	{
		continue;
	}
	
	//extract connection info
	if (it->hasPrefix("Host:")) host = it->suffix(':').trim();
	if (it->hasPrefix("Port:")) port = it->suffix(':').trim();
	if (it->hasPrefix("User:")) user = it->suffix(':').trim();
	if (it->hasPrefix("Password:")) password = it->suffix(':').trim();
	if (it->hasPrefix("DB:")) db = it->suffix(':').trim();
}

if (do_tests)
{
	CHECK((void connect(const std::string& db, const std::string& user, const std::string& password, const std::string& host = "localhost", UnsignedInt port=3306, const std::string& QTDBDriver = DB_PLUGIN ) throw(InvalidQuery)))
	  DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EXCEPTION(DBConnection::InvalidQuery,con.connect("doesnotexist",user,password,host, port.toInt()))
	RESULT

	CHECK(std::string DBName() const)
		DBConnection con;
	  TEST_EQUAL(con.DBName(),"");
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.DBName(),db);
	RESULT

	CHECK(bool isConnected() const)
		DBConnection con;
	  TEST_EQUAL(con.isConnected(),false);
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.isConnected(),true);
	RESULT

	CHECK(void disconnect())
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.disconnect();
	  TEST_EQUAL(con.DBName(),"");
	  TEST_EQUAL(con.isConnected(),false);
	RESULT

	CHECK((void executeQuery(const std::string& query) throw(InvalidQuery, NotConnected)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.executeQuery("DROP TABLE IF EXISTS Dummy");
	  con.executeQuery("CREATE TABLE Dummy (id int,text varchar(5))");
	RESULT
	
	CHECK(const std::string& lastQuery() const)
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.lastQuery(),"");
	  con.executeQuery("INSERT INTO Dummy values (5,'bla'),(4711,'bluff')");
	  TEST_EQUAL(con.lastQuery(),"INSERT INTO Dummy values (5,'bla'),(4711,'bluff')");
	RESULT

	CHECK(std::string lastError() const)
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.lastError(),"");
		try
		{
	  	con.executeQuery("INSERT INTOOO Dummy values (5,lsdkh,sdfjsdf)");
	  }
		catch(...)
		{
			
		}
	  TEST_EQUAL(con.lastError()=="",false);
	RESULT

	CHECK(QSqlQuery& lastResult())
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		con.executeQuery("SELECT * FROM Dummy");
		TEST_EQUAL(con.lastResult().size(),2);
		con.lastResult().first();
 		TEST_EQUAL(con.lastResult().value(0).toString(),"5");
 		TEST_EQUAL(con.lastResult().value(1).toString(),"bla");
 		con.lastResult().next();
 		TEST_EQUAL(con.lastResult().value(0).toString(),"4711");
 		TEST_EQUAL(con.lastResult().value(1).toString(),"bluff");
	RESULT

	CHECK(String getStringValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL("bla",con.getStringValue("Dummy","text","5"));
		TEST_EQUAL("bluff",con.getStringValue("Dummy","text","4711"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getStringValue("Dummy2","text56","4711"))
		//TODO test ConversionError 
	RESULT

	CHECK(SignedInt getIntValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL(5,con.getIntValue("Dummy","id","5"));
		TEST_EQUAL(4711,con.getIntValue("Dummy","id","4711"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getIntValue("Dummy2","text56","4711"))
		//TODO test ConversionError 
	RESULT

	CHECK(double getDoubleValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_REAL_EQUAL(5.0,con.getIntValue("Dummy","id","5"));
		TEST_REAL_EQUAL(4711.0,con.getIntValue("Dummy","id","4711"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getIntValue("Dummy2","text56","4711"))
		//TODO test ConversionError 
	RESULT

	CHECK(UnsignedInt getId(const std::string& table, const std::string& column, const std::string& value) throw (InvalidQuery,NotConnected))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL(5,con.getId("Dummy","text","bla"));
		TEST_EQUAL(4711,con.getId("Dummy","text","bluff"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getId("Dummy2","text56","4711"))		
	RESULT

	CHECK((void render( std::ostream& out=std::cout, const std::string& separator=" | ", const std::string& line_begin="", const std::string& line_end="\n")))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		con.executeQuery("SELECT * FROM Dummy");
		stringstream s;
		con.render(s,"|",">","<");
		TEST_EQUAL(s.str(),">id|text<>5|bla<>4711|bluff<")
		stringstream s2;
		con.render(s2,"x","","; ");
		TEST_EQUAL(s2.str(),"idxtext; 5xbla; 4711xbluff; ")
	RESULT

	CHECK((template<class StringListType> void executeQueries(const StringListType& queries) throw(InvalidQuery, NotConnected)))
	  vector<String> qs;
	  qs.push_back("DROP TABLE IF EXISTS Dummy");
	  qs.push_back("CREATE TABLE Dummy (id int,text varchar(5))");
	  qs.push_back("INSERT INTO Dummy values (1,'bla'),(2,'bluff')");
	  qs.push_back("INSERT INTO Dummy values (3,'bla2'),(4,'bluff2')");
	  qs.push_back("DELETE FROM Dummy where id>2");

		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.executeQueries(qs);
	  
	  con.executeQuery("SELECT * FROM Dummy");
		stringstream s2;
		con.render(s2,"x","",";");
		TEST_EQUAL(s2.str(),"idxtext;1xbla;2xbluff;")	  	  
	RESULT

	CHECK(Deleting table 'Dummy')
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.executeQuery("DROP TABLE IF EXISTS Dummy");
	RESULT
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



