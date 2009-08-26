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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>

#include <QtSql/QSqlQuery>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DBConnection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DBConnection* ptr = 0;
START_SECTION((DBConnection()))
	ptr = new DBConnection();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~DBConnection()))
	delete ptr;
END_SECTION

// check for credentials
// if they are not present, abort the test (successfully)
bool do_tests=true;
TextFile credentials;
try
{
	credentials.load(String(OPENMS_BINARY_PATH) + "/source/TEST/DB_credentials.txt",true);
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
START_SECTION((void connect(const String &db, const String &user, const String &password, const String &host="localhost", UInt port=3306, const String &QTDBDriver=DB_PLUGIN, const String &connection_name="OpenMS_default_connection")))
	  DBConnection con;
	  TEST_EXCEPTION(DBConnection::InvalidQuery,con.connect("doesnotexist",user,password,host, port.toInt()))
	  con.connect(db,user,password,host, port.toInt());
END_SECTION

START_SECTION((String DBName() const))
		DBConnection con;
	  TEST_EQUAL(con.DBName(),"");
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.DBName(),db);
END_SECTION

START_SECTION((bool isConnected() const))
		DBConnection con;
	  TEST_EQUAL(con.isConnected(),false);
	  con.connect(db,user,password,host, port.toInt());
	  TEST_EQUAL(con.isConnected(),true);
END_SECTION

START_SECTION((void disconnect()))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.disconnect();
	  TEST_EQUAL(con.DBName(),"");
	  TEST_EQUAL(con.isConnected(),false);
END_SECTION

START_SECTION((QSqlQuery executeQuery(const String &query, bool first=false)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.executeQuery("DROP TABLE IF EXISTS Dummy");
	  con.executeQuery("CREATE TABLE Dummy (id int,text varchar(5),number float )");
	  QSqlQuery result = con.executeQuery("INSERT INTO Dummy values (5,'bla','45.11'),(4711,'bluff','471.123')");
	  TEST_EQUAL(result.numRowsAffected(),2)
END_SECTION

START_SECTION((String getStringValue(const String &table, const String &column, const String &id)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL(con.getStringValue("Dummy","text","5"),"bla");
		TEST_EQUAL(con.getStringValue("Dummy","text","4711"),"bluff");
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getStringValue("Dummy2","text56","4711"))
		TEST_EXCEPTION(Exception::ConversionError, con.getStringValue("Dummy","id","sdfsdfsdf"))
END_SECTION

START_SECTION((Int getIntValue(const String &table, const String &column, const String &id)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL(5,con.getIntValue("Dummy","id","5"));
		TEST_EQUAL(4711,con.getIntValue("Dummy","id","4711"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getIntValue("Dummy2","text56","4711"))
		TEST_EXCEPTION(Exception::ConversionError, con.getIntValue("Dummy","text","sdfsdf"))
END_SECTION

START_SECTION((double getDoubleValue(const String &table, const String &column, const String &id)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_REAL_SIMILAR(45.11,con.getDoubleValue("Dummy","number","5"));
		TEST_REAL_SIMILAR(471.123,con.getDoubleValue("Dummy","number","4711"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getDoubleValue("Dummy2","text56","4711"))
		TEST_EXCEPTION(Exception::ConversionError, con.getDoubleValue("Dummy","text","sdfsdf"))
END_SECTION

START_SECTION((UInt getId(const String &table, const String &column, const String &value)))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		TEST_EQUAL(5,con.getId("Dummy","text","bla"));
		TEST_EQUAL(4711,con.getId("Dummy","text","bluff"));
		TEST_EXCEPTION(DBConnection::InvalidQuery, con.getId("Dummy2","text56","4711"))
END_SECTION

START_SECTION((void render(QSqlQuery &result, std::ostream &out=std::cout, const String &separator=" | ", const String &line_begin="", const String &line_end="\n")))
		DBConnection con;
		con.connect(db,user,password,host, port.toInt());
		QSqlQuery result = con.executeQuery("SELECT * FROM Dummy");
		stringstream s;
		con.render(result,s,"|",">","<");
		TEST_EQUAL(s.str(),">id|text|number<>5|bla|45.11<>4711|bluff|471.123<")
		stringstream s2;
		con.render(result,s2,"x","","; ");
		TEST_EQUAL(s2.str(),"idxtextxnumber; 5xblax45.11; 4711xbluffx471.123; ")
END_SECTION

START_SECTION((template<class StringListType> void executeQueries(const StringListType& queries)))
	  vector<String> qs;
	  qs.push_back("DROP TABLE IF EXISTS Dummy");
	  qs.push_back("CREATE TABLE Dummy (id int,text varchar(5))");
	  qs.push_back("INSERT INTO Dummy values (1,'bla'),(2,'bluff')");
	  qs.push_back("INSERT INTO Dummy values (3,'bla2'),(4,'bluff2')");
	  qs.push_back("DELETE FROM Dummy where id>2");

		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
	  con.executeQueries(qs);

	  QSqlQuery result = con.executeQuery("SELECT * FROM Dummy");
		stringstream s2;
		con.render(result,s2,"x","",";");
		TEST_EQUAL(s2.str(),"idxtext;1xbla;2xbluff;")
END_SECTION


START_SECTION((UInt getAutoId()))
		DBConnection con;
	  con.connect(db,user,password,host, port.toInt());
		con.executeQuery("DROP TABLE IF EXISTS Dummy");
	  con.executeQuery("CREATE TABLE `Dummy` (`id` INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY) TYPE = MYISAM ;");
	  con.executeQuery("INSERT INTO `Dummy` ( `id` ) VALUES ( NULL );");
	  TEST_EQUAL(con.getAutoId(),1)
	  con.executeQuery("INSERT INTO `Dummy` ( `id` ) VALUES ( NULL );");
	  TEST_EQUAL(con.getAutoId(),2)
END_SECTION

//remove Dummy table in the end
DBConnection con;
con.connect(db,user,password,host, port.toInt());
con.executeQuery("DROP TABLE IF EXISTS Dummy");
}
else
{
	ADD_MESSAGE("skipped")
}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



