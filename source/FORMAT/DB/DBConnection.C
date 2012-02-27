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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


//OpenMS includes
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

//QT
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlError>
#include <QtCore/QVariant>
#include <QtSql/QSqlRecord>

using namespace std;

namespace OpenMS
{

	DBConnection::DBConnection()
		: connection_name_()
	{
		
	}
	
	DBConnection::~DBConnection()
	{
		disconnect();
	}
	
	void DBConnection::connect(const String& db, const String& user, const String& password, const String& host,UInt port,const String& QTDBDriver, const String& connection_name)
	{
		connection_name_ = connection_name.toQString();

		QSqlDatabase::addDatabase(QTDBDriver.c_str(), connection_name_);
		
		getDB_().setHostName(host.c_str());
		getDB_().setUserName(user.c_str());
		getDB_().setDatabaseName(db.c_str());
		getDB_().setPassword(password.c_str());
		getDB_().setPort(port);
		if (!getDB_().open())
		{
			//construct query
			String query = String("Connecting to DB '") + db + "' ( host: '"+host+"', port: '"+String(port)+"', user: '"+user+"', password: '"+password+"')";
			//sore error
			String error = getDB_().lastError().databaseText();
			//close connection
			disconnect();
						
			throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,error);
		}
	}
	
	bool DBConnection::isConnected() const
	{
		return getDB_().isOpen();
	}
	
	void DBConnection::disconnect()
	{
		getDB_().close();
		QSqlDatabase::removeDatabase(connection_name_);
	}

	QSqlQuery DBConnection::executeQuery(const String& query, bool first)
	{
		QSqlDatabase db_handle = getDB_();
		
		if (!db_handle.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		QSqlQuery result(db_handle);
		
		//execute the query
		if (!result.exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle.lastError().text() );
		}
		if (first) result.first();
		
		return result;
	}
	
	UInt DBConnection::getId(const String& table, const String& column, const String& value)
	{
		String query = String("SELECT id FROM ") + table + " WHERE " + column + "='" + value + "' LIMIT 1";
		return executeQuery(query, true).value(0).toInt();
	}
	
	String DBConnection::DBName() const 
	{ 
		QSqlDatabase db_handle = getDB_();
		if (!db_handle.isOpen()) return "";
		return db_handle.databaseName();
	}
	
	void DBConnection::render(QSqlQuery& result, ostream& out, const String& separator, const String& line_begin, const String& line_end)
	{
		
	 	//if res is still empty, do nothing
	 	if (result.size()!=0)
	 	{
	 		Int col_count = result.record().count();
	 		
	 		out << line_begin;
	 		
	 		//render field names
	 		for (Int i = 0; i < col_count; i++) 
	 		{
	    	if (i!=0)
	    	{
	    		out << separator;
	    	} 
	
	      out << result.record().fieldName(i).toStdString();
	    } 
	    out << line_end;
			
			//set back internal pointer
			result.first();
			
			//render lines
	 		while(result.isValid())
	 		{
	 			out << line_begin;
		 		for (Int j = 0; j < col_count; j++) 
		 		{ 
		      if (j!=0)
	    		{
	    			out << separator;
	    		}
		      out << result.value(j).toString().toAscii().data();
		    } 			
	 			out << line_end;
	 			result.next();	
	 		}
	 	}
	}

	Int DBConnection::getIntValue(const String& table, const String& column, const String& id)
	{
		String query = String("SELECT ") + column + " FROM " + table + " WHERE id='" + id + "'";

		QSqlQuery result = executeQuery(query, true);

		if (!result.value(0).canConvert(QVariant::Int))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to int failed!");
		}
		return result.value(0).toInt();
	}

	double DBConnection::getDoubleValue(const String& table, const String& column, const String& id)
	{
		String query = String("SELECT ") + column + " FROM " + table + " WHERE id='" + id + "'";

		QSqlQuery result = executeQuery(query, true);

		if (!result.value(0).canConvert(QVariant::Double))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to double failed!");
		}
		return result.value(0).toDouble();
	}

	String DBConnection::getStringValue(const String& table, const String& column, const String& id)
	{
		String query = String("SELECT ") + column + " FROM " + table + " WHERE id='" + id + "'";

		QSqlQuery result = executeQuery(query, true);

		if (!result.value(0).canConvert(QVariant::String))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to double failed!");
		}
		return result.value(0).toString();
	}
	
	UInt DBConnection::getAutoId()
	{
		QSqlQuery result = executeQuery("SELECT LAST_INSERT_ID()", true);
		if (!result.value(0).canConvert(QVariant::Int))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to int failed in DBConnection::getAutoId()!");
		}
		return result.value(0).toInt();
	}
	
	
	// ---------------------------------------------------------------
	//
	// Excpetions
	//
	// ---------------------------------------------------------------
	
	DBConnection::InvalidQuery::InvalidQuery(const char* file, Int line, const char* function, const String& sql_query, const String& sql_error) 
		:	BaseException(file, line, function, "Invalid Query", "an SQL query failed")
	{
		what_ = String("Query '")+sql_query +"' failed: '"+sql_error+"'";
    Exception::GlobalExceptionHandler::getInstance().setMessage(what_);
	}
	
	DBConnection::InvalidQuery::~InvalidQuery() throw()
	{
	}
	
	DBConnection::NotConnected::NotConnected(const char* file, Int line, const char* function) 
		:	BaseException(file, line, function, "Not Connected", "the DBConnection was accessed but it is not connected to a SQL database")
	{
	}
	
	DBConnection::NotConnected::~NotConnected() throw()
	{
	}
	
}

