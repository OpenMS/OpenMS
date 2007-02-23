// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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


//OpenMS includes
#include <OpenMS/FORMAT/DB/DBConnection.h>

//QT
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlError>
#include <QtCore/QVariant>
#include <QtSql/QSqlRecord>

using namespace std;

namespace OpenMS
{

	DBConnection::DBConnection()
		:db_handle_(),
		lir_(0)
	{
		
	}
	
	DBConnection::~DBConnection()
	{
		if (db_handle_.isOpen())
		{
			disconnect();
		}
		QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
	}
	
	void DBConnection::connect(const string& db, const string& user, const string& password, const string& host,UnsignedInt port,const string& QTDBDriver) throw(InvalidQuery)
	{
		db_handle_ = QSqlDatabase::addDatabase(QTDBDriver.c_str(),QSqlDatabase::defaultConnection);
		db_handle_.setHostName(host.c_str());
		db_handle_.setUserName(user.c_str());
		db_handle_.setDatabaseName(db.c_str());
		db_handle_.setPassword(password.c_str());
		db_handle_.setPort(port);
		if (!db_handle_.open())
		{
			//construct query
			String query = "Connecting to DB ";
			query = query + db + "( host: '"+host+"' ,port: '"+String(port)+"' ,user: '"+user+"' ,password: '"+password+"')";
			//sore error
			string error = db_handle_.lastError().databaseText().toAscii().data();
			//close connection
			QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
			
			throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,error);
		}
		else
		{
			lir_ = new QSqlQuery();
		}
	}
	
	bool DBConnection::isConnected() const
	{
		return db_handle_.isOpen();
	}
	
	void DBConnection::disconnect()
	{
		db_handle_.close();
	}

	void DBConnection::executeQuery(const string& query, QSqlQuery& result) throw(InvalidQuery,NotConnected)
	{
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		//execute the query
		if (!result.exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
		result.first();
	}
	
	UnsignedInt DBConnection::getId(const std::string& table, const std::string& column, const std::string& value) throw(InvalidQuery,NotConnected)
	{
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		String query;
		query = query + "SELECT id FROM " + table + " WHERE " + column + "='" + value + "' LIMIT 1";
		
		//execute the query
		if (!lir_->exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
		
		lir_->first();
		return lir_->value(0).toInt();
	}
	
	std::string DBConnection::DBName() const 
	{ 
		if (!db_handle_.isOpen()) return "";
		return db_handle_.databaseName().toAscii().data();
	}

	QSqlQuery& DBConnection::executeQuery_(const string& query) throw(InvalidQuery,NotConnected)
	{
		//check if there is a connection active
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		if (!lir_->exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
	  return *lir_;
	}
	
		void DBConnection::render(QSqlQuery& result, ostream& out, const string& separator, const string& line_begin, const string& line_end)
	{
		
	 	//if res is still empty, do nothing
	 	if (result.size()!=0)
	 	{
	 		SignedInt col_count = result.record().count();
	 		
	 		out << line_begin;
	 		
	 		//render field names
	 		for (SignedInt i = 0; i < col_count; i++) 
	 		{
	    	if (i!=0)
	    	{
	    		out << separator;
	    	} 
	
	      out << result.record().fieldName(i).toAscii().data() ;
	    } 
	    out << line_end;
			
			//set back internal pointer
			result.first();
			
			//render lines
	 		while(result.isValid())
	 		{
	 			out << line_begin;
		 		for (SignedInt j = 0; j < col_count; j++) 
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

	SignedInt DBConnection::getIntValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError)
	{
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		String query;
		query = query + "SELECT " + column + " FROM " + table + " WHERE id='" + id + "'";
		
		//execute the query
		if (!lir_->exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
		
		lir_->first();
		if (!lir_->value(0).canConvert(QVariant::Int))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to int failed!");
		}
		return lir_->value(0).toInt();
	}

	double DBConnection::getDoubleValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError)
	{
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		String query;
		query = query + "SELECT " + column + " FROM " + table + " WHERE id='" + id + "'";
		
		//execute the query
		if (!lir_->exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
		
		lir_->first();
		if (!lir_->value(0).canConvert(QVariant::Double))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to double failed!");
		}
		return lir_->value(0).toDouble();
	}

	String DBConnection::getStringValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError)
	{
		if (!db_handle_.isOpen())
		{
			throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		String query;
		query = query + "SELECT " + column + " FROM " + table + " WHERE id='" + id + "'";
		
		//execute the query
		if (!lir_->exec(query.c_str()))
		{
	    throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,query,db_handle_.lastError().text().toAscii().data() );
		}
		
		lir_->first();
		if (!lir_->value(0).canConvert(QVariant::String))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to double failed!");
		}
		return lir_->value(0).toString().toAscii().data();
	}
	
	UnsignedInt DBConnection::getAutoId()
	{
		executeQuery_("SELECT LAST_INSERT_ID()");
		lir_->first();
		if (!lir_->value(0).canConvert(QVariant::Int))
		{
			throw Exception::ConversionError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Conversion of QVariant to int failed in DBConnection::getAutoId()!");
		}
		return lir_->value(0).toInt();
	}
	
	
	// ---------------------------------------------------------------
	//
	// Excpetions
	//
	// ---------------------------------------------------------------
	
	DBConnection::InvalidQuery::InvalidQuery(const char* file, SignedInt line, const char* function, string sql_query, string sql_error) throw()
		:	Base(file, line, function, "Invalid Query", "an SQL query failed")
	{
		what_ = String("Query '")+sql_query +"' failed: '"+sql_error+"'";
		Exception::globalHandler.setMessage(what_);
	}
	
	DBConnection::InvalidQuery::~InvalidQuery() throw()
	{
	}
	
	DBConnection::NotConnected::NotConnected(const char* file, SignedInt line, const char* function) throw()
		:	Base(file, line, function, "Not Connected", "the DBConnection was accessed but it is not connected to a SQL database")
	{
	}
	
	DBConnection::NotConnected::~NotConnected() throw()
	{
	}
	
}

