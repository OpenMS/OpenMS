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

#ifndef OPENMS_FORMAT_DB_DBCONNECTION_H
#define OPENMS_FORMAT_DB_DBCONNECTION_H

#include <OpenMS/config.h>

//QT includes
#include <QtSql/QSqlDatabase>

//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>

#include <iostream>

namespace OpenMS
{
  /** 
  	@brief A class for connecting to a SQL database
    
    @note Do not use '*' in SELECT statments. The order of result columns is not definded then! Read the QT documentation for details.
   	
    @ingroup DatabaseIO
  */
  class OPENMS_DLLAPI DBConnection
  {
    public:
      /**
      	@brief Exception in case of an invalid Query
      	
      	
      	
        @ingroup Exceptions
      */
      class OPENMS_DLLAPI InvalidQuery
            : public Exception::BaseException
      {
        public:
          InvalidQuery(const char* file, Int line, const char*  function, const String& sql_query, const String& sql_error) ;
          ~InvalidQuery() throw();
      };

      /**
       	@brief Exception in case of trying to execute a query without having established a database connection
      	
      	
      	
        @ingroup Exceptions
      */
      class OPENMS_DLLAPI NotConnected
            : public Exception::BaseException
      {
        public:
          NotConnected(const char* file, Int line, const char*  function) ;
          ~NotConnected() throw();
      };
      
      
      
      /// Default constructor
      DBConnection();

      ///Destructor
     	~DBConnection();
			
      /**
      	@brief Connects to a SQL database.
      	
	      @param db the database name
	      @param user the login of the server
	      @param password the password for the user
	      @param host the host where the server is running (default: "localhost")
	      @param port the port where the server is listening (default: 3306)
	      @param QTDBDriver the QT database driver used for the connection (default: "QMYSQL", unless you canged it in configure. See "OpenMS/include/OpenMS/config.h")
	      @param connection_name Name of the connection (needed for several concurrent connections only)
      
      	@exception InvalidQuery is thrown if the database connection could not be opened
      */
      void connect(const String& db, const String& user, const String& password, const String& host = "localhost", UInt port=3306, const String& QTDBDriver = DB_PLUGIN, const String& connection_name="OpenMS_default_connection");
			
			/// returns if a connection is established.
			bool isConnected() const;
			
      /**
				@brief disconnects from the SQL database

				@note In order to disconnect, all queries running on a database must be closed. See the QSqlDatabase::removeDatabase(...) documentation.
			*/
      void disconnect();

      /**
      	@brief Executes a query and returns the result
      	
      	The internal pointer of the returned result is positioned @b before the first row.
      	
      	@param query an SQL query
      	@param first if @em true, the internal pointer of the returned result is positioned to the first record.

				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
				
      */
      QSqlQuery executeQuery(const String& query, bool first = false);

			/**
				@brief Returns a single field of a table as an integer
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset
				
				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
				@exception Exception::ConversionError is thrown if the value could not be converted to the requested type
			*/
			Int getIntValue(const String& table, const String& column, const String& id);

			/**
				@brief Returns a single field of a table as a double
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset

				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
				@exception Exception::ConversionError is thrown if the value could not be converted to the requested type

			*/
			double getDoubleValue(const String& table, const String& column, const String& id);

			/**
				@brief Returns a single field of a table as string
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset

				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
				@exception Exception::ConversionError is thrown if the value could not be converted to the requested type
			*/
			String getStringValue(const String& table, const String& column, const String& id);

			/**
				@brief Looks up the ID for a specific entry in an table 
				
				If several entries in the table have the desired value in the column, the first one is returned. 
		
				@param table The table to look the field up
				@param column The column of the table
				@param value The value the selected @p column has

				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
			*/
			UInt getId(const String& table, const String& column, const String& value);
			
			/// Returns the last auto_increment ID of the SQL database
			UInt getAutoId();
			
			/// Returns the name of the connected DB
      String DBName() const;
      
			/**
				@brief Dumps a query result in table format into a stream.
				
				To dump a result as HTML table, use render(result, cout,"&lt;/td&gt;&lt;td&gt;","&lt;tr&gt;&lt;td&gt;","&lt;/td&gt;&lt;/tr&gt;");
				
				@param result The result to render
				@param out The output stream to use
				@param separator The string between the fields
				@param line_begin The string at the beginning of each line
				@param line_end The string at the end of each line
			*/
			void render(QSqlQuery& result, std::ostream& out=std::cout , const String& separator=" | " , const String& line_begin="" , const String& line_end="\n");


			/**
				@brief Executes all SQL queries from an container.
				
				Each line has to be a query or empty.
				
				@param queries A STL-compliant container of OpenMS String objects 

				@exception InvalidQuery is thrown if an invalid SQL query was given
				@exception NotConnected if there is no database connection
			*/
			template <class StringListType>
			void executeQueries(const StringListType& queries);

    private:
			
      /// Name (handle) of the connection
      QString connection_name_;
			
			/// Retruns the current database connection defined by connection_name_
			inline QSqlDatabase getDB_() const
			{
				return QSqlDatabase::database(connection_name_,false);
			}
			
  };



	//---------------------------------------------------------------
	//  Implementation of the inline / template functions
	//---------------------------------------------------------------

	template <class StringListType>
	void DBConnection::executeQueries(const StringListType& queries)
	{
		String line;
		for( typename StringListType::const_iterator it = queries.begin(); it != queries.end(); ++it)
		{
			line = *it;
			line.trim();
			if (line != "")
			{
				executeQuery(line);
			}
		}
	}
}

#endif
