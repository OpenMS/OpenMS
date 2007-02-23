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
  	
    @note DB support is based on QT, and needs a qapplication. If you do not use a GUI create a background QApplication with:
    
    @code
      QApplication app(argc,argv,false);
    @endcode
    
    @note Do not use '*' in SELECT statments. The order of result columns is not definded then! Read the QT documentation for details.
    
    @todo make DBConnection a singleton as several concurrent DBConnections do not work! (Thomas S.)
    
    @ingroup DatabaseIO
  */
  class DBConnection
  {
    public:
      /**
      	@brief Exception in case of an invalid Query
      	
      	
      	
        @ingroup Exceptions
      */
      class InvalidQuery
            : public Exception::Base
      {
        public:
          InvalidQuery(const char* file, SignedInt line, const char*  function, std::string sql_query, std::string sql_error) throw();
          ~InvalidQuery() throw();
      };

      /**
       	@brief Exception in case of trying to execute a query without having established a database connection
      	
      	
      	
        @ingroup Exceptions
      */
      class NotConnected
            : public Exception::Base
      {
        public:
          NotConnected(const char* file, SignedInt line, const char*  function) throw();
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
	      @param QTDBDriver the QT database driver used for the connection (default: "QMYSQL3", unless you canged it in configure. See "OpenMS/include/OpenMS/config.h")
      */
      void connect(const std::string& db, const std::string& user, const std::string& password, const std::string& host = "localhost", UnsignedInt port=3306, const std::string& QTDBDriver = DB_PLUGIN ) throw(InvalidQuery);
			
			/// returns if a connection is established.
			bool isConnected() const;
			
      ///disconnects from the SQL database.
      void disconnect();

      /**
      	@brief Executes a query
      	
      	@note Make sure that the result was created after connecting to the DB! Otherwise the query fails!
      	
      	@param query the query itself
      	@param result the results are written to this object
      */
      void executeQuery(const std::string& query, QSqlQuery& result) throw(InvalidQuery, NotConnected);

			/**
				@brief Returns a single field of a table as an integer
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset
			*/
			SignedInt getIntValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError);

			/**
				@brief Returns a single field of a table as a double
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset
			*/
			double getDoubleValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError);

			/**
				@brief Returns a single field of a table as string
				
				The table has to contain an <tt>id</tt> column.
		
				@param table The table to look the field up
				@param column The column of the table
				@param id The id of the dataset
			*/
			String getStringValue(const std::string& table, const std::string& column, const std::string& id) throw (InvalidQuery,NotConnected,Exception::ConversionError);

			/**
				@brief Looks up the ID for a specific entry in an table 
				
				If several entries in the table have the desired value in the column, the first one is returned. 
		
				@param table The table to look the field up
				@param column The column of the table
				@param value The value the selected @p column has
			*/
			UnsignedInt getId(const std::string& table, const std::string& column, const std::string& value) throw (InvalidQuery,NotConnected);
			
			/// Returns the last auto_increment ID of the SQL database
			UnsignedInt getAutoId();
			
			/// Returns the name of the connected DB
      std::string DBName() const;
      
			/**
				@brief Dumps a query result in table format into a stream.
				
				To dump a result as HTML table, use render(result, cout,"&lt;/td&gt;&lt;td&gt;","&lt;tr&gt;&lt;td&gt;","&lt;/td&gt;&lt;/tr&gt;");
				
				@param result The result to render
				@param out The output stream to use
				@param separator The string between the fields
				@param line_begin The string at the beginning of each line
				@param line_end The string at the end of each line
			*/
			void render(QSqlQuery& result, std::ostream& out=std::cout , const std::string& separator=" | " , const std::string& line_begin="" , const std::string& line_end="\n");


			/**
				@brief Executes all SQL queries from an container.
				
				Each line has to be a query or empty.
				
				@param queries A STL-compliant container of OpenMS String objects 
			*/
			template <class StringListType>
			void executeQueries(const StringListType& queries) throw(InvalidQuery,NotConnected);

    private:
			
			/**
				@brief Executes internal queries.
				
				This method does not change the last query and last result 
			*/
			QSqlQuery& executeQuery_(const std::string& query) throw(InvalidQuery,NotConnected);
      
      /// The real database handle
      QSqlDatabase db_handle_;

      /**
      	@brief A pointer to the result of the last query.
      	
      	Used for internal subqueries.
      */
      QSqlQuery* lir_;

  };



	//---------------------------------------------------------------
	//  Implementation of the inline / template functions
	//---------------------------------------------------------------

	template <class StringListType>
	void DBConnection::executeQueries(const StringListType& queries) throw(InvalidQuery,NotConnected)
	{
		String line;
		for( typename StringListType::const_iterator it = queries.begin(); it != queries.end(); ++it)
		{
			line = *it;
			line.trim();
			if (line != "")
			{
				executeQuery_(line);
			}
		}
	}
}

#endif
