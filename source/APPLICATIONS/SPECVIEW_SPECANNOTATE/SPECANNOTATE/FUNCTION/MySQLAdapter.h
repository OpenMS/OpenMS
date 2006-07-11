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
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef MYSQLADAPTER_H
#define MYSQLADAPTER_H


//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>

#ifndef OPENMS_CONFIG_H
#include <OpenMS/config.h>
#endif

#include "../config_specannotate.h"
#include "../config_specannotate.h"

//std and STL includes
#include <string>
#include <iostream>
#include <vector>
#include <map>

//QT includes
#include <qsqldatabase.h>


namespace OpenMS
  {

  /** This class is just a very simple wrapper around the QT-class \c QSqlDatabase
   *  It makes it possible to have multiple instances of QSqlDatabase in one program.
   *  In this application this is very crucial, since both of its threads need database-access and cannot use a single 
   *  instance of QSqlDatabase together.
   */
  /*
  class SQLDatabase: public QSqlDatabase
  {
  public:
    SQLDatabase(const QString& type, const QString& name, QObject * parent=0, const char * objname=0 ) : 
      QSqlDatabase( type, name, parent, objname ) {};
      //SQLDatabase( QSqlDriver * driver, QObject * parent = 0, const char * objname = 0 ) :
      //QSqlDatabase ( driver, parent, objname ) {};
      //    ~SQLDatabase() {};
  };
  */


  /** 
		This class provides a generic interface for MySQL database access.
		It uses the \c QT library to handle the database
		OpenMS Exceptions are used.
		This class contains code from the OpenMS class \c OpenMS/FORMAT/DBAdapter.h
		(OpenMS/FORMAT/DBAdapter.h was simply "stripped" to provide a GENERIC database Interface)
  	
  	@todo replace by DBConnection 
  */
  class MySQLAdapter
    {

    public:

      //! Exception in case of an invalid Query
    class InvalidQuery : public OpenMS::Exception::Base
        {
        public:
          InvalidQuery(const char* file, int line, const char* function, std::string mysql_error) throw();
          ~InvalidQuery() throw();
        };

      //! Exception in case of trying to execute a query without having established a database connection
    class NotConnected : public OpenMS::Exception::Base
        {
        public:
          NotConnected(const char* file, int line, const char* function) throw();
          ~NotConnected() throw();
        };

      //! Exception in case of trying to get an unary Result which is not unary
    class NoUnaryResult : public OpenMS::Exception::Base
        {
        public:
          NoUnaryResult(const char* file, int line, const char* function) throw();
          ~NoUnaryResult() throw();
        };


      //! default constructor
      MySQLAdapter();

      //! copy constructor.<br> (note that the database connection is not copied. you have to connect the new object)
      MySQLAdapter(const MySQLAdapter& adapter);


      //! destructor
      virtual ~MySQLAdapter();

      //! assignment operator<br> (note that the database connection is resetted. you have to connect the object again.)
      MySQLAdapter& operator = (const MySQLAdapter& adapter);

      /**  connects the MySQLAdapter to a MySQL database.<br><br>
           @param user the login of the server
           @param password the password for the user
           @param host the host where the server is running (default: "localhost")
      */
      void connect(const char* user, const char* password, const char* host = "localhost",
                   const std::string& QTDBDriver = QTDATABASEDRIVER) throw(InvalidQuery);

      /** executes a query.<br>
          the sting is stored in last_query).<br>
          the query result can be accessed via last_result.<br><br>
          @param query the query std::string
      */
      void executeQuery(std::string query_string) throw(InvalidQuery, NotConnected);


      //! accessor for the last query result
      QSqlQuery& lastResult();

      //! accessor for the query std::string used in the last query
      const std::string& lastQuery();

      /** creates a new DB.<br><br>
          @param db the database name
      */
      void createDB(std::string db) throw(InvalidQuery,NotConnected);

      /** creates a new DB.<br><br>
          @param db the database name
      */
      void deleteDB(std::string db) throw(InvalidQuery,NotConnected);

      /** selects the DB for use DB.<br><br>
          @param db the database name
      */
      void selectDB(std::string db) throw(InvalidQuery,NotConnected);

      //! returns the last error message from the MySQL DB.<br><br>
      std::string error();

      /** executes all SQL queries in the a text file.<br>
          Each line has to be a query or empty (without whitespaces).<br>
          @param filename the file name
      */
      void executeScript(const char* filename) throw(InvalidQuery,NotConnected);

      /** dumps the result in table form into a stream.<br>
          if result equals 0, last_result is taken as input.<br><br>
          to dump a result to cout, simply use render(result);<br>
          to dump a result as HTML table, use render(result,cout,"&lt;/td>&lt;td>","&lt;tr>&lt;td>","&lt;/td>&lt;/tr>");
      @param result the query result to render
      @param out the output stream to use
      @param separator std::string between the fields
      @param line_begin std::string at the beginning of each line
      @param line_end sting at the end of each line
      */
      void render(std::ostream& out=std::cout , const std::string& separator=" | " , const std::string& line_begin="" , const std::string& line_end="\n");


      /** returns one (unary) string as result
       *   this function can be used, if a query should yield only one value, like e.g. in \c AminoAcid::getName()
       */
      std::string getUnaryResult() throw(NoUnaryResult);

      //!writes unary result into result and returns true, if there is one. if there is none, it returns false. if there are multiple: exception
      bool ifGetUnaryResult(std::string& result) throw(NoUnaryResult);

    private:

      /**
         the database handle
      */
      //SQLDatabase* db_handle_;
      QSqlDatabase* db_handle_;


      /**
         contains the std::string of the last query
      */
      std::string last_query_;

      /**
         a pointer to the result of the last query.<br>
         the result will be overwritten by the next query executed with executeQuery.
      */
      QSqlQuery* lr_; //former last_result_

      /**
         a pointer to the result of the last query.<br>
         the result will be overwritten by the next query executed with executeQuery_.
      */
      QSqlQuery* lir_; //former last_internal_result_

      /**
         executes a query
      */
      QSqlQuery& executeQuery_(const std::string& query_string) throw(InvalidQuery, NotConnected);

    };
}

#endif  //MYSQLADAPTER
