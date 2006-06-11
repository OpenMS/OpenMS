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
// $Id: MySQLAdapter.C,v 1.3 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------



// std + STL includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>

//QT
#include <qsqlrecord.h>


#include "MySQLAdapter.h"


using namespace OpenMS;
using namespace std;


string MySQLAdapter::getUnaryResult() throw(NoUnaryResult)
{
  if ((lr_->isActive()) && (lr_->size() == 1) && (lr_->first()))
    {
      if (lr_->value(0).isNull())
	{
	  return "0";
	}
      else
	{
	  std::ostringstream st;
	  st << (lr_->value(0).toString());
	  return std::string(st.str());
	}
    }
  else
    {
      cerr << "Throwing exception NoUnaryResult because of result: " << endl;
      render();
      cerr << "As result of query:" << endl;
      cerr << last_query_;

      throw NoUnaryResult(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
}


bool MySQLAdapter::ifGetUnaryResult(string& result) throw(NoUnaryResult)
{
  if ((!lr_->isActive()) || (lr_->size() > 1))
    {
      cerr << "Throwing exception NoUnaryResult because of result: " << endl;
      render();
      cerr << "As result of query:" << endl;
      cerr << last_query_;

      throw NoUnaryResult(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
  else if ((lr_->isActive()) && (lr_->size() < 1))
    {
      return false;
    }
  else if ((lr_->isActive()) && (lr_->size() == 1) && (lr_->first()))
    {
      if (lr_->value(0).isNull())
	{
	  result = "0";
	}
      else
	{
	  std::ostringstream st;
	  st << (lr_->value(0).toString());
	  result = std::string(st.str());
	}
      return true;
    }
  return false;
}


void MySQLAdapter::render(ostream& out , const string& separator , const string& line_begin , const string& line_end)
{
  //if res is still empty, do nothing
  if (lr_->size()!=0)
    {
      SignedInt col_count = db_handle_->record(*lr_).count();

      //render field names
      for (SignedInt i = 0; i < col_count; i++)
        {
          if (i!=0)
            {
              out << separator;
            }

          out  << line_begin << db_handle_->record(*lr_).fieldName(i) ;
        }
      out << line_end;

      //render lines
      while(lr_->next())
        {
          for (SignedInt j = 0; j < col_count; j++)
            {
              if (j!=0)
                {
                  out << separator;
                }
              out << line_begin << lr_->value(0).toString();
            }
          out << line_end;
        }
    }
}


void MySQLAdapter::connect(const char* user, const char* password, const char* host, const string& QTDBDriver) throw(InvalidQuery)
{
  if (!(QSqlDatabase::contains("db_handle_")))
    {
      db_handle_ = QSqlDatabase::addDatabase(QTDBDriver, "db_handle_");
    }
  else
    {
      db_handle_ = QSqlDatabase::database("db_handle_");
    }

  //db_handle_ = new SQLDatabase(QTDBDriver, "db_handle_");

  db_handle_->setHostName(host);
  db_handle_->setUserName(user);
  db_handle_->setDatabaseName(DATABASE);
  db_handle_->setPassword(password);
  if (!db_handle_->open())
    {
      string tmp = db_handle_->lastError().databaseText().ascii();
      QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
      db_handle_=0;
      throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,tmp);
    }
  else
    {
      lir_ = new QSqlQuery(db_handle_);
      lr_ = new QSqlQuery(db_handle_);
    }
}


void MySQLAdapter::createDB(string db) throw(InvalidQuery,NotConnected)
{
  executeQuery_("CREATE DATABASE "+db);
}


void MySQLAdapter::deleteDB(string db) throw(InvalidQuery,NotConnected)
{
  executeQuery_("DROP DATABASE "+db);
}


void MySQLAdapter::selectDB(string db) throw(InvalidQuery,NotConnected)
{
  executeQuery_("USE "+db);
}


string MySQLAdapter::error()
{
  return db_handle_->lastError().text();
}


void MySQLAdapter::executeScript(const char* filename) throw(InvalidQuery,NotConnected)
{
  ifstream infile;
  char line[10000];
  infile.open (filename);
  while (infile.good())
    {
      infile.getline (line,10000, ';');
      if (strlen(line)!=0)
        executeQuery_(line);
    }
  infile.close();
}


QSqlQuery& MySQLAdapter::executeQuery_(const string& query_string) throw(InvalidQuery, NotConnected)
{
  //check if there is a connection active
  if (db_handle_==0)
    {
      throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

  db_handle_ = QSqlDatabase::database("db_handle_");
  delete lir_;
  lir_ = new QSqlQuery(db_handle_);

  if (!lir_->exec(query_string))
    {
      cerr <<"Invalid Query: "<< query_string << endl;
      throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,lir_->lastError().text() );
      //throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,db_handle_->lastError().text() );
    }
  return *lir_;
}


void MySQLAdapter::executeQuery(string query_string) throw(InvalidQuery,NotConnected)
{
  //set last query string
  last_query_=query_string;

  if (db_handle_==0)
    {
      throw NotConnected(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

  //new Instance of QSql
  db_handle_ = QSqlDatabase::database("db_handle_");
  delete lr_;
  lr_ = new QSqlQuery(db_handle_);

  //execute the query
  if (!lr_->exec(query_string))
    {
      cerr << "Invalid Query: " << query_string << endl;
      throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,lr_->lastError().text() );
      //throw InvalidQuery(__FILE__, __LINE__, __PRETTY_FUNCTION__,db_handle_->lastError().text() );
    }
}


MySQLAdapter::MySQLAdapter() : db_handle_(0),last_query_(""),lr_(0)
{}


MySQLAdapter::~MySQLAdapter()
{
  /*
  while(QSqlDatabase::contains(QSqlDatabase::defaultConnection))
    {
      QSqlDatabase::removeDatabase(QSqlDatabase::defaultConnection);
    }
    
  while(QSqlDatabase::contains("db_handle_"))
    {
      QSqlDatabase::removeDatabase("db_handle_");
    }
  */
  //delete db_handle_;

  //other stuff deleted by QT!
}


QSqlQuery& MySQLAdapter::lastResult()
{
  return *lr_;
}


const string& MySQLAdapter::lastQuery()
{
  return last_query_;
}


MySQLAdapter::MySQLAdapter(const MySQLAdapter& rhs) : db_handle_(0), last_query_(rhs.last_query_)
{
  lr_ = new QSqlQuery(*(rhs.lr_));
  lir_ = new QSqlQuery(*(rhs.lir_));
}


MySQLAdapter& MySQLAdapter::operator = (const MySQLAdapter& rhs)
{
  //check for self assignment!
  if (this==&rhs)
    {
      return *this;
    }

  //copy
  last_query_ = rhs.last_query_;
  lr_ = new QSqlQuery(*(rhs.lr_));
  lir_ = new QSqlQuery(*(rhs.lir_));

  //return reference to this
  return *this;
}




MySQLAdapter::InvalidQuery::InvalidQuery(const char* file, int line, const char* function, string mysql_error)
throw()
    :	Base(file, line, function, "Invalid Query", "a MySQL query failed")
{
  what_ = mysql_error.c_str();
  OpenMS::Exception::globalHandler.setMessage(what_);
}


MySQLAdapter::InvalidQuery::~InvalidQuery() throw()
{}


MySQLAdapter::NotConnected::NotConnected(const char* file, int line, const char* function)
throw()
    :	Base(file, line, function, "Not Connected", "the Adapter is not connected to a MySQL database")
{}


MySQLAdapter::NotConnected::~NotConnected() throw()
{}


MySQLAdapter::NoUnaryResult::NoUnaryResult(const char* file, int line, const char* function)
throw()
    :	Base(file, line, function, "NoUnaryResult", "the result of previous query is not unique")
{}


MySQLAdapter::NoUnaryResult::~NoUnaryResult() throw()
{}

