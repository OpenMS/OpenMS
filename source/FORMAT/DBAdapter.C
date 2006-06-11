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
// $Id: DBAdapter.C,v 1.32 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


//OpenMS includes
#include <OpenMS/FORMAT/DBAdapter.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/COMPARISON/CLUSTERING/Cluster.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFactory.h>

//STL
#include <typeinfo>

using namespace std;

namespace OpenMS
{
	
	DBAdapter::DBAdapter(DBConnection& db_con)
		: PersistenceManager(),
			db_con_(db_con),
			signature_to_tablenameid_(),
 			tableid_to_signature_(),
			tableid_to_tablenamefields_(),
			high_id_(0),
			low_id_(0),
			query_(),
			last_id_(0),
			read_ids_()
	{
		//register all classes that correspond to a Table: streamName -> (Table name / select statement)
		//The Object have to read the primitives in the order defined in the select statement
		
		//MSExperiment<>
		registerType_(RTTI::getStreamName< MSExperiment<> >() , "MSExperiment", 1, "id");
		//Spectrum<>
		registerType_(RTTI::getStreamName< MSSpectrum<> >() , "Spectrum", 2, "id, MS_Level, Retention_Start , Retention_Stop , Retention_Time");	
		//DPeak<1> 
		registerType_(RTTI::getStreamName< DPeak<1> >() , "Peak", 3, "id, Intensity, mz");
		//DPickedPeak<1>
		registerType_(RTTI::getStreamName< DPickedPeak<1> >() , "PrecursorInfo", 4, "id, Intensity, Charge, mz");	
	}
	
	DBAdapter::~DBAdapter()
	{
		
	}

  void DBAdapter::registerType_(const string& stream_name, const string& table, UnsignedInt table_id, const string& fields)
  {
  	signature_to_tablenameid_.insert( make_pair( stream_name , make_pair( table, table_id ) ) );
  	tableid_to_signature_.insert( make_pair( table_id, stream_name ) );
  	tableid_to_tablenamefields_.insert( make_pair( table_id, make_pair( table, fields) ) );
  }

	void DBAdapter::readId(UID id)
	{
		read_ids_.clear();
		read_ids_.push_back(id);
	}
	
	UID DBAdapter::generateID_(UnsignedInt table_id) throw(DBConnection::InvalidQuery,DBConnection::NotConnected)
	{	
		stringstream query;
		UnsignedInt new_id = 0;
		
		//fetch a new high id if no object was created so far 
		//or if the all low id's for a high id are used up
		if (low_id_==0)
		{
			db_con_.executeQuery("INSERT into _HighID values('')");
			db_con_.executeQuery("SELECT LAST_INSERT_ID()");
			last_query_ = &(db_con_.lastResult());
			last_query_->first();
			high_id_=last_query_->value(0).toInt();
			//cout << "New high ID: " << high_id_ << " ####################### "<<endl;
			high_id_ = high_id_ << 8;
		}
	
		// compose high/low id
		new_id=high_id_ + low_id_;
		// includent  low id
		low_id_++;
		// set to 0 if 256
		low_id_%=256;
		
		//store the id and table id in the mapping table
		query << "INSERT into _ID_Table values('" << new_id << "','" << table_id << "')";
	  db_con_.executeQuery(query.str());
		return new_id;
	}

	void DBAdapter::writeHeader (const char* signature, const char* /*name*/, const PersistentObject* object )
	{
		//cout << "writeHeader: " << signature << ""  << object<< endl;
		HashMap<string , pair<string , UnsignedInt> >::iterator it(signature_to_tablenameid_.find(signature));
		//new table
		if (it!=signature_to_tablenameid_.end())
		{
			//cout << "NEW QUERY "<< endl;
			//begin the new query
			const_cast<PersistentObject*>(object)->setPersistenceId(generateID_(it->second.second));
			last_id_ = object->getPersistenceId();
			query_.str("");
			query_ << "INSERT INTO " << it->second.first << " SET id='"<< object->getPersistenceId()<<"' ";
			uid_to_tablename_.insert( make_pair( object->getPersistenceId() , it->second.first ) );
		}
	}

	void DBAdapter::writeTrailer (const char* /*name*/)
	{
		if (current_->getPersistenceId()==last_id_)
		{
			// -- add foreign key if needed --
			const PersistentObject* parent(parents_[current_]);
			//look for the real parent, has not to be the direct parent object (as not all objects correspond to a table)
			//only those with a valid persistence ID are real parents
			while (true)
			{
				//no parent -> no foreign key
				if (parent ==0)
				{
					break;
				}
				// found the real parent	
				if (parent->getPersistenceId()!=0 && parent->getPersistenceId()!=last_id_)
				{
					//cout << "Parent: "<<parent<<" id: "<<parent->getPersistenceId() <<" table: "<<uid_to_tablename_[parent->getPersistenceId()]<<endl;
					//look up the table name of the parent
					query_ << ", " << uid_to_tablename_[parent->getPersistenceId()] << "_id= '" << parent->getPersistenceId() << "' ";
					break;
				};
				
				parent = parents_[parent];
			}
	
			//cout <<"QUERY: "<< query_.str() << endl;		
			db_con_.executeQuery(query_.str());
		}
	}

  void DBAdapter::writePrimitiveHeader (const char* /*signature*/, const char* name)
	{
		query_ << ", " << name << "=";
	}
	
	void DBAdapter::writePrimitiveTrailer()
	{
		
	}
				
	void DBAdapter::put(const SignedInt value)
	{
		query_ << "'" << value << "'";
	}
	
	void DBAdapter::put(const UnsignedInt value)
	{
		query_ << "'" << value << "'";
	}
	
	void DBAdapter::put(const double value)
	{
		query_ << "'" << value << "'";
	}

	void DBAdapter::put(const string& value)
	{
		query_ << "'" << value << "'";
	}

	void DBAdapter::clear()
	{
		//clean up
		query_.str("");
		uid_to_tablename_.clear();
		last_id_ = 0;
	}

	bool DBAdapter::getObjectHeader(String& stream_name)
	{
		
		stringstream ss, add;
		UID object_id = *(read_ids_.begin());
		
		//cout << "getObjectHeader " << object_id << " -- objects left: "<< (read_ids_.size()-1) << endl;
		
		// look up the table_id
		ss << "SELECT _Tables_id from _ID_Table WHERE id='" << object_id << "'";
		db_con_.executeQuery(ss.str());
		
		// no such entry in the ID_Tables_ table.
		if (db_con_.lastResult().size()==0)
		{
			return false;
		}

		db_con_.lastResult().first();		
		UnsignedInt table_id(db_con_.lastResult().value(0).toInt());
		
		//return the stream name
		stream_name = tableid_to_signature_[table_id];
		
		//load the data
		ss.str("");
		pair < string , string >* tmp_ptr(&(tableid_to_tablenamefields_[table_id]));
		
		//Object dependant operation (size of member arrays and stuff like that)
		switch (table_id)
		{
			case 1: //MSExperiment
				ss << "SELECT id FROM Spectrum WHERE MSExperiment_id='" << object_id << "'";
				db_con_.executeQuery(ss.str());
				ss.str("");
				last_query_ = &(db_con_.lastResult());
				while (last_query_->next())
				{
					read_ids_.push_back(last_query_->value(0).toInt());
					//cout << "-- add SPECTRUM: "<< last_query_->value(0).toInt() << endl;
				}
				
				// add spectrum count to the query
				add << ", " << last_query_->size();
				break;

			case 2: //Spectrum
				ss << "SELECT id FROM PrecursorInfo WHERE Spectrum_id='" << object_id << "'";
				db_con_.executeQuery(ss.str());
				db_con_.lastResult().first();	
				read_ids_.push_back(db_con_.lastResult().value(0).toInt());
				//cout << "-- add PrecursorInfo: "<< db_con_.lastResult().value(0).toInt() << endl;
				ss.str("");
				ss << "SELECT id FROM Peak WHERE Spectrum_id='" << object_id << "'";
				db_con_.executeQuery(ss.str());
				ss.str("");
				last_query_ = &(db_con_.lastResult());
				while (last_query_->next())
				{
					read_ids_.push_back(last_query_->value(0).toInt());
					//cout << "-- add PEAK: "<< last_query_->value(0).toInt() << endl;
				}
				
				// add spectrum count to the query
				add << ", " << last_query_->size();
				break;
		}
		
		ss << "SELECT "<< tmp_ptr->second << add.str() <<" FROM "<< tmp_ptr->first << " WHERE id='" << object_id << "'";	
		//cout << ss.str() << endl;
		db_con_.executeQuery(ss.str());
		last_query_ = &(db_con_.lastResult());
		last_query_->first();
		// set the field counter to -1 => first will be 0
		field_id_=-1; 
		
		//remove the read element
		read_ids_.pop_front();
		
		return true;
	}

	bool DBAdapter::objectsToDeserialize()
	{
		return (!read_ids_.empty());
	}

	bool DBAdapter::checkPrimitiveHeader(const char* /*stream_name*/, const char* /*name*/)
	{
		//cout << "PRIMITIVE HEDER: "<<name<<endl;
		return true;
	}
	
	bool DBAdapter::checkPrimitiveTrailer()
	{
		return true;
	}
	
	void DBAdapter::get(double& d)
	{
		d = last_query_->value(++field_id_).toDouble();
		//cout << "PRIMITIVE: double " << last_query_->value(field_id_).toDouble() <<endl;
	}
	
	void DBAdapter::get(UnsignedInt& i)
	{
		i = last_query_->value(++field_id_).toInt();
		//cout << "PRIMITIVE: UnsignedInt " << last_query_->value(field_id_).toInt()<<endl;
	}
	
	void DBAdapter::get(SignedInt& i)
	{
		i = last_query_->value(++field_id_).toInt();
		//cout << "PRIMITIVE: SignedInt " << last_query_->value(field_id_).toInt() << endl;
	}
	
	void DBAdapter::get(string& s)
	{
		s = last_query_->value(++field_id_).toString().ascii();
		//cout << "PRIMITIVE: string " << last_query_->value(field_id_).toString().ascii() << endl;
	}
	
	void DBAdapter::get(UID& id)
	{
		id = last_query_->value(++field_id_).toInt();
		//cout << "PRIMITIVE: UID " << last_query_->value(field_id_).toInt() <<endl;
	}

	bool DBAdapter::checkObjectReferenceHeader(const char* /*type_name*/, const char* /*name*/)
	{
		//cout << "checkObjectReferenceHeader: " << name << " " << type_name << endl;
		return true;
	}
	
	MSExperiment<>* DBAdapter::loadMSExperiment(UID id)
	{
		readId(id);
		PersistentObject* ptr;
		*this >> ptr;
		return dynamic_cast< MSExperiment<>* >(ptr);
	}

	MSSpectrum<>* DBAdapter::loadSpectrum(UID id)
	{
		readId(id);
		PersistentObject* ptr;
		*this >> ptr;
		return dynamic_cast< MSSpectrum<>* >(ptr);
	}
	
} //namespace


