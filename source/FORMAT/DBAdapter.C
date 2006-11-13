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


//OpenMS includes
#include <OpenMS/FORMAT/DBAdapter.h>

using namespace std;

namespace OpenMS
{
	
	DBAdapter::DBAdapter(DBConnection& db_con)
		: db_con_(db_con)
	{

	}
	
	DBAdapter::~DBAdapter()
	{
		
	}
	
	void DBAdapter::loadMetaInfo_(UID id, MetaInfoInterface& info)
	{
		//cout << "Loading Meta: " << id << endl;
		if (id==0)
		{
			info.clearMetaInfo();
			return;
		}
		
		QSqlQuery result;
		stringstream query;
		
		query << "SELECT Type-1,Name,Value FROM META_TypeNameValue WHERE fid_MetaInfo='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		while(result.isValid())
		{
			switch(result.value(0).asInt())
			{
				case 0: //string
					info.setMetaValue(result.value(1).asString().ascii(),result.value(2).asString().ascii());
					break;
				case 1: //double
					info.setMetaValue(result.value(1).asString().ascii(),result.value(2).asDouble());
					break;
				case 2: //int
					info.setMetaValue(result.value(1).asString().ascii(),result.value(2).asInt());
					break;
				default:
					throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"DBAdapter","Unknown META_TypeNameValue:type in DBAdapter!");
			}
			result.next();
		}
	}
	
	UID DBAdapter::storeMetaInfo_(const String& parent_table, UID parent_id, const MetaInfoInterface& info)
	{
		/* All MetaInfo data in META_TypeNameValue needs to be connected to entries in parent tables in an 1:n fashion.
		   The additional table META_MetaInfo forces links between META_TypeNameValue and the parent table to have
		     a globally unique ID.
		 	These links logically consist of the following structures:
		  - the parent table has a field called fid_MetaInfo (foreign key to id-field in table META_MetaInfo)
		  - table MetaInfo links parent_table in a 1:n-fashion to META_TypeNameValue:
		  	parent_table.fid_MetaInfo       -> META_MetaInfo.id
		  	META_TypeNameValue.fid_MetaInfo -> META_MetaInfo.id
		  
		  This is implemented in the following way:
		  - saving meta information needs the parent table's name and the id of the associated entry of the parent table
		  - As soon as meta information is to be saved:
		  	- the MetaInfo-ID is retrieved (if already existing)
		  	- all existing TypeNameValue entries are deleted
		  	- references (MetaInfo.id and parent_table.fid_MetaInfo) are deleted or created (if necessary)
		  	- new TypeNameValue entries are created 
		*/

		//init
		QSqlQuery result;
		stringstream query;
		
		query << "SELECT fid_MetaInfo FROM " << parent_table << " WHERE id='" << parent_id << "' AND fid_MetaInfo IS NOT NULL";
		db_con_.executeQuery(query.str(),result);
		
		
		//metainfo present => delete values of it
		UID meta_id = 0;
		if (result.size()==1)
		{
			// information exists in DB -> clearing content to replace with more recent information
			// cout << "MetaInfo for entry '" << parent_id << "' in table '" << parent_table << "' already exists => deleting TypeNameValues..." << endl;
			meta_id = result.value(0).asInt();
			query.str("");
			query << "DELETE FROM META_TypeNameValue WHERE fid_MetaInfo='" << meta_id << "'";
			db_con_.executeQuery(query.str(),result);
		}
		
		//connection between metainfo and object
		if (info.isMetaEmpty() && meta_id!=0)
		{
			// information exists in DB, but new information is empty -> reference to old information needs to be deleted
			// cout << "Nothing to save for entry '" << parent_id << "' in table '" << parent_table << "' => clearing reference..." << endl;
			query.str("");
			query << "DELETE FROM META_MetaInfo WHERE id='" << meta_id << "'";
			db_con_.executeQuery(query.str(),result);
			query.str("");
			query << "UPDATE " << parent_table << " SET fid_MetaInfo=NULL WHERE id='" << parent_id << "'";
			db_con_.executeQuery(query.str(),result);
		}
		
		if ((!info.isMetaEmpty()) && meta_id==0)
		{
			// information does not exist in DB, but there is new information to save -> create reference
			// cout << "Unsaved MetaInfo to save, saving to entry '" << parent_id << "' in table '" << parent_table << "' => creating reference..." << endl;
			db_con_.executeQuery("INSERT INTO META_MetaInfo () VALUES ()",result);
			meta_id = db_con_.getAutoId();
			query.str("");
			query << "UPDATE " << parent_table << " SET fid_MetaInfo='" << meta_id << "' WHERE id='" << parent_id << "'";
			db_con_.executeQuery(query.str(),result);
		}

		if (!info.isMetaEmpty())
		{
			// there is new information to save -> create contents
			// cout << "writing meta values..." << endl;
			
			query.str("");
			query << "INSERT INTO META_TypeNameValue (fid_MetaInfo,Name,Type,Value) VALUES ";
			vector<string> keys;
			info.getKeys(keys);
			const DataValue* val = 0;
			vector<string>::iterator it = keys.begin();
			while (true)
			{
				query << "('" << meta_id << "','" << *it;
				val = &info.getMetaValue(*it);
				// cout << *it << "=" << *val << endl;
				switch (val->valueType())
				{
					case DataValue::STRVALUE:
						query << "','string','" << string(*val);
						break;
					case DataValue::INTVALUE:
						query << "','int','" << int(*val);
						break;
					case DataValue::SHOVALUE:
						query << "','int','" << short(*val);
						break;
					case DataValue::LONVALUE:
						query << "','int','" << long(*val);
						break;
					case DataValue::DOUVALUE:
						query << "','double','" << double(*val);
						break;
					case DataValue::FLOVALUE:
						query << "','double','" << float(*val) ;
						break;
					default:
						throw Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__,"DBAdapter","Unknown DataValue type in DBAdapter!");
				}
				query << "')";
				++it;
				if (it==keys.end())
				{
					break;
				}
				query << ",";
			}
			db_con_.executeQuery(query.str(),result);
		}
		
		return meta_id;
	}
} //namespace


