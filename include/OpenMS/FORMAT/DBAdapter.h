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

#ifndef OPENMS_FORMAT_DBADAPTER_H
#define OPENMS_FORMAT_DBADAPTER_H

//#include <OpenMS/config.h>

//OpenMS includes
#include <OpenMS/FORMAT/PersistenceManager.h>
#include <OpenMS/FORMAT/DBConnection.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

//QT includes
#include <qsqlquery.h>

//std and STL includes
#include <string>

namespace OpenMS
{	
  /** 
  	@brief A class for accessing and storing data in a SQL database
  	
    It can be used to create objects from the DB or store them in the DB.
    In order to write an object call operator<<(const PersistentObject&) of the PersistenceManager.
    
    @ingroup DatabaseIO
    
    @note This class will be reimplemented in the next version of OpenMS. We do not recommend using it!
    
    @todo speed up by implementing a not so generic version (Marc)
    @todo add setup and clear method for the DB (Marc)
    @todo do not create a new entry when the object has a DB id already (Marc)
    @todo support for metadata (Marc)
    @todo loadMSSpectrum(id, spectrum) loadMSExperiment(id,experiment) (Marc)
  */
  class DBAdapter
        : public PersistenceManager
  {
    public:

			/// Dimension order of a HPLC-MS experiment
			enum DimensionId 
			{ 
				MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ, 
				RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT 
			};

			/// Default constructor
      DBAdapter(DBConnection& db_con);

      /// Destructor
      ~DBAdapter();

			/**
				@brief Sets the UID to read from the DB.
				
			*/
			void readId(UID id);
			
			/// Reads a MSExperiment (for convenience)
			MSExperiment<>* loadMSExperiment(UID id);

			/// Reads a MSSpectrum (for convenience)
			MSSpectrum<>* loadSpectrum(UID id);
			
			/** 
				@name Layer 0 Methods of the PersistenceManager
			*/
			//@{
			//writing
			virtual void writeTrailer (const char* name);
			virtual void writeHeader (const char* signature, const char* name, const PersistentObject* object );
      virtual void writePrimitiveHeader (const char* signature, const char* name);
			virtual void writePrimitiveTrailer();
			virtual void put(const SignedInt value);
			virtual void put(const UnsignedInt value);
			virtual void put(const double value);
			virtual void put(const std::string& value);
			//reading
			virtual bool getObjectHeader(String& stream_name);
			virtual bool objectsToDeserialize();
			virtual bool checkPrimitiveHeader(const char* type_name, const char* name);
			virtual bool checkPrimitiveTrailer();
			virtual void get(double& d);
			virtual void get(UnsignedInt& i);
			virtual void get(SignedInt& i);
			virtual void get(string& s);
			virtual void get(UID& id);
			virtual bool checkObjectReferenceHeader(const char* type_name, const char* name);
			//both
			virtual void clear();
			//@}	
						
    private:
    	///reference to the DB connection handed over in the constructor
    	DBConnection& db_con_;

    	///Each object type that corresponds to a Table has to be registered in the constructor with this method
    	void registerType_(const std::string& stream_name, const std::string& table, UnsignedInt table_id, const std::string& fields);

 			///Mapping of signatures to a table name /table_id pair
 			HashMap<std::string , std::pair<std::string , UnsignedInt> > signature_to_tablenameid_;

 			///Mapping of table id to a signature
 			HashMap<UnsignedInt, std::string > tableid_to_signature_;

 			///Mapping of table id to a table name / field list pair
 			HashMap<UnsignedInt, std::pair < std::string , std::string > > tableid_to_tablenamefields_;
    	
			/** 
				@name Private members and functions for writing
			*/
			//@{		
      ///high id for id generation
      UnsignedInt high_id_;

      ///low id for id generation
      UnsignedInt low_id_;
    	
      /**
      	@brief Generates a unique DB ID
      	
      	@param id the id of table the new object is stored in 
      */
      UID generateID_(UnsignedInt id) throw(DBConnection::InvalidQuery, DBConnection::NotConnected);
  
  		///The query constructed for the current object
  		std::stringstream query_;
  		
  		///Mapping of UID to table name
 			HashMap<UID , std::string > uid_to_tablename_;
 			
 			/// Last ID that was assigned from the DB
 			UID last_id_;
			//@}

			/** 
				@name Private members and functions for reading
			*/
			//@{		
      /// UID to read from the DB
      std::list<UID> read_ids_;
 			
 			///The id of the field for the current primitive
 			SignedInt field_id_;
 			
 			///Pointer to the last query-result
 			QSqlQuery* last_query_; 
 			
			//@}
  };

}
#endif
