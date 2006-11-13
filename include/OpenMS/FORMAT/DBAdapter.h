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
#include <OpenMS/FORMAT/DBConnection.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

//QT includes
#include <qsqlquery.h>
#include <qdatetime.h>

//std and STL includes
#include <string>

namespace OpenMS
{	
  /** 
  	@brief A class for accessing and storing data in a SQL database
  	
    It can be used to create objects from the DB or store them in the DB.
    
    @ingroup DatabaseIO
    
    @todo add setup and clear method for the DB (Marc)
    @todo test: update function (Marc)
    @todo add support for missing classes / add missing members to classes(Marc)
  */
 
  class DBAdapter
  {
    public:
			/// Constructor
      DBAdapter(DBConnection& db_con);

      /// Destructor
      ~DBAdapter();

			/// Stores a MSExperiment
			template <class ExperimentType>
			void storeExperiment(ExperimentType& exp);

			/// Reads a MSExperiment
			template <class ExperimentType>
			void loadExperiment(UID id, ExperimentType& exp);

			/// Reads a MSSpectrum
			template <class SpectrumType>
			void loadSpectrum(UID id, SpectrumType& spec);

    private:
    	/// Reference to the DB connection handed over in the constructor
    	DBConnection& db_con_;
			
			/// Not implemented, thus private
			DBAdapter();
			
			/**
				@brief Stores, updates or deletes MetaInfo data
				
				@return the id of the new MetaInfo table raw
			*/
			UID storeMetaInfo_(const String& parent_table, UID parent_id, const MetaInfoInterface& info);
			
			/**
				@brief Loads MetaInfo data
				
			*/
			void loadMetaInfo_(UID id, MetaInfoInterface& info);
  };


//------------------------------------------- IMPLEMENTATION OF TEMPLATE METHODS ----------------------------------

	template <class ExperimentType>
	void DBAdapter::storeExperiment(ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String end;              // end of the query that is added afer all fields
		String tmp;              // temporary data
		bool new_entry;          // stores if the current object is already in the DB
		QSqlQuery result;        // place to store the query results in
		int parent_id;           // stores parent_id of meta information
		
		//----------------------------------------------------------------------------------------
		//------------------------------- EXPERIMENTAL SETTINGS ---------------------------------- 
		//----------------------------------------------------------------------------------------
		
		new_entry = (exp.getPersistenceId()==0);
		if (new_entry)
		{
			query << "INSERT INTO META_MSExperiment SET ";
			end = "";
		}
		else
		{
			query << "UPDATE META_MSExperiment SET ";
			end = " WHERE id='" + String(exp.getPersistenceId()) + "'";
		}
		//type
		query << "Type=" << (1u+exp.getType());
		//date
		exp.getDate().get(tmp);
		query << ",Date='" << tmp << "'";
		//description
		query << ",Description='" << exp.getComment() << "'";
		
		query << end;
		db_con_.executeQuery(query.str(),result);
		if (new_entry)
		{
			exp.setPersistenceId(db_con_.getAutoId());
		}
		
		storeMetaInfo_("META_MSExperiment",exp.getPersistenceId(), exp);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- SPECTRUM ---------------------------------------- 
		//----------------------------------------------------------------------------------------
		for (typename ExperimentType::Iterator exp_it = exp.begin(); exp_it != exp.end(); ++exp_it)
		{
			query.str("");
			new_entry = (exp_it->getPersistenceId()==0);
			if (new_entry)
			{
				query << "INSERT INTO DATA_Spectrum SET ";
				end = "";
			}
			else
			{
				query << "UPDATE DATA_Spectrum SET ";
				end = " WHERE id='" + String (exp_it->getPersistenceId()) + "'";
			}
			//FC (MSExperiment)
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "'";
			//type
			query << ",Type=" << (1u+exp_it->getType());
			//RT
			query << ",RetentionTime='" << exp_it->getRetentionTime() << "'";
			//MS-Level
			query << ",MSLevel='" << exp_it->getMSLevel() << "'";
			//Description
			query << ",Description='" << exp_it->getComment() << "'";
			
			//TODO: MassType (average/monoisotopic)
			//TODO: TIC
			
			query << end;
			db_con_.executeQuery(query.str(),result);
			if (new_entry)
			{
				exp_it->setPersistenceId(db_con_.getAutoId());
			}
			storeMetaInfo_("DATA_Spectrum",exp_it->getPersistenceId(), *exp_it);
			
			//----------------------------------------------------------------------------------------
			//-------------------------------------- PRECURSOR --------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			if (exp_it->getMSLevel()>1)
			{
				query.str("");
				if (new_entry)
				{
					query << "INSERT INTO DATA_Precursor SET fid_Spectrum='" << exp_it->getPersistenceId() << "',";
					end = "";
				}
				else
				{
					// first determine parent_id for meta information.
					// TODO: determination of parent_id can be optimized as soon as precursor becomes a PersistentObject
					query << "SELECT id FROM DATA_Precursor WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
					db_con_.executeQuery(query.str(),result);
					result.first();
					parent_id = result.value(0).asInt();
					
					query.str("");
					query << "UPDATE DATA_Precursor SET ";
					end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
				}
				//Intensity
				query << "Intensity='" << exp_it->getPrecursorPeak().getIntensity() << "'";
				//mz
				query << ",mz='" << exp_it->getPrecursorPeak().getPosition()[0] << "'";
				//charge
				query << ",Charge='" << exp_it->getPrecursorPeak().getCharge() << "'";
				//activation method
				query << ",ActivationMethod=" << (1u+exp_it->getPrecursor().getActivationMethod());		
				//activation energy unit
				query << ",ActivationEnergyUnit=" << (1u+exp_it->getPrecursor().getActivationEnergyUnit());		
				//activation energy
				query << ",ActivationEnergy='" << exp_it->getPrecursor().getActivationEnergy() << "'";
				
				query << end;
				db_con_.executeQuery(query.str(),result);
				if (new_entry) parent_id = db_con_.getAutoId();
				storeMetaInfo_("DATA_Precursor",parent_id, exp_it->getPrecursor());
				//TODO store persistence ID => Precusor class a persistent object
			}
			
			//----------------------------------------------------------------------------------------
			//---------------------------------------- PEAKS ----------------------------------------- 
			//----------------------------------------------------------------------------------------
			query.str("");
			tmp = String(exp_it->getPersistenceId());
			db_con_.executeQuery("DELETE FROM DATA_Peak WHERE fid_Spectrum='" + tmp + "'",result);
			query << "INSERT INTO DATA_Peak (fid_Spectrum,Intensity,mz) VALUES ";
			tmp = String("('")+tmp+"','";
			for (typename ExperimentType::SpectrumType::Iterator spec_it = exp_it->begin(); spec_it != exp_it->end(); ++spec_it)
			{
				//FC (Spectrum)
				query << tmp;
				//Intensity
				query << spec_it->getIntensity() << "','";
				//mz
				query << spec_it->getPosition()[0] << "'),";
			}
			db_con_.executeQuery(String(query.str()).substr(0,-1),result);
			// TODO: call storeMetaInfo_() and loadMetaInfo_() for each Peak
			// storeMetaInfo_("DATA_Peak", parent_id, spec_it);
			
			//----------------------------------------------------------------------------------------
			//---------------------------------- INSTRUMENT SETTINGS --------------------------------- 
			//----------------------------------------------------------------------------------------
			
			const InstrumentSettings & settings = exp_it->getInstrumentSettings();
			
			query.str("");
			
			if (new_entry)
			{
				query << "INSERT INTO META_InstrumentSettings SET fid_Spectrum=" << exp_it->getPersistenceId() << ",";
				end = "";
			}
			else
			{
				query << "SELECT id FROM META_InstrumentSettings WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				db_con_.executeQuery(query.str(),result);
				result.first();
				parent_id = result.value(0).asInt();
					
				query.str("");
				query << "UPDATE META_InstrumentSettings SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}

			query << "MZRangeBegin=" << settings.getMzRangeStart() << ",";
			query << "MZRangeEnd=" << settings.getMzRangeStop() << ",";
			query << "Polarity=" << (1u+settings.getPolarity()) << ",";
			query << "ScanMode=" << (1u+settings.getScanMode());
			query << end;
			
			db_con_.executeQuery(query.str(),result);
		
			if (new_entry) parent_id = db_con_.getAutoId();
			storeMetaInfo_("META_InstrumentSettings", parent_id, exp_it->getInstrumentSettings());
		}
	}

	template <class ExperimentType>
	void DBAdapter::loadExperiment(UID id, ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String tmp;              // temporary data
		QSqlQuery result;        // place to store the query results in
		
		query << "SELECT Type-1,Date,fid_MetaInfo,Description FROM META_MSExperiment WHERE id='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		//Experiment meta info
		exp.setType((ExperimentalSettings::ExperimentType)(result.value(0).asInt()));
		if (result.value(1).asDate().isValid())
		{
			Date d;
			d.set(result.value(1).asDate().toString(Qt::ISODate).ascii());
			exp.setDate(d);
		}
		exp.setComment(result.value(3).asString().ascii());
		loadMetaInfo_(result.value(2).asInt(),exp);
		
		//spectra
		query.str("");
		query << "SELECT id FROM DATA_Spectrum WHERE fid_MSExperiment='" << id << "' ORDER BY id ASC";
		
		db_con_.executeQuery(query.str(),result);
		exp.resize(result.size());
		UnsignedInt i = 0;
		result.first();
		while (result.isValid())
		{
			loadSpectrum(result.value(0).asInt(),exp[i]);
			++i;
			result.next();
		}
		
		//id
		exp.setPersistenceId(id);
	}

	template <class SpectrumType>
	void DBAdapter::loadSpectrum(UID id, SpectrumType& spec)
	{
		spec = SpectrumType();
		
		std::stringstream query; // query to build
		QSqlQuery result;        // place to store the query results in
		InstrumentSettings settings; // stores settings that are read from DB
		
		query << "SELECT Type-1,RetentionTime,MSLevel,Description,fid_MetaInfo FROM DATA_Spectrum WHERE id='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		//Spectrum meta info
		spec.setType((SpectrumSettings::SpectrumType)(result.value(0).asInt()));
		spec.setRetentionTime(result.value(1).toDouble());
		spec.setMSLevel(result.value(2).asInt());
		spec.setComment(result.value(3).asString().ascii());
		loadMetaInfo_(result.value(4).asInt(),spec);
		
		// Instrument settings
		query.str("");
		query << "SELECT MZRangeBegin, MZRangeEnd, Polarity-1, ScanMode-1, fid_MetaInfo FROM META_InstrumentSettings WHERE id='" << id << "'";
		db_con_.executeQuery(query.str(),result);
		result.first();
		
		settings.setMzRangeStart(result.value(0).toDouble());
		settings.setMzRangeStop(result.value(1).toDouble());
		settings.setPolarity((IonSource::Polarity) (result.value(2).asInt()));
		settings.setScanMode((InstrumentSettings::ScanMode) (result.value(3).asInt()));
		spec.setInstrumentSettings(settings);
		loadMetaInfo_(result.value(4).asInt(),spec.getInstrumentSettings());

		//precursor
		
		if(spec.getMSLevel()>1)
		{
			query.str("");
			query << "SELECT mz,Intensity,Charge,ActivationMethod-1,ActivationEnergyUnit-1,ActivationEnergy,fid_MetaInfo FROM DATA_Precursor WHERE fid_Spectrum='" << id << "'";
			db_con_.executeQuery(query.str(),result);
			result.first();
			spec.getPrecursorPeak().getPosition()[0] = (result.value(0).toDouble());
			spec.getPrecursorPeak().setIntensity(result.value(1).toDouble());
			spec.getPrecursorPeak().setCharge(result.value(2).asInt());
			spec.getPrecursor().setActivationMethod((Precursor::ActivationMethod)(result.value(3).asInt()));
			spec.getPrecursor().setActivationEnergyUnit((Precursor::EnergyUnits)(result.value(4).asInt()));
			spec.getPrecursor().setActivationEnergy(result.value(5).toDouble());
			loadMetaInfo_(result.value(6).asInt(),spec.getPrecursor());
		}
		
		//Peaks
		
		query.str("");
		query << "SELECT mz,Intensity FROM DATA_Peak WHERE fid_Spectrum='" << id << "' ORDER BY mz ASC";
		db_con_.executeQuery(query.str(),result);

		typename SpectrumType::PeakType p;
		result.first();
		while(result.isValid())
		{
			p.getPosition()[0] = result.value(0).toDouble();
			p.setIntensity(result.value(1).toDouble());
			spec.push_back(p);
			result.next();
		}
		
		//id
		spec.setPersistenceId(id);
	}
}
#endif
