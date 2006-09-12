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
    @todo add missing members that are only in the DB schema (Marc)
    @todo intensive testing: empty map, all members set, update only, ... (Marc)
    @todo add support for missing classes and MetaInfo(Marc)
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
    	///reference to the DB connection handed over in the constructor
    	DBConnection& db_con_;
			
			/// Not implemented, thus private
			DBAdapter();
			
			UID storeMetaInfo(const String& parent_table, UID parent_id, const MetaInfoInterface& info);
  };

	template <class ExperimentType>
	void DBAdapter::storeExperiment(ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String end;              // end of the query that is added afer all fields
		String tmp;              // temporary data
		bool new_entry;          // stores if the current object is already in the DB
		
		//----------------------------------------------------------------------------------------
		//------------------------------- EXPERIMENTAL SETTINGS ---------------------------------- 
		//----------------------------------------------------------------------------------------
		
		new_entry = (exp.getPersistenceId()==0);
		if (new_entry)
		{
			query << "INSERT INTO META_MSExperiment SET ";
			end = "'";
		}
		else
		{
			query << "UPDATE META_MSExperiment SET ";
			end = "' WHERE id='" + String(exp.getPersistenceId()) + "'";
		}
		//type
		query << "Type=" << (1u+exp.getType());
		//date
		exp.getDate().get(tmp);
		query << ",Date='" << tmp;
		
		//TODO: Description
		
		query << end;
		db_con_.executeQuery(query.str());
		if (new_entry)
		{
			exp.setPersistenceId(db_con_.getAutoId());
		}
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
				end = "'";
			}
			else
			{
				query << "UPDATE DATA_Spectrum SET ";
				end = "' WHERE id='" + String(exp.getPersistenceId()) + "'";
			}
			//FC (MSExperiment)
			query << "fid_MSExperiment='" << exp.getPersistenceId();
			//type
			query << "',Type=" << (1u+exp_it->getType());
			//RT
			query << ",RetentionTime='" << exp_it->getRetentionTime();
			//MS-Level
			query << "',MSLevel='" << exp_it->getMSLevel();
			
			//TODO: Description
			//TODO: MassType (average/monoisotopic)
			//TODO: TIC
			
			query << end;
			db_con_.executeQuery(query.str());
			if (new_entry)
			{
				exp_it->setPersistenceId(db_con_.getAutoId());
			}
			
			//----------------------------------------------------------------------------------------
			//-------------------------------------- PRECURSOR --------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			if (exp_it->getMSLevel()>1)
			{
				query.str("");
				if (new_entry)
				{
					query << "INSERT INTO DATA_Precursor SET fid_Spectrum='" << exp_it->getPersistenceId();
					end = "'";
				}
				else
				{
					query << "UPDATE DATA_Precursor SET ";
					end = "' WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
				}
				//Intensity
				query << "',Intensity='" << exp_it->getPrecursorPeak().getIntensity();
				//mz
				query << "',mz='" << exp_it->getPrecursorPeak().getPosition()[0];
				//charge
				query << "',Charge='" << exp_it->getPrecursorPeak().getCharge();				
				//activation method
				query << "',ActivationMethod=" << (1u+exp_it->getPrecursor().getActivationMethod());		
				//activation energy unit
				query << ",ActivationEnergyUnit=" << (1u+exp_it->getPrecursor().getActivationEnergyUnit());		
				//activation energy
				query << ",ActivationEnergy='" << exp_it->getPrecursor().getActivationEnergy();		
				
				query << end;
				db_con_.executeQuery(query.str());
			}
			//TODO MSLevel change from 2 to 1 => delete entry
			
			//----------------------------------------------------------------------------------------
			//---------------------------------------- PEAKS ----------------------------------------- 
			//----------------------------------------------------------------------------------------
			query.str("");
			tmp = String(exp_it->getPersistenceId());
			db_con_.executeQuery("DELETE FROM DATA_Peak WHERE fid_Spectrum='" + tmp + "'");
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
			db_con_.executeQuery(String(query.str()).substr(0,-1));					
		}
	}

	template <class ExperimentType>
	void DBAdapter::loadExperiment(UID id, ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String tmp;              // temporary data
		
		query << "SELECT Type-1,Date FROM META_MSExperiment WHERE id='" << id << "'";
		db_con_.executeQuery(query.str());
		//scope for res...
		{ 		
			QSqlQuery& res(db_con_.lastResult());
			res.first();
			
			//Experiment meta info
			exp.setType((ExperimentalSettings::ExperimentType)(res.value(0).toInt()));
			if (res.value(1).asDate().isValid())
			{
				Date d;
				d.set(res.value(1).asDate().toString(Qt::ISODate).ascii());
				exp.setDate(d);
			}
		}
		
		//spectra
		query.str("");
		query << "SELECT id FROM DATA_Spectrum WHERE fid_MSExperiment='" << id << "' ORDER BY id ASC";
		db_con_.executeQuery(query.str());
		//scope for res...
		{ 
			QSqlQuery res(db_con_.lastResult()); //copy as it is overwritten otherwise
			exp.resize(res.size());
			UnsignedInt i = 0;
			res.first();
			while (res.isValid())
			{
				loadSpectrum(res.value(0).toInt(),exp[i]);
				++i;
				res.next();
			}
		}
		
		//id
		exp.setPersistenceId(id);
	}

	template <class SpectrumType>
	void DBAdapter::loadSpectrum(UID id, SpectrumType& spec)
	{
		std::stringstream query; // query to build
		
		query << "SELECT Type-1,RetentionTime,MSLevel FROM DATA_Spectrum WHERE id='" << id << "'";
		db_con_.executeQuery(query.str());
		//scope for res...
		{
			QSqlQuery& res(db_con_.lastResult());
			res.first();
			
			//Spectrum meta info
			spec.setType((SpectrumSettings::SpectrumType)(res.value(0).toInt()));
			spec.setRetentionTime(res.value(1).toDouble());
			spec.setMSLevel(res.value(2).toInt());
		}
		
		//precursor
		
		if(spec.getMSLevel()>1)
		{
			query.str("");
			query << "SELECT mz,Intensity,Charge,ActivationMethod-1,ActivationEnergyUnit-1,ActivationEnergy FROM DATA_Precursor WHERE fid_Spectrum='" << id << "'";		
			db_con_.executeQuery(query.str());
			QSqlQuery& res(db_con_.lastResult());
			res.first();
			spec.getPrecursorPeak().getPosition()[0] = (res.value(0).toDouble());
			spec.getPrecursorPeak().setIntensity(res.value(1).toDouble());
			spec.getPrecursorPeak().setCharge(res.value(2).toInt());
			spec.getPrecursor().setActivationMethod((Precursor::ActivationMethod)(res.value(3).toInt()));
			spec.getPrecursor().setActivationEnergyUnit((Precursor::EnergyUnits)(res.value(4).toInt()));
			spec.getPrecursor().setActivationEnergy(res.value(5).toDouble());
		}
		
		//Peaks
		
		query.str("");
		query << "SELECT mz,Intensity FROM DATA_Peak WHERE fid_Spectrum='" << id << "' ORDER BY mz ASC";
		db_con_.executeQuery(query.str());
		//scope for res...
		{ 
			QSqlQuery& res(db_con_.lastResult());
			typename SpectrumType::PeakType p;
			res.first();
			while(res.isValid())
			{
				p.getPosition()[0] = res.value(0).toDouble();
				p.setIntensity(res.value(1).toDouble());
				spec.push_back(p);
				res.next();
			}
		}
		
		//id
		spec.setPersistenceId(id);
	}
}
#endif
