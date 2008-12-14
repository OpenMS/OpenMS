// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_DB_DBADAPTER_H
#define OPENMS_FORMAT_DB_DBADAPTER_H

//OpenMS includes
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>

//QT includes
#include <QtSql/QSqlQuery>
#include <QtCore/QVariant>
#include <QtCore/QDate>

//std and STL includes
#include <string>
#include <map>

namespace OpenMS
{	
  /** 
  	@brief A class for accessing and storing data in a SQL database
  	
    It can be used to create objects from the DB or store them in the DB.
		
		@todo Check if METADATA implementation is complete - all members and MetaInfoInterface (HiWi)
		
    @ingroup DatabaseIO
  */
 
  class OPENMS_DLLAPI DBAdapter
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

			template <class SpectrumType>
			/// Reads a MSSpectrum
			void loadSpectrum(UID id, SpectrumType& spec);

      /// Mutable access to the options for loading/storing 
      PeakFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing 
      const PeakFileOptions& getOptions() const;

			/**
				@brief Returns true if the DB is up-to-date (Checks the version in ADMIN_Version table).
				
				@param warning if this is set, a warning is issued to stderr if the db is not up-to-date. 
			*/
			bool checkDBVersion(bool warning);
			
			/// Deletes all tables in the database and creates a new OpenMS database
			void createDB();
			
    private:
    	/// Reference to the DB connection handed over in the constructor
    	DBConnection& db_con_;
			
			/// Not implemented
			DBAdapter();
			
			/**
				@brief Stores, updates or deletes MetaInfo data
				
				@return the id of the new META_MetaInfo table row
			*/
			UID storeMetaInfo_(const String& parent_table, UID parent_id, const MetaInfoInterface& info);	
			/// Overload for Peak1D that does nothing
			UID storeMetaInfo_(const String& parent_table, UID parent_id, const Peak1D& peak);	
			///Overloaded method for RichPeak1D, which is both a MetaInfoInterface and a Peak1D
			UID storeMetaInfo_(const String& parent_table, UID parent_id, const RichPeak1D& peak);	
			
			///Loads MetaInfo data from database
			void loadMetaInfo_(UID id, MetaInfoInterface& info);
			///Overloaded method for Peak1D, which does nothing
			void loadMetaInfo_(UID id, Peak1D& peak);
			///Overloaded method for RichPeak1D, which is both a MetaInfoInterface and a Peak1D
			void loadMetaInfo_(UID id, RichPeak1D& peak);
			
			/**
				@brief Conditionally deletes MetaInfo data from database
				
			*/
			void deleteMetaInfo_(const String& parent_table, const String& condition);
			
			/**
				@brief Stores, updates or deletes file information
				
				@return the id of the new META_File table row
			*/
			UID storeFile_(const String& parent_table, UID parent_id, const SourceFile& file);
			
			/**
				@brief Loads file information
				
			*/
			void loadFile_(UID id, SourceFile& file);
  
  			/**
				@brief Stores, updates or deletes sample information
				
				@return the id of the new META_Sample table row
			*/
			UID storeSample_(const Sample& sample, UID exp_id, UID parent_id);
			
			/**
				@brief Loads sample information
				
			*/
			void loadSample_(UID id, Sample& sample);

			PeakFileOptions options_;
};


//------------------------------------------- IMPLEMENTATION OF TEMPLATE METHODS ----------------------------------

	template <class ExperimentType>
	void DBAdapter::storeExperiment(ExperimentType& exp)
	{
		std::stringstream query; // query to build
		String end;              // end of the query that is added afer all fields
		String tmp;              // temporary data
		bool new_entry(false);          // stores if the current object is already in the DB
		QSqlQuery result;        // place to store the query results in
		int parent_id(-1);           // stores parent_id of meta information
		UID acquisition_info_id(0); // stores id of acquisition_info
		UID meta_id(0);             // stores MetaInfo id of meta information that was just stored
		UID meta_parent_id(0);      // stores parent ID of MetaInfo that will be stored

		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		
		
		//----------------------------------------------------------------------------------------
		//------------------------------- EXPERIMENT --------------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		query.str("");
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
		//date
		//TODO make sure date and time are stored/loaded correctly
		query << "Date='" << exp.getDateTime().get() << "'";
		//description
		query << ",Description='" << exp.getComment() << "'";
		
		query << end;
		result = db_con_.executeQuery(query.str());
		if (new_entry)
		{
			exp.setPersistenceId(db_con_.getAutoId());
		}
		
		storeMetaInfo_("META_MSExperiment", exp.getPersistenceId(), exp);
		
		//----------------------------------------------------------------------------------------
		//------------------------------- PROTEIN IDENTIFICATIONS / HITS-------------------------- 
		//----------------------------------------------------------------------------------------

		std::vector<ProteinIdentification>& pi = exp.getProteinIdentifications();
		
		// first delete all old values, no matter whether we're updating or not
		// do we have to delete the MetaInfo as well?
		query.str("");
		query << "DELETE FROM ID_ProteinIdentification WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		result = db_con_.executeQuery(query.str());

		for (std::vector<ProteinIdentification>::const_iterator pi_it = pi.begin(); pi_it != pi.end(); pi_it++)
		{
			query.str("");
			query << "INSERT INTO ID_ProteinIdentification SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "'";
			query << ",SearchEngine='" << pi_it->getSearchEngine() << "'";
			query << ",SearchEngineVersion='" << pi_it->getSearchEngineVersion() << "'";
			query << ",Date='" << pi_it->getDateTime().get() << "'";
			query << ",ScoreType='" << pi_it->getScoreType() << "'";
			query << ",HigherScoreBetter='" << pi_it->isHigherScoreBetter() << "'";
			query << ",SignificanceThreshold='" << pi_it->getSignificanceThreshold() << "'";
			
			result = db_con_.executeQuery(query.str());
			parent_id = db_con_.getAutoId();
			
			storeMetaInfo_("ID_ProteinIdentification", parent_id, *pi_it);
			//needs getter and setter methods first
			//storeFile_("ID_ProteinIdentification", parent_id, pi_it->getSourceFile());

			for (std::vector<ProteinHit>::const_iterator ph_it = pi_it->getHits().begin(); ph_it != pi_it->getHits().end(); ph_it++)
			{
				query.str("");
				query << "INSERT INTO ID_ProteinHit SET ";
				query << "fid_ProteinIdentification='" << parent_id << "'";
				query << ",Score='" << ph_it->getScore() << "'";
				query << ",Accession='" << ph_it->getAccession() << "'";
				query << ",Sequence='" << ph_it->getSequence() << "'";
				
				result = db_con_.executeQuery(query.str());
				meta_id = db_con_.getAutoId();
				
				storeMetaInfo_("ID_ProteinHit", meta_id, *ph_it);
			}
		}
		
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- SAMPLE ------------------------------------------ 
		//----------------------------------------------------------------------------------------
		
		query.str("");
		deleteMetaInfo_("META_Sample", "fid_MSExperiment=" + String(exp.getPersistenceId()));
		// this also deletes all references in META_SampleTreatment, META_Digestion and META_Modification by constraint
		query << "DELETE FROM META_Sample WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		storeSample_(exp.getSample(), exp.getPersistenceId(), 0);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- CONTACTPERSON ----------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const std::vector<ContactPerson>& contacts = exp.getContacts();
		
		query.str("");
		deleteMetaInfo_("META_ContactPerson", "fid_MSExperiment=" + String(exp.getPersistenceId()));
		query << "DELETE FROM META_ContactPerson WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		result = db_con_.executeQuery(query.str());
			
		for (std::vector<ContactPerson>::const_iterator contact_it = contacts.begin(); contact_it != contacts.end(); contact_it++)
		{
			query.str("");
			query << "INSERT INTO META_ContactPerson SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "'";
			query << ",PreName='" << contact_it->getFirstName() << "'";
			query << ",LastName='" << contact_it->getLastName() << "'";
			query << ",Affiliation='" << contact_it->getInstitution() << "'";
			query << ",Email='" << contact_it->getEmail() << "'";
			query << ",Comment='" << contact_it->getContactInfo() << "'";
			
			result = db_con_.executeQuery(query.str());
			parent_id = db_con_.getAutoId();
			
			storeMetaInfo_("META_ContactPerson", parent_id, *contact_it);
		}
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- HPLC -------------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const HPLC& hplc = exp.getHPLC();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_HPLC SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_HPLC WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());
			result.first();
			parent_id = result.value(0).toInt();
			
			query.str("");
			query << "UPDATE META_HPLC SET ";
			end = " WHERE fid_MSExperiment='" + String (exp.getPersistenceId()) + "'";
		}
		
		query << "InstrumentName='" << hplc.getInstrument() << "'";
		query << ",ColumnName='" << hplc.getColumn() << "'";
		query << ",Description='" << hplc.getComment() << "'";
		query << ",Flux=" << hplc.getFlux();
		query << ",Pressure=" << hplc.getPressure();
		query << ",Temperature=" << hplc.getTemperature();
		
		query << end;
		result = db_con_.executeQuery(query.str());
		
		if (new_entry)
		{
			parent_id = db_con_.getAutoId();
		}
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- GRADIENT* --------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const Gradient& gradient = exp.getHPLC().getGradient();
		const std::vector <String>& eluents = gradient.getEluents();
		const std::vector <Int>& time = gradient.getTimepoints();
		const std::vector< std::vector< UInt > >& percentages = gradient.getPercentages();
		std::stringstream query_eluents, query_time, query_percentages;
		UID eluents_id(0), time_id(0);
		
		// this also deletes all references in META_GradientPercentage by constraint
		query.str("");
		query << "DELETE FROM META_GradientEluent WHERE fid_HPLC=" << parent_id;
		result = db_con_.executeQuery(query.str());
		query.str("");
		query << "DELETE FROM META_GradientTime WHERE fid_HPLC=" << parent_id;
		result = db_con_.executeQuery(query.str());
		
		if (! eluents.empty())
		{
			query_eluents.str("");
			query_eluents << "INSERT INTO META_GradientEluent (fid_HPLC, Name) VALUES ";
			for (std::vector<String>::const_iterator eluents_it = eluents.begin(); eluents_it != eluents.end(); eluents_it++)
			{
				query_eluents << "(";
				query_eluents << parent_id;
				query_eluents << ",'" << *eluents_it << "'";
				query_eluents << "),";
			}
			result = db_con_.executeQuery(String(query_eluents.str()).substr(0,-1));
			eluents_id = db_con_.getAutoId();
		}
		
		if (! time.empty())
		{
			query_time.str("");
			query_time << "INSERT INTO META_GradientTime (fid_HPLC, Time) VALUES ";
			for (std::vector<Int>::const_iterator time_it = time.begin(); time_it != time.end(); time_it++)
			{
				query_time << "(";
				query_time << parent_id;
				query_time << "," << *time_it;
				query_time << "),";
			}
			result = db_con_.executeQuery(String(query_time.str()).substr(0,-1));
			time_id = db_con_.getAutoId();
		}
		
		if (! percentages.empty() && ! eluents.empty() && ! time.empty())
		{
			query_percentages.str("");
			query_percentages << "INSERT INTO META_GradientPercentage (fid_GradientEluent, fid_GradientTime, Percentage) VALUES ";
			int i = 0;
			// iterate over eluents
			for (std::vector< std::vector< UInt> >::const_iterator percent_outer_it = percentages.begin(); percent_outer_it != percentages.end(); percent_outer_it++)
			{
				int j = 0;
				// iterate over timepoints
				for (std::vector< UInt>::const_iterator percent_inner_it = (*percent_outer_it).begin(); percent_inner_it != (*percent_outer_it).end(); percent_inner_it++)
				{
					query_percentages << "(";
					query_percentages << eluents_id + i;
					query_percentages << "," << time_id + j;
					query_percentages << "," << *percent_inner_it;
					query_percentages << "),";
					j++;
				}
				i++;
			}
			result = db_con_.executeQuery(String(query_percentages.str()).substr(0,-1));
		}
			
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- INSTRUMENT -------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		const Instrument& instrument = exp.getInstrument();
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_MSInstrument SET ";
			query << "fid_MSExperiment='" << exp.getPersistenceId() << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_MSInstrument WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());
			result.first();
			parent_id = result.value(0).toInt();
			
			query.str("");
			query << "UPDATE META_MSInstrument SET ";
			end = " WHERE fid_MSExperiment='" + String (exp.getPersistenceId()) + "'";
		}
		
		query << "Model='" << instrument.getModel() << "'";
		query << ",Vendor='" << instrument.getVendor() << "'";
		query << ",Description='" << instrument.getCustomizations() << "'";
		
		query << end;
		result = db_con_.executeQuery(query.str());
		
		if (new_entry)
		{
			parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_MSInstrument", parent_id, instrument);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- IONDETECTOR ------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		//TODO Hack: Multiple ion detectors should be supported
		exp.getInstrument().getIonDetectors().resize(1);
		const IonDetector& detector = exp.getInstrument().getIonDetectors()[0];
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_IonDetector SET ";
			query << "fid_MSInstrument='" << parent_id << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
			result = db_con_.executeQuery(query.str());
			result.first();
			meta_parent_id = result.value(0).toInt();
			
			query.str("");
			query << "UPDATE META_IonDetector SET ";
			end = " WHERE fid_MSInstrument='" + String (parent_id) + "'";
		}
		
		query << "AcquisitionMode=" << (1u+detector.getAcquisitionMode());
		query << ",Type=" << (1u+detector.getType());
		query << ",Resolution=" << detector.getResolution();
		query << ",ADCSamplingFrequency=" << detector.getADCSamplingFrequency();
		
		query << end;
		result = db_con_.executeQuery(query.str());
		
		if (new_entry)
		{
			meta_parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_IonDetector", meta_parent_id, detector);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- IONSOURCE --------------------------------------- 
		//----------------------------------------------------------------------------------------
		
		//TODO Hack: Multiple ion sources should be supported
		exp.getInstrument().getIonSources().resize(1);
		const IonSource& source = exp.getInstrument().getIonSources()[0];
		query.str("");
		
		if (new_entry)
		{
			query << "INSERT INTO META_IonSource SET ";
			query << "fid_MSInstrument='" << parent_id << "',";
			end = "";
		}
		else
		{
			query << "SELECT id FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
			result = db_con_.executeQuery(query.str());
			
			query.str("");
			query << "UPDATE META_IonSource SET ";
			end = " WHERE fid_MSInstrument='" + String (parent_id) + "'";
			
			result.first();
			meta_parent_id = result.value(0).toInt();
		}
		
		query << "InletType=" << (1u+source.getInletType());
		query << ",IonizationMethod=" << (1u+source.getIonizationMethod());
		query << ",IonizationMode=" << (1u+source.getPolarity());
		
		query << end;
		result = db_con_.executeQuery(query.str());
		
		if (new_entry)
		{
			meta_parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_IonSource", meta_parent_id, source);
		
		//----------------------------------------------------------------------------------------
		//-------------------------------------- MASSANALYZER ------------------------------------ 
		//----------------------------------------------------------------------------------------
		
		const std::vector<MassAnalyzer>& analyzers = exp.getInstrument().getMassAnalyzers();
		query.str("");
		
		deleteMetaInfo_("META_MassAnalyzer", "fid_MSInstrument=" + String(parent_id));
		query << "DELETE FROM META_MassAnalyzer WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
			
		for (std::vector<MassAnalyzer>::const_iterator analyzer_it = analyzers.begin(); analyzer_it != analyzers.end(); analyzer_it++)
		{
			query.str("");
			query << "INSERT INTO META_MassAnalyzer SET ";
			query << "fid_MSInstrument='" << parent_id << "'";
			query << ",Accuracy=" << analyzer_it->getAccuracy();
			query << ",FinalMSExponent=" << analyzer_it->getFinalMSExponent();
			query << ",IsolationWidth=" << analyzer_it->getIsolationWidth();
			query << ",MagneticFieldStrength=" << analyzer_it->getMagneticFieldStrength();
			query << ",ReflectronState=" << (1u+analyzer_it->getReflectronState());
			query << ",Resolution=" << analyzer_it->getResolution();
			query << ",ResolutionMethod=" << (1u+analyzer_it->getResolutionMethod());
			query << ",ResolutionType=" << (1u+analyzer_it->getResolutionType());
			query << ",ScanDirection=" << (1u+analyzer_it->getScanDirection());
			query << ",ScanLaw=" << (1u+analyzer_it->getScanLaw());
			query << ",ScanRate=" << analyzer_it->getScanRate();
			query << ",ScanTime=" << analyzer_it->getScanTime();
			query << ",TOFPathLength=" << analyzer_it->getTOFTotalPathLength();
			query << ",Type=" << (1u+analyzer_it->getType());
			
			result = db_con_.executeQuery(query.str());
			storeMetaInfo_("META_MassAnalyzer", db_con_.getAutoId(), *analyzer_it);
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
			query << ",RetentionTime='" << exp_it->getRT() << "'";
			//MS-Level
			query << ",MSLevel='" << exp_it->getMSLevel() << "'";
			//Description
			query << ",Description='" << exp_it->getComment() << "'";
			
			//TODO: MassType (average/monoisotopic)
			//TODO: TIC
			
			
			query << end;
			result = db_con_.executeQuery(query.str());
			if (new_entry)
			{
				exp_it->setPersistenceId(db_con_.getAutoId());
			}
			storeFile_("DATA_Spectrum", exp_it->getPersistenceId(), exp_it->getSourceFile());
			meta_id = storeMetaInfo_("DATA_Spectrum", exp_it->getPersistenceId(), *exp_it);
			

			//----------------------------------------------------------------------------------------
			//------------------------------- PEPTIDE IDENTIFICATIONS / HITS-------------------------- 
			//----------------------------------------------------------------------------------------
	
			std::vector<PeptideIdentification>& pei = exp_it->getPeptideIdentifications();

			// first delete all old values, no matter whether we're updating or not
			// do we have to delete the MetaInfo as well?
			query.str("");
			query << "DELETE FROM ID_PeptideIdentification WHERE fid_Spectrum='";
			query << exp_it->getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());
			
			for (std::vector<PeptideIdentification>::const_iterator pei_it = pei.begin(); pei_it != pei.end(); pei_it++)
			{
				query.str("");
				query << "INSERT INTO ID_PeptideIdentification SET ";
				query << "fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				query << ",SignificanceThreshold='" << pei_it->getSignificanceThreshold() << "'";
				query << ",ScoreType='" << pei_it->getScoreType() << "'";
				query << ",HigherScoreBetter='" << pei_it->isHigherScoreBetter() << "'";
				
				result = db_con_.executeQuery(query.str());
				parent_id = db_con_.getAutoId();
				
				storeMetaInfo_("ID_PeptideIdentification", parent_id, *pei_it);
				//needs getter and setter methods first
				//storeFile_("ID_PeptideIdentification", parent_id, pei_it->getSourceFile());
	
				for (std::vector<PeptideHit>::const_iterator peh_it = pei_it->getHits().begin(); peh_it != pei_it->getHits().end(); peh_it++)
				{
					query.str("");
					query << "INSERT INTO ID_PeptideHit SET ";
					query << "fid_Identification='" << parent_id << "'";
					query << ",Score='" << peh_it->getScore() << "'";
					query << ",charge='" << peh_it->getCharge() << "'";
					query << ",Sequence='" << peh_it->getSequence() << "'";
					query << ",AABefore='" << peh_it->getAABefore() << "'";
					query << ",AAAfter='" << peh_it->getAAAfter() << "'";
					
					result = db_con_.executeQuery(query.str());
					meta_id = db_con_.getAutoId();
					
					storeMetaInfo_("ID_PeptideHit", meta_id, *peh_it);
				}
			}
	
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
					result = db_con_.executeQuery(query.str());
					result.first();
					parent_id = result.value(0).toInt();
					
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
				//window size
				query << ",WindowSize='" << exp_it->getPrecursor().getWindowSize() << "'";
				
				query << end;
				result = db_con_.executeQuery(query.str());
				if (new_entry) parent_id = db_con_.getAutoId();
				storeMetaInfo_("DATA_Precursor",parent_id, exp_it->getPrecursor());
				//TODO store persistence ID => Precusor class a persistent object
			}
			
			//----------------------------------------------------------------------------------------
			//---------------------------------------- PEAKS ----------------------------------------- 
			//----------------------------------------------------------------------------------------
			query.str("");
			deleteMetaInfo_("DATA_Peak", "fid_Spectrum=" + String(exp_it->getPersistenceId()));
			result = db_con_.executeQuery("DELETE FROM DATA_Peak WHERE fid_Spectrum=" + String(exp_it->getPersistenceId()));
			if (exp_it->size()!=0)
			{
				query << "INSERT INTO DATA_Peak (fid_Spectrum,Intensity,mz) VALUES ";
				tmp = "(" + String(exp_it->getPersistenceId()) + ",'";
				for (typename ExperimentType::SpectrumType::Iterator spec_it = exp_it->begin(); spec_it != exp_it->end(); ++spec_it)
				{
					//Foreign Key (Spectrum)
					query << tmp;
					//Intensity
					query << spec_it->getIntensity() << "','";
					//mz
					query << spec_it->getPosition() << "'),";
				}
				result = db_con_.executeQuery(String(query.str()).substr(0,-1));
			}
			// We know that all inserted peaks have IDs beginning from last_insert_id() (= ID of first inserted entry
			// of last insert operation), so we can insert Meta Information without actually fetching the ID
			UID insert_id = db_con_.getAutoId();
			for (typename ExperimentType::SpectrumType::Iterator spec_it = exp_it->begin(); spec_it != exp_it->end(); ++spec_it)
			{
				storeMetaInfo_("DATA_Peak", insert_id, *spec_it);
				insert_id++;
			}
			
			//----------------------------------------------------------------------------------------
			//---------------------- METAINFODESCRIPTION / METADATAARRAYS  --------------------------- 
			//----------------------------------------------------------------------------------------

			const typename ExperimentType::SpectrumType::MetaDataArrays& meta_data_arrays = exp_it->getMetaDataArrays();
			
			for (typename ExperimentType::SpectrumType::MetaDataArrays::const_iterator mdarrays_it = meta_data_arrays.begin(); mdarrays_it != meta_data_arrays.end(); ++mdarrays_it)
			{
				// first check if there is already an entry in META_MetaInfoDescription for this spectrum and this name
				// We cannot simply delete all entries for the spectrum because this might leave unreferenced META_TNVs 
				// and META_MetaInfos in the database.
				query.str("");
				query << "SELECT id FROM META_MetaInfoDescription WHERE fid_Spectrum=";
				query << exp_it->getPersistenceId();
				query << " AND Name='" << mdarrays_it->getName() << "'";
				
				result = db_con_.executeQuery(query.str());
				
				query.str("");
				
				if (result.size() > 0)
				{
					parent_id = result.value(0).toInt();
					new_entry = false;
					query << "UPDATE META_MetaInfoDescription SET ";
					end  = " WHERE fid_Spectrum=" + String(exp_it->getPersistenceId());
					end += " AND Name='" + mdarrays_it->getName() + "'";
				}
				else
				{
					new_entry = true;
					query << "INSERT INTO META_MetaInfoDescription SET ";
					query << "fid_Spectrum=" << exp_it->getPersistenceId() << ", ";
					query << "Name='" << mdarrays_it->getName() << "',";
					end = "";
				}
				
				query << "Description='" << mdarrays_it->getComment() << "'";
				query << end;
				
				result = db_con_.executeQuery(query.str());
				if (new_entry)
				{
					parent_id = db_con_.getAutoId();
				}
				
				storeFile_("META_MetaInfoDescription", parent_id, mdarrays_it->getSourceFile());
				storeMetaInfo_("META_MetaInfoDescription", parent_id, *mdarrays_it);
				
				// store meta data contained in the MetaDataArrays
				query.str("");
				query << "DELETE FROM DATA_PeakMetaData WHERE fid_MetaInfoDescription=";
				query << parent_id;
				result = db_con_.executeQuery(query.str());
				
				query.str("");
				query << "SELECT id FROM DATA_Peak WHERE fid_Spectrum=" << exp_it->getPersistenceId();
				result = db_con_.executeQuery(query.str());
				result.first();
				
				query.str("");
				query << "INSERT INTO DATA_PeakMetaData (fid_Peak,fid_MetaInfoDescription,Value) VALUES ";
				for (typename ExperimentType::SpectrumType::MetaDataArray::const_iterator meta_array_it = mdarrays_it->begin(); meta_array_it != mdarrays_it->end(); meta_array_it++)
				{
					if(result.isValid())
					{
						query << "(" << result.value(0).toInt() << "," << parent_id << "," << *meta_array_it << "),";
						result.next();
					}
					else
					{
						break;
					}
				}
				result = db_con_.executeQuery(String(query.str()).substr(0,-1));
			}
			
			
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
				result = db_con_.executeQuery(query.str());
				result.first();
				parent_id = result.value(0).toInt();
					
				query.str("");
				query << "UPDATE META_InstrumentSettings SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}
			//TODO handle several scan windows
			if (settings.getScanWindows().size()>0)
			{
				query << "MZRangeBegin=" << settings.getScanWindows()[0].begin << ",";
				query << "MZRangeEnd=" << settings.getScanWindows()[0].end << ",";
			}
			else
			{
				query << "MZRangeBegin=0.0, MZRangeEnd=0.0,";
			}
			query << "Polarity=" << (1u+settings.getPolarity()) << ",";
			//TODO handle all types!!!
			query << "ScanMode=" << (1u+settings.getScanMode());
			query << end;
			
			result = db_con_.executeQuery(query.str());
		
			if (new_entry) parent_id = db_con_.getAutoId();
			storeMetaInfo_("META_InstrumentSettings", parent_id, exp_it->getInstrumentSettings());
			
			//----------------------------------------------------------------------------------------
			//--------------------------------- ACQUISITIONINFO -------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			const AcquisitionInfo & info = exp_it->getAcquisitionInfo();
			
			query.str("");
			
			if (new_entry)
			{
				query << "INSERT INTO META_AcquisitionInfo SET fid_Spectrum=" << exp_it->getPersistenceId() << ",";
				end = "";
			}
			else
			{
				query << "SELECT id FROM META_AcquisitionInfo WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				result = db_con_.executeQuery(query.str());
				result.first();
				acquisition_info_id = result.value(0).toInt();
					
				query.str("");
				query << "UPDATE META_AcquisitionInfo SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}

			query << "MethodOfCombination='" << info.getMethodOfCombination() << "'";
			query << end;
			
			result = db_con_.executeQuery(query.str());
			if (new_entry)
			{
				acquisition_info_id = db_con_.getAutoId();
			}
		
			//----------------------------------------------------------------------------------------
			//----------------------------------- ACQUISITION ---------------------------------------- 
			//----------------------------------------------------------------------------------------
			
			query.str("");
			deleteMetaInfo_("META_Acquisition", "fid_AcquisitionInfo=" + String(parent_id));
			query << "DELETE FROM META_Acquisition WHERE fid_AcquisitionInfo=" << parent_id;
			result = db_con_.executeQuery(query.str());
				
			for (std::vector<Acquisition>::const_iterator info_it = info.begin(); info_it != info.end(); info_it++)
			{
				query.str("");
				query << "INSERT INTO META_Acquisition SET fid_AcquisitionInfo=" << acquisition_info_id << ",";
				query << "Number=" << info_it->getNumber();
				
				result = db_con_.executeQuery(query.str());
				parent_id = db_con_.getAutoId();
				
				storeMetaInfo_("META_Acquisition", parent_id, *info_it);
			}
		}
	}

	template <class ExperimentType>
	void DBAdapter::loadExperiment(UID id, ExperimentType& exp)
	{
		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		std::stringstream query; // query to build
		String tmp;              // temporary data
		QSqlQuery result, sub_result;        // place to store the query results in
		UID parent_id;					 // holds ID of parent data set
		
		query << "SELECT Type-1,Date,fid_MetaInfo,Description FROM META_MSExperiment WHERE id='" << id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		//Experiment meta info
		if (result.value(1).toDate().isValid())
		{
			DateTime d;
			d.set(result.value(1).toDateTime().toString(Qt::ISODate));
			exp.setDateTime(d);
		}
		exp.setComment(result.value(3).toString());
		loadMetaInfo_(result.value(2).toInt(),exp);
		
		std::vector<ProteinIdentification> pi_vec;
	  std::vector<ProteinHit> ph_vec;
		ProteinIdentification pi;
	  ProteinHit ph;
		
		query.str("");
		query << "SELECT id, SearchEngine, SearchEngineVersion, Date, ScoreType, HigherScoreBetter, SignificanceThreshold, fid_MetaInfo, fid_File	FROM ID_ProteinIdentification WHERE fid_MSExperiment='" << id << "'";

		result = db_con_.executeQuery(query.str());
		
		result.first();
		while(result.isValid())
		{
			parent_id = result.value(0).toInt();
			pi.setSearchEngine(result.value(1).toString());
			pi.setSearchEngineVersion(result.value(2).toString());
			pi.setDateTime(DateTime(result.value(3).toDateTime()));
			pi.setScoreType(result.value(4).toString());
			pi.setHigherScoreBetter(result.value(5).toInt());
			pi.setSignificanceThreshold(result.value(6).toDouble());
			
			loadMetaInfo_(result.value(7).toInt(), pi);
			//needs getter and setter methods first
			//			loadFile_(result.value(8).toInt(), pi.getSourceFile());
			
			query.str("");
			query << "SELECT Score, Accession, Sequence, fid_MetaInfo FROM ID_ProteinHit WHERE fid_ProteinIdentification='" << parent_id << "'";
			sub_result = db_con_.executeQuery(query.str());
		
			sub_result.first();
			while(sub_result.isValid())
			{
				ph.setScore(sub_result.value(0).toDouble());
				ph.setAccession(sub_result.value(1).toString());
				ph.setSequence(sub_result.value(2).toString());

				loadMetaInfo_(sub_result.value(3).toInt(), ph);

				ph_vec.push_back(ph);
				
				sub_result.next();
			}

			pi.setHits(ph_vec);

			pi_vec.push_back(pi);

			result.next();
		}
		
		exp.setProteinIdentifications(pi_vec);

		// Sample
		Sample sample;
		query.str("");
		// finding root of recursive sample tree
		query << "SELECT id FROM META_Sample WHERE fid_MSExperiment='" << id << "' AND fid_Sample IS NULL";
		result = db_con_.executeQuery(query.str());
		result.first();
		loadSample_ (result.value(0).toInt(), sample);
		exp.setSample(sample);
		
		// ContactPerson
		ContactPerson contact;
		query.str("");
		query << "SELECT PreName,LastName,Affiliation,Email,Comment,fid_MetaInfo FROM META_ContactPerson WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		result.first();
		while(result.isValid())
		{
			contact.setFirstName(result.value(0).toString());
			contact.setLastName(result.value(1).toString());
			contact.setInstitution(result.value(2).toString());
			contact.setEmail(result.value(3).toString());
			contact.setContactInfo(result.value(4).toString());
			loadMetaInfo_(result.value(5).toInt(),contact);
			result.next();
			exp.getContacts().push_back(contact);
		}
		
		// HPLC
		query.str("");
		query << "SELECT id,InstrumentName,ColumnName,Description,Flux,Pressure,Temperature FROM META_HPLC WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		parent_id=result.value(0).toInt();
		exp.getHPLC().setInstrument(result.value(1).toString());
		exp.getHPLC().setColumn(result.value(2).toString());
		exp.getHPLC().setComment(result.value(3).toString());
		exp.getHPLC().setFlux(result.value(4).toInt());
		exp.getHPLC().setPressure(result.value(5).toInt());
		exp.getHPLC().setTemperature(result.value(6).toInt());
		
		// Gradient*
		// I tried taking the big query apart in order to skip the double join, but this leads to
		// the problem of saving all requested keys in a vector in order to request the percentages (complex).
		// I'll still preserve the code in order to optimize in the future. Maybe I was just being blind. ;-)
		
		String last_name;
		bool timepoints_done = false;
		query.str("");
		/*
		query << "SELECT id,Name FROM META_GradientEluent WHERE fid_HPLC=" << parent_id;
		result = db_con_.executeQuery(query.str());
		result.first();
		while(result.isValid())
		{
			exp.getHPLC().getGradient().addEluent(result.value(0).toString());
			result.next();
		}
		
		query.str("");
		query << "SELECT id,Time FROM META_GradientTime WHERE fid_HPLC=" << parent_id;
		result = db_con_.executeQuery(query.str());
		result.first();
		while(result.isValid())
		{
			exp.getHPLC().getGradient().addTimepoint(result.value(0).toInt());
			result.next();
		}
		*/
		
		query << "SELECT Name,Time,Percentage FROM META_GradientEluent, META_GradientTime, META_GradientPercentage WHERE META_GradientEluent.fid_HPLC=" << parent_id << " AND fid_GradientEluent=META_GradientEluent.id AND fid_GradientTime=META_GradientTime.id";
		result = db_con_.executeQuery(query.str());		
		result.first();
		
		if (result.isValid())
		{
			last_name = result.value(0).toString();
			exp.getHPLC().getGradient().addEluent(last_name);
		}
		
		while(result.isValid())
		{
			if (result.value(0).toString() != last_name.toQString())
			{
				exp.getHPLC().getGradient().addEluent(result.value(0).toString());
				timepoints_done = true;
			}
			
			if (timepoints_done == false)
			{
				exp.getHPLC().getGradient().addTimepoint(result.value(1).toInt());
			}
			
			exp.getHPLC().getGradient().setPercentage(result.value(0).toString(), result.value(1).toInt(), result.value(2).toInt());
			
			last_name = result.value(0).toString();
			result.next();
		}
		
		// Instrument
		query.str("");
		query << "SELECT id,Model,Vendor,Description,fid_MetaInfo FROM META_MSInstrument WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		parent_id = result.value(0).toInt();
		exp.getInstrument().setModel(result.value(1).toString());
		exp.getInstrument().setVendor(result.value(2).toString());
		exp.getInstrument().setCustomizations(result.value(3).toString());
		loadMetaInfo_(result.value(4).toInt(),exp.getInstrument());
		
		// IonDetector
		query.str("");
		query << "SELECT AcquisitionMode-1,Type-1,Resolution,ADCSamplingFrequency,fid_MetaInfo FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		//TODO Hack: Multiple ion detectors should be supported
		exp.getInstrument().getIonDetectors().resize(1);
		exp.getInstrument().getIonDetectors()[0].setAcquisitionMode((IonDetector::AcquisitionMode) result.value(0).toInt());
		exp.getInstrument().getIonDetectors()[0].setType((IonDetector::Type) result.value(1).toInt());
		exp.getInstrument().getIonDetectors()[0].setResolution(result.value(2).toDouble());
		exp.getInstrument().getIonDetectors()[0].setADCSamplingFrequency(result.value(3).toDouble());
		loadMetaInfo_(result.value(4).toInt(),exp.getInstrument().getIonDetectors()[0]);
		
		// IonSource
		query.str("");
		query << "SELECT InletType-1,IonizationMethod-1,IonizationMode-1,fid_MetaInfo FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		//TODO Hack: Multiple ion sources should be supported
		exp.getInstrument().getIonSources().resize(1);
		exp.getInstrument().getIonSources()[0].setInletType((IonSource::InletType) result.value(0).toInt());
		exp.getInstrument().getIonSources()[0].setIonizationMethod((IonSource::IonizationMethod) result.value(1).toInt());
		exp.getInstrument().getIonSources()[0].setPolarity((IonSource::Polarity)(Int) result.value(2).toDouble());
		loadMetaInfo_(result.value(3).toInt(),exp.getInstrument().getIonSources()[0]);
		
		// MassAnalyzers
		MassAnalyzer analyzer;
		std::vector<MassAnalyzer> analyzers;
		query.str("");
		query << "SELECT Accuracy,FinalMSExponent,IsolationWidth,MagneticFieldStrength,ReflectronState-1,Resolution,ResolutionMethod-1,ResolutionType-1,ScanDirection-1,ScanLaw-1,ScanRate,ScanTime,TOFPathLength,Type-1,fid_MetaInfo FROM META_MassAnalyzer WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		
		result.first();
		while(result.isValid())
		{
			analyzer.setAccuracy(result.value(0).toDouble());
			analyzer.setFinalMSExponent(result.value(1).toInt());
			analyzer.setIsolationWidth(result.value(2).toDouble());
			analyzer.setMagneticFieldStrength(result.value(3).toDouble());
			analyzer.setReflectronState((MassAnalyzer::ReflectronState) result.value(4).toInt());
			analyzer.setResolution(result.value(5).toDouble());
			analyzer.setResolutionMethod((MassAnalyzer::ResolutionMethod) result.value(6).toInt());
			analyzer.setResolutionType((MassAnalyzer::ResolutionType) result.value(7).toInt());
			analyzer.setScanDirection((MassAnalyzer::ScanDirection) result.value(8).toInt());
			analyzer.setScanLaw((MassAnalyzer::ScanLaw) result.value(9).toInt());
			analyzer.setScanRate(result.value(10).toDouble());
			analyzer.setScanTime(result.value(11).toDouble());
			analyzer.setTOFTotalPathLength(result.value(12).toDouble());
			analyzer.setType((MassAnalyzer::AnalyzerType) result.value(13).toInt());
			loadMetaInfo_(result.value(14).toInt(), analyzer);
			
			analyzers.push_back(analyzer);
			result.next();
		}
		exp.getInstrument().setMassAnalyzers(analyzers);

		//id
		exp.setPersistenceId(id);

		// if we don't have to load the spectra, we're already done
		if (options_.getMetadataOnly())
		{
			return;
		}

		//spectra
		query.str("");
		query << "SELECT id FROM DATA_Spectrum WHERE fid_MSExperiment=" << id;
		if (options_.hasRTRange())
		{
			query << " AND RetentionTime > " << options_.getRTRange().min() << " AND RetentionTime < " << options_.getRTRange().max();
		}
		if (options_.hasMSLevels())
		{
			const std::vector<int>& levels = options_.getMSLevels();
			query << " AND (";
			for (std::vector<int>::const_iterator it = levels.begin(); it != levels.end(); it++)
			{
				query << "MSLevel=" << *it;
				if (it+1 != levels.end())
				{
					query << " OR ";
				}
			}
			query << ")";
		}
		query << " ORDER BY id ASC";
		
		result = db_con_.executeQuery(query.str());
		exp.resize(result.size());
		UInt i = 0;
		result.first();
		while (result.isValid())
		{
			loadSpectrum(result.value(0).toInt(), exp[i]);
			++i;
			result.next();
		}
	}

	template <class SpectrumType>
	void DBAdapter::loadSpectrum(UID id, SpectrumType& spec)
	{
		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION --------------------------------------- 
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;
		
		spec = SpectrumType();
		
		std::stringstream query;     // query to build
		QSqlQuery result, sub_result;// place to store the query results in
		InstrumentSettings settings; // stores settings that are read from DB
		UID parent_id;               // stores parent_id of Acquisition
		
		query << "SELECT Type-1,RetentionTime,MSLevel,Description,fid_MetaInfo,fid_File FROM DATA_Spectrum WHERE id='" << id << "'";
		result = db_con_.executeQuery(query.str());
		result.first();
		
		//Spectrum meta info
		spec.setType((SpectrumSettings::SpectrumType)(result.value(0).toInt()));
		spec.setRT(result.value(1).toDouble());
		spec.setMSLevel(result.value(2).toInt());
		spec.setComment(result.value(3).toString());
		loadMetaInfo_(result.value(4).toInt(),spec);
		loadFile_(result.value(5).toInt(),spec.getSourceFile());
		
		// Instrument settings
		query.str("");
		query << "SELECT MZRangeBegin, MZRangeEnd, Polarity-1, ScanMode-1, fid_MetaInfo FROM META_InstrumentSettings WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str());
		result.first();

		InstrumentSettings::ScanWindow window;
		window.begin = result.value(0).toDouble();
		window.end = result.value(1).toDouble();
		settings.getScanWindows().push_back(window);
		settings.setPolarity((IonSource::Polarity) (result.value(2).toInt()));
		settings.setScanMode((InstrumentSettings::ScanMode) (result.value(3).toInt()));
		spec.setInstrumentSettings(settings);
		loadMetaInfo_(result.value(4).toInt(),spec.getInstrumentSettings());


    // PeptideIdentification / PeptideHits

		std::vector<PeptideIdentification> pei_vec;
	  std::vector<PeptideHit> peh_vec;
		PeptideIdentification pei;
	  PeptideHit peh;
		
		query.str("");
		query << "SELECT id, SignificanceThreshold, ScoreType, HigherScoreBetter, fid_MetaInfo, fid_File	FROM ID_PeptideIdentification WHERE fid_Spectrum='" << id << "'";

		result = db_con_.executeQuery(query.str());
		
		result.first();
		while(result.isValid())
		{
			parent_id = result.value(0).toInt();
			pei.setSignificanceThreshold(result.value(1).toDouble());
			pei.setScoreType(result.value(2).toString());
			pei.setHigherScoreBetter(result.value(3).toInt());
			
			loadMetaInfo_(result.value(4).toInt(), pei);
			//needs getter and setter methods first
			//			loadFile_(result.value(5).toInt(), pei.getSourceFile());
			
			query.str("");
			query << "SELECT Score, Sequence, Charge, AABefore, AAAfter, fid_MetaInfo FROM ID_PeptideHit WHERE fid_Identification='" << parent_id << "'";
			sub_result = db_con_.executeQuery(query.str());
		
			sub_result.first();
			while(sub_result.isValid())
			{
				peh.setScore(sub_result.value(0).toDouble());
				peh.setSequence(String(sub_result.value(1).toString()));
				peh.setCharge(sub_result.value(2).toInt());
				peh.setAABefore(sub_result.value(3).toString().toStdString()[0]);
				peh.setAAAfter(sub_result.value(4).toString().toStdString()[0]);

				loadMetaInfo_(sub_result.value(5).toInt(), peh);

				peh_vec.push_back(peh);
				
				sub_result.next();
			}

			pei.setHits(peh_vec);

			pei_vec.push_back(pei);

			result.next();
		}

		spec.setPeptideIdentifications(pei_vec);

		// AcquisitionInfo
		query.str("");
		query << "SELECT id, MethodOfCombination FROM META_AcquisitionInfo WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str());
		result.first();
		
		spec.getAcquisitionInfo().setMethodOfCombination(result.value(1).toString());
		parent_id = result.value(0).toInt();
		
		// Acquisition
		query.str("");
		query << "SELECT Number,fid_MetaInfo FROM META_Acquisition WHERE fid_AcquisitionInfo='" << parent_id << "' ORDER BY id ASC";
		result = db_con_.executeQuery(query.str());

		Acquisition acquisition;
		
		result.first();
		while(result.isValid())
		{
			acquisition.setNumber(result.value(0).toInt());
			loadMetaInfo_(result.value(1).toInt(), acquisition);
			spec.getAcquisitionInfo().push_back(acquisition);
			result.next();
		}
		
		// MetaInfoDescription/MetaDataArrays
		query.str("");
		query << "SELECT Name, Description, fid_MetaInfo, fid_File ";
		query << "FROM META_MetaInfoDescription WHERE fid_Spectrum=" << id;

		result = db_con_.executeQuery(query.str());
		result.first();
		
		while(result.isValid())
		{
			typename SpectrumType::MetaDataArray meta_array;
			meta_array.setName(result.value(0).toString());
			meta_array.setComment(result.value(1).toString());
			loadMetaInfo_(result.value(2).toInt(), meta_array);
			loadFile_(result.value(3).toInt(),meta_array.getSourceFile());
			
			spec.getMetaDataArrays().push_back(meta_array);
			result.next();
		}

		//precursor
		if(spec.getMSLevel()>1)
		{
			query.str("");
			query << "SELECT mz,Intensity,Charge,ActivationMethod-1,ActivationEnergyUnit-1,ActivationEnergy,WindowSize,fid_MetaInfo FROM DATA_Precursor WHERE fid_Spectrum='" << id << "'";
			result = db_con_.executeQuery(query.str());
			result.first();
			spec.getPrecursorPeak().getPosition()[0] = (result.value(0).toDouble());
			spec.getPrecursorPeak().setIntensity(result.value(1).toDouble());
			spec.getPrecursorPeak().setCharge(result.value(2).toInt());
			spec.getPrecursor().setActivationMethod((Precursor::ActivationMethod)(result.value(3).toInt()));
			spec.getPrecursor().setActivationEnergyUnit((Precursor::EnergyUnits)(result.value(4).toInt()));
			spec.getPrecursor().setActivationEnergy(result.value(5).toDouble());
			spec.getPrecursor().setWindowSize(result.value(6).toDouble());
			loadMetaInfo_(result.value(7).toInt(),spec.getPrecursor());
		}
		
		// Peaks/MetaDataArrays
		query.str("");
		query << "SELECT mz,Intensity,fid_MetaInfo,id FROM DATA_Peak WHERE fid_Spectrum='" << id << "' ";
		if (options_.hasMZRange())
		{
			query << " AND mz > " << options_.getMZRange().min() << " AND mz < " << options_.getMZRange().max();
		}
		if (options_.hasIntensityRange())
		{
			query << " AND Intensity > " << options_.getIntensityRange().min() << " AND Intensity < " << options_.getIntensityRange().max();
		}
		query << " ORDER BY mz ASC";
		result = db_con_.executeQuery(query.str());

		typename SpectrumType::PeakType p;
		result.first();
		while(result.isValid())
		{
			p.setPosition(result.value(0).toDouble());
			p.setIntensity(result.value(1).toDouble());
			loadMetaInfo_(result.value(2).toInt(), p);
			spec.push_back(p);		
			for (typename SpectrumType::MetaDataArrays::iterator mdarrays_it = spec.getMetaDataArrays().begin(); mdarrays_it != spec.getMetaDataArrays().end(); mdarrays_it++)
			{
				query.str("");
				query << "SELECT id FROM META_MetaInfoDescription WHERE Name='";
				query << mdarrays_it->getName() << "' AND fid_Spectrum=" << id;
				sub_result = db_con_.executeQuery(query.str());
				sub_result.first();
				query.str("");
				query << "SELECT Value FROM DATA_PeakMetaData WHERE fid_Peak=";
				query << result.value(3).toInt() << " AND fid_MetaInfoDescription=" << sub_result.value(0).toInt();
				sub_result = db_con_.executeQuery(query.str());
				sub_result.first();
				mdarrays_it->push_back(sub_result.value(0).toDouble());
			}
			result.next();
		}
		
		//id
		spec.setPersistenceId(id);
	}
}

#endif
