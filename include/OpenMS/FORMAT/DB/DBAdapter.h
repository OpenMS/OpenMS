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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DB_DBADAPTER_H
#define OPENMS_FORMAT_DB_DBADAPTER_H

//OpenMS includes
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

//QT includes
#include <QtSql/QSqlQuery>
#include <QtCore/QVariant>
#include <QtCore/QDate>

//std and STL includes
#include <string>
#include <map>

namespace OpenMS
{
	class DBConnection;

  /**
  	@brief A class for accessing and storing data in a SQL database

    It can be used to create objects from the DB or store them in the DB.

		@todo Add DataProcessing to MetaInfoDescription (Hiwi, Mathias)
		@todo Check if test is really complete (Hiwi, Mathias)
		@todo Check that all values are quoted (Hiwi, Mathias)
		@todo Implement StringDataArray and IntegerDataArray of spectrum (Hiwi, Mathias)

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

		//----------------------------------------------------------------------------------------
		//------------------------------- CHECK DB VERSION ---------------------------------------
		//----------------------------------------------------------------------------------------
		if (!checkDBVersion(true)) return;



		//----------------------------------------------------------------------------------------
		//------------------------------- store EXPERIMENT ---------------------------------------
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
		//ExperimentalSettings FractionIdentifier
		query << ",FractionIdentifier='" << exp.getFractionIdentifier()<< "'";

		query << end;
		result = db_con_.executeQuery(query.str());
		if (new_entry)
		{
			exp.setPersistenceId(db_con_.getAutoId());
		}

		storeMetaInfo_("META_MSExperiment", exp.getPersistenceId(), exp);

		//----------------------------------------------------------------------------------------
		//----------- store	PROTEIN IDENTIFICATIONS/HITS / SEARCHPARAMETERS-----------------------
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

			//save the ProteinHits
			for (std::vector<ProteinHit>::const_iterator ph_it = pi_it->getHits().begin(); ph_it != pi_it->getHits().end(); ph_it++)
			{
				query.str("");
				query << "INSERT INTO ID_ProteinHit SET ";
				query << "fid_ProteinIdentification='" << parent_id << "'";
				query << ",Score='" << ph_it->getScore() << "'";
				query << ",Accession='" << ph_it->getAccession() << "'";
				query << ",Sequence='" << ph_it->getSequence() << "'";
				query << ",Rank='" << ph_it->getRank() << "'";

				result = db_con_.executeQuery(query.str());
				meta_id = db_con_.getAutoId();

				storeMetaInfo_("ID_ProteinHit", meta_id, *ph_it);
			}

			//save the searchparameters
			query.str("");
			query << "DELETE FROM ID_SearchParameters WHERE fid_ProteinIdentification='" << parent_id << "'";
			result = db_con_.executeQuery(query.str());

			query.str("");
			query << "INSERT INTO ID_SearchParameters SET ";
			query << "fid_ProteinIdentification='" << parent_id << "'";
			query << ",DB='" << pi_it->getSearchParameters().db << "'";
			query << ",DBVersion='" << pi_it->getSearchParameters().db_version << "'";
			query << ",Taxonomy='" << pi_it->getSearchParameters().taxonomy << "'";
			query << ",Charges='" << pi_it->getSearchParameters().charges << "'";
			query << ",MassType='" << (1u+pi_it->getSearchParameters().mass_type) << "'";
			query << ",Enzyme='" << (1u+pi_it->getSearchParameters().enzyme) << "'";
			query << ",MissedCleavages='" << pi_it->getSearchParameters().missed_cleavages << "'";
			query << ",PeakMassTolerance='" << pi_it->getSearchParameters().peak_mass_tolerance << "'";
			query << ",PrecursorTolerance='" << pi_it->getSearchParameters().precursor_tolerance << "'";

			result = db_con_.executeQuery(query.str());

			meta_id = db_con_.getAutoId();
			storeMetaInfo_("ID_SearchParameters", meta_id, pi_it->getSearchParameters());

			//if modifications then save them
			query.str("");
			query << "DELETE FROM ID_FixedModifications WHERE fid_SearchParameters='" << meta_id << "'";
			result = db_con_.executeQuery(query.str());

			for (std::vector<String>::const_iterator mod_it = pi_it->getSearchParameters().fixed_modifications.begin(); mod_it != pi_it->getSearchParameters().fixed_modifications.end(); mod_it++)
			{

				query.str("");
				query << "INSERT INTO ID_FixedModifications SET ";
				query << "fid_SearchParameters='" << meta_id << "'";
				query << ",name='" << *mod_it << "'";

				result = db_con_.executeQuery(query.str());
			}
			for (std::vector<String>::const_iterator mod_it = pi_it->getSearchParameters().variable_modifications.begin(); mod_it != pi_it->getSearchParameters().variable_modifications.end(); mod_it++)
			{
				query.str("");
				query << "INSERT INTO ID_VariableModifications SET ";
				query << "fid_SearchParameters='" << meta_id << "'";
				query << ",name='" << *mod_it << "'";

				result = db_con_.executeQuery(query.str());
			}
		}


		//----------------------------------------------------------------------------------------
		//-------------------------------- store SAMPLE ------------------------------------------
		//----------------------------------------------------------------------------------------

		query.str("");
		deleteMetaInfo_("META_Sample", "fid_MSExperiment=" + String(exp.getPersistenceId()));
		// this also deletes all references in META_SampleTreatment, META_Digestion and META_Modification by constraint
		query << "DELETE FROM META_Sample WHERE fid_MSExperiment='" << exp.getPersistenceId() << "'";
		storeSample_(exp.getSample(), exp.getPersistenceId(), 0);

		//----------------------------------------------------------------------------------------
		//-------------------------------- store CONTACTPERSON -----------------------------------
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
		//-------------------------------------- store HPLC --------------------------------------
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
			result = db_con_.executeQuery(query.str(),true);
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
		//------------------------- store GRADIENTEluent/Time/Percentage -------------------------
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
		//--------------------------------- store INSTRUMENT -------------------------------------
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
			result = db_con_.executeQuery(query.str(),true);
			parent_id = result.value(0).toInt();

			query.str("");
			query << "UPDATE META_MSInstrument SET ";
			end = " WHERE fid_MSExperiment='" + String (exp.getPersistenceId()) + "'";
		}

		query << "Model='" << instrument.getModel() << "'";
		query << ",Vendor='" << instrument.getVendor() << "'";
		query << ",Description='" << instrument.getCustomizations() << "'";
		query << ",IonOpticsType='" << (1u+instrument.getIonOptics()) << "'";

		query << end;
		result = db_con_.executeQuery(query.str());

		if (new_entry)
		{
			parent_id = db_con_.getAutoId();
		}

		storeMetaInfo_("META_MSInstrument", parent_id, instrument);

		deleteMetaInfo_("META_Software", "SoftwareApplicator='META_MSInstrument' AND fid_SoftwareApplicator=" + String(parent_id));
		query.str("");
		query << "DELETE FROM META_Software WHERE fid_SoftwareApplicator='" << parent_id << "' AND SoftwareApplicator='META_MSInstrument'";
		result = db_con_.executeQuery(query.str());
		query.str("");
		query << "INSERT INTO META_Software SET ";
		query << "fid_SoftwareApplicator='" << parent_id << "'";
		query << ",SoftwareApplicator='META_MSInstrument'";
		query << ",Name='" << instrument.getSoftware().getName()<< "'";
		query << ",Version='" << instrument.getSoftware().getVersion()<< "'";
		result = db_con_.executeQuery(query.str());

		UID software_id = db_con_.getAutoId();
		storeMetaInfo_("META_Software", software_id, instrument.getSoftware());

		//----------------------------------------------------------------------------------------
		//----------------------------------- store IONDETECTOR ----------------------------------
		//----------------------------------------------------------------------------------------
		const std::vector<IonDetector>& detectors = exp.getInstrument().getIonDetectors();
		query.str("");

		deleteMetaInfo_("META_IonDetector", "fid_MSInstrument=" + String(parent_id));
		query << "DELETE FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());

		for (std::vector<IonDetector>::const_iterator detectors_it = detectors.begin(); detectors_it != detectors.end(); detectors_it++)
		{
			query.str("");
			query << "INSERT INTO META_IonDetector SET ";
			query << "fid_MSInstrument='" << parent_id << "'";
			query << ",AcquisitionMode=" << (1u+detectors_it->getAcquisitionMode());
			query << ",Type=" << (1u+detectors_it->getType());
			query << ",Resolution=" << detectors_it->getResolution();
			query << ",ADCSamplingFrequency=" << detectors_it->getADCSamplingFrequency();
			query << ",InstrumentOrder=" << (detectors_it->getOrder());

			result = db_con_.executeQuery(query.str());
			storeMetaInfo_("META_IonDetector", db_con_.getAutoId(), *detectors_it);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------------- store IONSOURCE ----------------------------------
		//----------------------------------------------------------------------------------------

		const std::vector<IonSource>& sources = exp.getInstrument().getIonSources();
		query.str("");

		deleteMetaInfo_("META_IonSource", "fid_MSInstrument=" + String(parent_id));
		query << "DELETE FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());

		for (std::vector<IonSource>::const_iterator sources_it = sources.begin(); sources_it != sources.end(); sources_it++)
		{
			query.str("");
			query << "INSERT INTO META_IonSource SET ";
			query << "fid_MSInstrument='" << parent_id << "'";
			query << ",InletType=" << (1u+sources_it->getInletType());
			query << ",IonizationMethod=" << (1u+sources_it->getIonizationMethod());
			query << ",IonizationMode=" << (1u+sources_it->getPolarity());
			query << ",InstrumentOrder=" << (sources_it->getOrder());

			result = db_con_.executeQuery(query.str());
			storeMetaInfo_("META_IonSource", db_con_.getAutoId(), *sources_it);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------------ store MASSANALYZER --------------------------------
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
			query << ",InstrumentOrder=" << analyzer_it->getOrder();

			result = db_con_.executeQuery(query.str());
			storeMetaInfo_("META_MassAnalyzer", db_con_.getAutoId(), *analyzer_it);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------------- store SPECTRUM -----------------------------------
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
			//nativeID
			query << ",NativeID='" << exp_it->getNativeID() << "'";

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
			//----------------------------- store PEPTIDE IDENTIFICATIONS / HITS ---------------------
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
			//------------------------------------ store PRECURSOR -----------------------------------
			//----------------------------------------------------------------------------------------

			// first delete all old values, no matter whether we're updating or not
			// do we have to delete the MetaInfo as well?
			query.str("");
			query << "DELETE FROM DATA_Precursor WHERE fid_Spectrum='";
			query << exp_it->getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());


			for (Size precs_it = 0; precs_it < exp_it->getPrecursors().size(); precs_it++)
			{
				query.str("");
				query << "INSERT INTO DATA_Precursor SET ";
				query << "fid_Spectrum='"+ String(exp_it->getPersistenceId()) + "'";
				//Intensity
				query << ",Intensity='" << exp_it->getPrecursors()[precs_it].getIntensity() << "'";
				//mz
				query << ",WindowMz='" << exp_it->getPrecursors()[precs_it].getMZ() << "'";
				//charge
				query << ",Charge='" << exp_it->getPrecursors()[precs_it].getCharge() << "'";
				//activation energy
				query << ",ActivationEnergy='" << exp_it->getPrecursors()[precs_it].getActivationEnergy() << "'";
				//IsolationWindowLow
				query << ",WindowLow='" << exp_it->getPrecursors()[precs_it].getIsolationWindowLowerOffset() << "'";
				//IsolationWindowUp
				query << ",WindowUp='" << exp_it->getPrecursors()[precs_it].getIsolationWindowUpperOffset() << "'";

				result = db_con_.executeQuery(query.str());
				parent_id = db_con_.getAutoId();
				storeMetaInfo_("DATA_Precursor",parent_id, exp_it->getPrecursors()[precs_it]);


				for (Size pcs_it = 0; pcs_it < exp_it->getPrecursors()[precs_it].getPossibleChargeStates().size(); ++pcs_it)
				{
					query.str("");
					query << "INSERT INTO DATA_PrecursorPCS SET ";
					query << "fid_Precursor='"+ String(parent_id) + "'";
					//PossibleChargeStates
					query << ",PossibleChargeStates='" << (exp_it->getPrecursors()[precs_it].getPossibleChargeStates()[pcs_it]) << "'";
					result = db_con_.executeQuery(query.str());
				}

				for (std::set<Precursor::ActivationMethod>::iterator am_it = exp_it->getPrecursors()[precs_it].getActivationMethods().begin(); am_it != exp_it->getPrecursors()[precs_it].getActivationMethods().end(); ++am_it)
				{
					query.str("");
					query << "INSERT INTO DATA_PrecursorAM SET ";
					query << "fid_Precursor='"+ String(parent_id) + "'";
					//activation method
					query << ",ActivationMethods=" << (1u+*(am_it));
					result = db_con_.executeQuery(query.str());
				}

			}

			//----------------------------------------------------------------------------------------
			//------------------------------------- store PRODUCTS -----------------------------------
			//----------------------------------------------------------------------------------------

			// first delete all old values, no matter whether we're updating or not
			// do we have to delete the MetaInfo as well?
			query.str("");
			query << "DELETE FROM DATA_Products WHERE fid_Spectrum='";
			query << exp_it->getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());

			for (Size precs_it = 0; precs_it < exp_it->getProducts().size(); precs_it++)
			{
				query.str("");
				query << "INSERT INTO DATA_Products SET ";
				query << "fid_Spectrum='"+ String(exp_it->getPersistenceId()) + "'";
				//mz
				query << ",WindowMz='" << exp_it->getProducts()[precs_it].getMZ() << "'";
				//IsolationWindowLow
				query << ",WindowLow='" << exp_it->getProducts()[precs_it].getIsolationWindowLowerOffset() << "'";
				//IsolationWindowUp
				query << ",WindowUp='" << exp_it->getProducts()[precs_it].getIsolationWindowUpperOffset() << "'";

				result = db_con_.executeQuery(query.str());
				parent_id = db_con_.getAutoId();
				storeMetaInfo_("DATA_Products",parent_id, exp_it->getProducts()[precs_it]);
			}

			//----------------------------------------------------------------------------------------
			//--------------------------------------- store PEAKS ------------------------------------
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
			//-------------------- store METAINFODESCRIPTION / METADATAARRAYS  -----------------------
			//----------------------------------------------------------------------------------------

			const typename ExperimentType::SpectrumType::FloatDataArrays& meta_data_arrays = exp_it->getFloatDataArrays();

			for (typename ExperimentType::SpectrumType::FloatDataArrays::const_iterator mdarrays_it = meta_data_arrays.begin(); mdarrays_it != meta_data_arrays.end(); ++mdarrays_it)
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
					query << "Name='" + mdarrays_it->getName() + "' "; //TODO
					end  = " WHERE fid_Spectrum=" + String(exp_it->getPersistenceId());
					end += " AND Name='" + mdarrays_it->getName() + "'";
				}
				else
				{
					new_entry = true;
					query << "INSERT INTO META_MetaInfoDescription SET ";
					query << "fid_Spectrum=" << exp_it->getPersistenceId() << ", ";
					query << "Name='" << mdarrays_it->getName() << "'";
					end = "";
				}

				query << end;

				result = db_con_.executeQuery(query.str());
				if (new_entry)
				{
					parent_id = db_con_.getAutoId();
				}

				storeMetaInfo_("META_MetaInfoDescription", parent_id, *mdarrays_it);

				// store meta data contained in the FloatDataArrays
				query.str("");
				query << "DELETE FROM DATA_PeakMetaData WHERE fid_MetaInfoDescription=";
				query << parent_id;
				result = db_con_.executeQuery(query.str());

				query.str("");
				query << "SELECT id FROM DATA_Peak WHERE fid_Spectrum=" << exp_it->getPersistenceId();
				result = db_con_.executeQuery(query.str(),true);

				query.str("");
				query << "INSERT INTO DATA_PeakMetaData (fid_Peak,fid_MetaInfoDescription,Value) VALUES ";
				for (typename ExperimentType::SpectrumType::FloatDataArray::const_iterator meta_array_it = mdarrays_it->begin(); meta_array_it != mdarrays_it->end(); meta_array_it++)
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
			//------------------------------- store INSTRUMENT SETTINGS ------------------------------
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
				result = db_con_.executeQuery(query.str(),true);
				parent_id = result.value(0).toInt();

				query.str("");
				query << "UPDATE META_InstrumentSettings SET ";
				end = " WHERE fid_Spectrum='" + String(exp_it->getPersistenceId()) + "'";
			}

			query << "Polarity=" << (1u+settings.getPolarity()) << ",";
			query << "ScanMode=" << (1u+settings.getScanMode()) << ",";
			query << "ZoomScan=" << (settings.getZoomScan());
			query << end;

			result = db_con_.executeQuery(query.str());

			if (new_entry) parent_id = db_con_.getAutoId();
			storeMetaInfo_("META_InstrumentSettings", parent_id, exp_it->getInstrumentSettings());

			//----------------------------------------------------------------------------------------
			//----------------------------------- store SCANWINDOWS---- ------------------------------
			//----------------------------------------------------------------------------------------

			//handle several scan windows
			const std::vector<ScanWindow>& wins = settings.getScanWindows();

			// first delete all old values, no matter whether we're updating or not
			// do we have to delete the MetaInfo as well?
			query.str("");
			query << "DELETE FROM META_ScanWindows WHERE fid_Spectrum='";
			query << exp_it->getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());

			for (Size wins_it = 0; wins_it < wins.size(); wins_it++)
			{
				query.str("");
				query << "INSERT INTO META_ScanWindows SET ";
				query << "fid_Spectrum='"+ String(exp_it->getPersistenceId()) + "'";
				query << ",MZRangeBegin=" << settings.getScanWindows()[wins_it].begin;
				query << ",MZRangeEnd=" << settings.getScanWindows()[wins_it].end;
				result = db_con_.executeQuery(query.str());
				storeMetaInfo_("META_ScanWindows",parent_id, settings.getScanWindows()[wins_it]);
				//~ parent_id = db_con_.getAutoId();
			}

			//----------------------------------------------------------------------------------------
			//------------------------------- store ACQUISITIONINFO ----------------------------------
			//----------------------------------------------------------------------------------------

			const AcquisitionInfo & info = exp_it->getAcquisitionInfo();

			query.str("");

			if (new_entry)
			{
				query << "INSERT INTO META_AcquisitionInfo SET fid_Spectrum='" << exp_it->getPersistenceId() << "',";
				end = "";
			}
			else
			{
				query << "SELECT id FROM META_AcquisitionInfo WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				result = db_con_.executeQuery(query.str(),true);
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
			//-----------------------------------store ACQUISITION -----------------------------------
			//----------------------------------------------------------------------------------------

			query.str("");
			deleteMetaInfo_("META_Acquisition", "fid_AcquisitionInfo='" + String(parent_id) + "'");
			query << "DELETE FROM META_Acquisition WHERE fid_AcquisitionInfo='" << parent_id << "'";
			result = db_con_.executeQuery(query.str());

			for (std::vector<Acquisition>::const_iterator info_it = info.begin(); info_it != info.end(); info_it++)
			{
				query.str("");
				query << "INSERT INTO META_Acquisition SET fid_AcquisitionInfo='" << acquisition_info_id << "',";
				query << "Number='" << info_it->getIdentifier() << "'";

				result = db_con_.executeQuery(query.str());
				parent_id = db_con_.getAutoId();

				storeMetaInfo_("META_Acquisition", parent_id, *info_it);
			}


			//----------------------------------------------------------------------------------------
			//-------------------------------------store DATAPROCESSING ----------------------------------
			//----------------------------------------------------------------------------------------
			const std::vector<DataProcessing>& processings = exp_it->getDataProcessing();

			deleteMetaInfo_("META_DataProcessing", "fid_Spectrum=" + String(exp_it->getPersistenceId()));
			query.str("");
			query << "DELETE FROM META_DataProcessing WHERE fid_Spectrum='" << exp_it->getPersistenceId() << "'";
			result = db_con_.executeQuery(query.str());

			for (std::vector<DataProcessing>::const_iterator processings_it = processings.begin(); processings_it != processings.end(); processings_it++)
			{
				query.str("");
				query << "INSERT INTO META_DataProcessing SET ";
				query << "fid_Spectrum='" << exp_it->getPersistenceId() << "'";
				query << ",CompletionTime='" << processings_it->getCompletionTime().get()<< "'";

				result = db_con_.executeQuery(query.str());

				UID dataprocessing_id = db_con_.getAutoId();
				storeMetaInfo_("META_DataProcessing", dataprocessing_id, *processings_it);

				deleteMetaInfo_("META_Software", "SoftwareApplicator='META_DataProcessing' AND fid_SoftwareApplicator=" + String(dataprocessing_id));
				query.str("");
				query << "DELETE FROM META_Software WHERE fid_SoftwareApplicator='" << dataprocessing_id << "' AND SoftwareApplicator='META_DataProcessing'";
				result = db_con_.executeQuery(query.str());
				query.str("");
				query << "INSERT INTO META_Software SET ";
				query << "fid_SoftwareApplicator='" << dataprocessing_id << "'";
				query << ",SoftwareApplicator='META_DataProcessing'";
				query << ",Name='" << processings_it->getSoftware().getName()<< "'";
				query << ",Version='" << processings_it->getSoftware().getVersion()<< "'";
				result = db_con_.executeQuery(query.str());

				UID software_id = db_con_.getAutoId();
				storeMetaInfo_("META_Software", software_id, processings_it->getSoftware());

				//----------------------------------------------------------------------------------------
				//------------------------------------ store PROCESSINGACTIONS ---------------------------
				//----------------------------------------------------------------------------------------
				for (std::set<DataProcessing::ProcessingAction>::const_iterator acts_it = processings_it->getProcessingActions().begin(); acts_it != processings_it->getProcessingActions().end(); acts_it++)
				{
					query.str("");
					query << "INSERT INTO META_ProcessingActions SET ";
					query << "ProcessingActionType='" << (1u+(*acts_it)) << "'";
					query << ",fid_DataProcessing='" << dataprocessing_id << "'";
					result = db_con_.executeQuery(query.str());
				}
			}

		}//iteration over experiments' spectra ends here
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

		//----------------------------------------------------------------------------------------
		//------------------------------- load EXPERIMENTs ---------------------------------------
		//----------------------------------------------------------------------------------------

		query << "SELECT Date,fid_MetaInfo,Description,FractionIdentifier FROM META_MSExperiment WHERE id='" << id << "'";
		result = db_con_.executeQuery(query.str(), true);
		
		//Experiment meta info
		try
		{
			DateTime d;
			d.set(result.value(0).toDateTime().toString(Qt::ISODate));
			exp.setDateTime(d);
		}
		catch(Exception::ParseError& )
		{
			//no nothing, the date is simply unset
		}
		exp.setComment(result.value(2).toString());
		exp.setFractionIdentifier(result.value(3).toString());
		loadMetaInfo_(result.value(1).toInt(),exp);

		std::vector<ProteinIdentification> pi_vec;
	  std::vector<ProteinHit> ph_vec;
		ProteinIdentification pi;
	  ProteinHit ph;

		query.str("");
		query << "SELECT id, SearchEngine, SearchEngineVersion, Date, ScoreType, HigherScoreBetter, SignificanceThreshold, fid_MetaInfo, fid_File	FROM ID_ProteinIdentification WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str());
		while(result.next())
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

			//load searchparameters
			query.str("");
			query << "SELECT id,DB,DBVersion,Taxonomy,Charges,MassType-1,Enzyme-1,MissedCleavages,PeakMassTolerance,PrecursorTolerance,fid_MetaInfo FROM ID_SearchParameters WHERE fid_ProteinIdentification='" << parent_id << "'";
			sub_result = db_con_.executeQuery(query.str(),true);

			UID sub_id = sub_result.value(0).toInt();
			ProteinIdentification::SearchParameters params (pi.getSearchParameters());
			params.db = sub_result.value(1).toString();
			params.db_version = sub_result.value(2).toString();
			params.taxonomy = sub_result.value(3).toString();
			params.charges = sub_result.value(4).toString();
			params.mass_type = ((ProteinIdentification::PeakMassType)sub_result.value(5).toInt());
			params.enzyme = ((ProteinIdentification::DigestionEnzyme)sub_result.value(6).toInt());
			params.missed_cleavages = sub_result.value(7).toInt();
			params.peak_mass_tolerance = sub_result.value(8).toDouble();
			params.precursor_tolerance = sub_result.value(9).toDouble();
			loadMetaInfo_(sub_result.value(10).toInt(), params);

			//fill modifications
			query.str("");
			query << "SELECT name FROM ID_VariableModifications WHERE fid_SearchParameters='" << sub_id << "'";
			sub_result = db_con_.executeQuery(query.str());
			while(sub_result.next())
			{
				params.variable_modifications.push_back(sub_result.value(0).toString());
			}
			query.str("");
			query << "SELECT name FROM ID_FixedModifications WHERE fid_SearchParameters='" << sub_id << "'";
			sub_result = db_con_.executeQuery(query.str());
			while(sub_result.next())
			{
				params.fixed_modifications.push_back(sub_result.value(0).toString());
			}
			pi.setSearchParameters(params);


			query.str("");
			query << "SELECT Score, Accession, Sequence, Rank, fid_MetaInfo FROM ID_ProteinHit WHERE fid_ProteinIdentification='" << parent_id << "'";
			sub_result = db_con_.executeQuery(query.str());
			while(sub_result.next())
			{
				ph.setScore(sub_result.value(0).toDouble());
				ph.setAccession(sub_result.value(1).toString());
				ph.setSequence(sub_result.value(2).toString());
				ph.setRank(sub_result.value(3).toInt());

				loadMetaInfo_(sub_result.value(4).toInt(), ph);

				ph_vec.push_back(ph);
			}

			pi.setHits(ph_vec);

			pi_vec.push_back(pi);
		}

		exp.setProteinIdentifications(pi_vec);

		// Sample
		Sample sample;
		query.str("");
		// finding root of recursive sample tree
		query << "SELECT id FROM META_Sample WHERE fid_MSExperiment='" << id << "' AND fid_Sample IS NULL";
		result = db_con_.executeQuery(query.str(),true);
		loadSample_ (result.value(0).toInt(), sample);
		exp.setSample(sample);

		// ContactPerson
		ContactPerson contact;
		query.str("");
		query << "SELECT PreName,LastName,Affiliation,Email,Comment,fid_MetaInfo FROM META_ContactPerson WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			contact.setFirstName(result.value(0).toString());
			contact.setLastName(result.value(1).toString());
			contact.setInstitution(result.value(2).toString());
			contact.setEmail(result.value(3).toString());
			contact.setContactInfo(result.value(4).toString());
			loadMetaInfo_(result.value(5).toInt(),contact);
			exp.getContacts().push_back(contact);
		}

		// HPLC
		query.str("");
		query << "SELECT id,InstrumentName,ColumnName,Description,Flux,Pressure,Temperature FROM META_HPLC WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str(),true);
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
		while(result.next())
		{
			exp.getHPLC().getGradient().addEluent(result.value(0).toString());
		}

		query.str("");
		query << "SELECT id,Time FROM META_GradientTime WHERE fid_HPLC=" << parent_id;
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			exp.getHPLC().getGradient().addTimepoint(result.value(0).toInt());
		}
		*/

		query << "SELECT Name,Time,Percentage FROM META_GradientEluent, META_GradientTime, META_GradientPercentage WHERE META_GradientEluent.fid_HPLC=" << parent_id << " AND fid_GradientEluent=META_GradientEluent.id AND fid_GradientTime=META_GradientTime.id";
		result = db_con_.executeQuery(query.str(),true);
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

		//----------------------------------------------------------------------------------------
		//------------------------------- load INSTRUMENT	----------------------------------------
		//----------------------------------------------------------------------------------------

		query.str("");
		query << "SELECT id,Model,Vendor,Description,IonOpticsType-1,fid_MetaInfo FROM META_MSInstrument WHERE fid_MSExperiment='" << id << "'";
		result = db_con_.executeQuery(query.str(),true);

		parent_id = result.value(0).toInt();
		exp.getInstrument().setModel(result.value(1).toString());
		exp.getInstrument().setVendor(result.value(2).toString());
		exp.getInstrument().setCustomizations(result.value(3).toString());
		exp.getInstrument().setIonOptics((Instrument::IonOpticsType) result.value(4).toInt());
		loadMetaInfo_(result.value(5).toInt(),exp.getInstrument());

		query.str("");
		query << "SELECT Name,Version,fid_MetaInfo, id FROM META_Software WHERE fid_SoftwareApplicator='" << result.value(0).toInt() << "' AND SoftwareApplicator = 'META_MSInstrument'";
		result = db_con_.executeQuery(query.str(),true);
		if(result.isValid())
		{
			Software sw;
			sw.setName(result.value(0).toString());
			sw.setVersion(result.value(1).toString());
			loadMetaInfo_(result.value(2).toInt(),sw);

			exp.getInstrument().setSoftware(sw);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------- load IONDETECTORS	--------------------------------------
		//----------------------------------------------------------------------------------------
		std::vector<IonDetector> detectors;
		query.str("");
		query << "SELECT AcquisitionMode-1,Type-1,Resolution,ADCSamplingFrequency,InstrumentOrder,fid_MetaInfo FROM META_IonDetector WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			IonDetector detector;
			detector.setAcquisitionMode((IonDetector::AcquisitionMode) result.value(0).toInt());
			detector.setType((IonDetector::Type) result.value(1).toInt());
			detector.setResolution(result.value(2).toDouble());
			detector.setADCSamplingFrequency(result.value(3).toDouble());
			detector.setOrder(result.value(4).toInt());
			loadMetaInfo_(result.value(5).toInt(),detector);

			detectors.push_back(detector);
		}
		exp.getInstrument().setIonDetectors(detectors);

		//----------------------------------------------------------------------------------------
		//------------------------------- load IONSOURCES	--------------------------------------
		//----------------------------------------------------------------------------------------
		std::vector<IonSource> sources;
		query.str("");
		query << "SELECT InletType-1,IonizationMethod-1,IonizationMode-1,InstrumentOrder,fid_MetaInfo FROM META_IonSource WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			IonSource source;
			source.setInletType((IonSource::InletType) result.value(0).toInt());
			source.setIonizationMethod((IonSource::IonizationMethod) result.value(1).toInt());
			source.setPolarity((IonSource::Polarity)(Int) result.value(2).toDouble());
			source.setOrder(result.value(3).toInt());
			loadMetaInfo_(result.value(4).toInt(),source);

			sources.push_back(source);
		}
		exp.getInstrument().setIonSources(sources);

		//----------------------------------------------------------------------------------------
		//------------------------------- load MASSANALYZERS	------------------------------------
		//----------------------------------------------------------------------------------------
		std::vector<MassAnalyzer> analyzers;
		query.str("");
		query << "SELECT Accuracy,FinalMSExponent,IsolationWidth,MagneticFieldStrength,ReflectronState-1,Resolution,ResolutionMethod-1,ResolutionType-1,ScanDirection-1,ScanLaw-1,ScanRate,ScanTime,TOFPathLength,Type-1,InstrumentOrder,fid_MetaInfo FROM META_MassAnalyzer WHERE fid_MSInstrument='" << parent_id << "'";
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			MassAnalyzer analyzer;
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
			analyzer.setOrder(result.value(14).toInt());
			loadMetaInfo_(result.value(15).toInt(), analyzer);

			analyzers.push_back(analyzer);
		}
		exp.getInstrument().setMassAnalyzers(analyzers);

		//id
		exp.setPersistenceId(id);

		// if we don't have to load the spectra, we're already done
		if (options_.getMetadataOnly())
		{
			return;
		}

		//----------------------------------------------------------------------------------------
		//----------------------------------- load SPECTRAS	--------------------------------------
		//----------------------------------------------------------------------------------------
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
		while (result.next())
		{
			loadSpectrum(result.value(0).toInt(), exp[i]);
			++i;
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

		query << "SELECT Type-1,NativeID, RetentionTime,MSLevel,Description,fid_MetaInfo,fid_File FROM DATA_Spectrum WHERE id='" << id << "'";
		result = db_con_.executeQuery(query.str(),true);

		//Spectrum meta info
		spec.setType((SpectrumSettings::SpectrumType)(result.value(0).toInt()));
		spec.setNativeID(result.value(1).toString());
		spec.setRT(result.value(2).toDouble());
		spec.setMSLevel(result.value(3).toInt());
		spec.setComment(result.value(4).toString());
		loadMetaInfo_(result.value(5).toInt(),spec);
		loadFile_(result.value(6).toInt(),spec.getSourceFile());

		//----------------------------------------------------------------------------------------
		//------------------------------- load INSTRUMENTSETTINGS	--------------------------------
		//----------------------------------------------------------------------------------------
		query.str("");
		query << "SELECT Polarity-1, ScanMode-1, ZoomScan, fid_MetaInfo FROM META_InstrumentSettings WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str(),true);

		settings.setPolarity((IonSource::Polarity) (result.value(0).toInt()));
		settings.setScanMode((InstrumentSettings::ScanMode) (result.value(1).toInt()));
		settings.setZoomScan(result.value(2).toBool());
		spec.setInstrumentSettings(settings);
		loadMetaInfo_(result.value(3).toInt(),spec.getInstrumentSettings());

		//----------------------------------------------------------------------------------------
		//------------------------------- load SCANWINDOWS	--------------------------------------
		//----------------------------------------------------------------------------------------
		query.str("");
		query << "SELECT MZRangeBegin,MZRangeEnd,fid_MetaInfo FROM META_ScanWindows WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			ScanWindow window;
			window.begin = result.value(0).toDouble();
			window.end = result.value(1).toDouble();
			loadMetaInfo_(result.value(2).toInt(),window);
			spec.getInstrumentSettings().getScanWindows().push_back(window);
		}

		//----------------------------------------------------------------------------------------
		//----------------- load PEPTIDEIDENTIFICATION / PEPTIDEHITS	----------------------------
		//----------------------------------------------------------------------------------------


		std::vector<PeptideIdentification> pei_vec;
	  std::vector<PeptideHit> peh_vec;
		PeptideIdentification pei;
	  PeptideHit peh;

		query.str("");
		query << "SELECT id, SignificanceThreshold, ScoreType, HigherScoreBetter, fid_MetaInfo, fid_File	FROM ID_PeptideIdentification WHERE fid_Spectrum='" << id << "'";

		result = db_con_.executeQuery(query.str());
		while(result.next())
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

			while(sub_result.next())
			{
				peh.setScore(sub_result.value(0).toDouble());
				peh.setSequence(String(sub_result.value(1).toString()));
				peh.setCharge(sub_result.value(2).toInt());
				peh.setAABefore(sub_result.value(3).toString().toStdString()[0]);
				peh.setAAAfter(sub_result.value(4).toString().toStdString()[0]);

				loadMetaInfo_(sub_result.value(5).toInt(), peh);

				peh_vec.push_back(peh);
			}

			pei.setHits(peh_vec);

			pei_vec.push_back(pei);
		}

		spec.setPeptideIdentifications(pei_vec);

		//----------------------------------------------------------------------------------------
		//------------------------------- load AQUISITIONINFO	--------------------------------------
		//----------------------------------------------------------------------------------------

		query.str("");
		query << "SELECT id, MethodOfCombination FROM META_AcquisitionInfo WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str(),true);

		spec.getAcquisitionInfo().setMethodOfCombination(result.value(1).toString());
		parent_id = result.value(0).toInt();

		//----------------------------------------------------------------------------------------
		//------------------------------- load AQUISITION ----------------------------------------
		//----------------------------------------------------------------------------------------
		query.str("");
		query << "SELECT Number,fid_MetaInfo FROM META_Acquisition WHERE fid_AcquisitionInfo='" << parent_id << "' ORDER BY id ASC";
		result = db_con_.executeQuery(query.str());

		while(result.next())
		{
			Acquisition acquisition;
			acquisition.setIdentifier(result.value(0).toString());
			loadMetaInfo_(result.value(1).toInt(), acquisition);
			spec.getAcquisitionInfo().push_back(acquisition);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------- load DATAPROCESSINGS -----------------------------------
		//----------------------------------------------------------------------------------------
		query.str("");
		query << "SELECT CompletionTime,fid_MetaInfo, id FROM META_DataProcessing WHERE fid_Spectrum='" << id << "'";
		result = db_con_.executeQuery(query.str());

		while(result.next())
		{
			DataProcessing processings;

			try
			{
				DateTime d;
				d.set(result.value(0).toDateTime().toString(Qt::ISODate));
				processings.setCompletionTime(d);
			}
			catch(Exception::ParseError& )
			{
				//no nothing, the date is simply unset
			}

			loadMetaInfo_(result.value(1).toInt(),processings);

			query.str("");
			query << "SELECT ProcessingActionType-1 FROM META_ProcessingActions WHERE fid_DataProcessing='" << result.value(2).toInt() << "'";
			sub_result = db_con_.executeQuery(query.str());
			while(sub_result.next())
			{
				processings.getProcessingActions().insert((DataProcessing::ProcessingAction)sub_result.value(0).toInt());
			}

			query.str("");
			query << "SELECT Name,Version,fid_MetaInfo, id FROM META_Software WHERE fid_SoftwareApplicator='" << result.value(2).toInt() << "' AND SoftwareApplicator = 'META_DataProcessing'";
			sub_result = db_con_.executeQuery(query.str(),true);
			if(sub_result.isValid())
			{
				Software sw;
				sw.setName(sub_result.value(0).toString());
				sw.setVersion(sub_result.value(1).toString());
				loadMetaInfo_(sub_result.value(2).toInt(),sw);

				processings.setSoftware(sw);
			}

			spec.getDataProcessing().push_back(processings);
		}

		//----------------------------------------------------------------------------------------
		//----------------- load METAINFODESCRIPTION/METADATAARRAYS	------------------------------
		//----------------------------------------------------------------------------------------

		query.str("");
		query << "SELECT Name, fid_MetaInfo FROM META_MetaInfoDescription WHERE fid_Spectrum=" << id;
		result = db_con_.executeQuery(query.str());

		while(result.next())
		{
			typename SpectrumType::FloatDataArray meta_array;
			meta_array.setName(result.value(0).toString());
			loadMetaInfo_(result.value(1).toInt(), meta_array);
			
			spec.getFloatDataArrays().push_back(meta_array);
		}

		//----------------------------------------------------------------------------------------
		//------------------------------- load PRECURSOR	----------------------------------------
		//----------------------------------------------------------------------------------------
		if(spec.getMSLevel()>1)
		{
			query.str("");
			query << "SELECT WindowMz,Intensity,Charge,ActivationEnergy,WindowLow,WindowUp,fid_MetaInfo,id FROM DATA_Precursor WHERE fid_Spectrum='" << id << "' HAVING Intensity IS NOT NULL";
			result = db_con_.executeQuery(query.str());
			spec.getPrecursors().resize(result.size());
			UInt res = 0;
			while(result.next())
			{
				spec.getPrecursors()[res].setMZ(result.value(0).toDouble());
				spec.getPrecursors()[res].setIntensity(result.value(1).toDouble());
				spec.getPrecursors()[res].setCharge(result.value(2).toInt());
				spec.getPrecursors()[res].setActivationEnergy(result.value(3).toDouble());
				spec.getPrecursors()[res].setIsolationWindowLowerOffset(result.value(4).toDouble());
				spec.getPrecursors()[res].setIsolationWindowUpperOffset(result.value(5).toDouble());
				loadMetaInfo_(result.value(6).toInt(),spec.getPrecursors()[res]);

				UID prec_id = result.value(7).toInt();
				query.str("");
				query << "SELECT PossibleChargeStates FROM DATA_PrecursorPCS WHERE fid_Precursor='" << prec_id << "' HAVING PossibleChargeStates IS NOT NULL";
				QSqlQuery subresult = db_con_.executeQuery(query.str());
				while(subresult.next())
				{
					spec.getPrecursors()[res].getPossibleChargeStates().push_back(subresult.value(0).toInt());
				}

				query.str("");
				query << "SELECT ActivationMethods-1 FROM DATA_PrecursorAM WHERE fid_Precursor='" << prec_id << "'";
				subresult = db_con_.executeQuery(query.str());
				std::set<Precursor::ActivationMethod> tmp_set;
				while(subresult.next())
				{
					tmp_set.insert((Precursor::ActivationMethod)(subresult.value(0).toInt()));
				}
				spec.getPrecursors()[res].setActivationMethods(tmp_set);
				++res;
			}
		}

		//----------------------------------------------------------------------------------------
		//------------------------------ load PRODUCTS	------------------------------------------
		//----------------------------------------------------------------------------------------
		query.str("");
		query << "SELECT WindowMz,WindowLow,WindowUp,fid_MetaInfo FROM DATA_Products WHERE fid_Spectrum='" << id << "'";
		result = db_con_.executeQuery(query.str());
		spec.getProducts().resize(result.size());
		UInt res=0;
		while(result.next())
		{
			spec.getProducts()[res].setMZ(result.value(0).toDouble());
			spec.getProducts()[res].setIsolationWindowLowerOffset(result.value(1).toDouble());
			spec.getProducts()[res].setIsolationWindowUpperOffset(result.value(2).toDouble());
			loadMetaInfo_(result.value(3).toInt(),spec.getProducts()[res]);
			++res;
		}
		
		//----------------------------------------------------------------------------------------
		//--------------------------- load PEAKS/METADATAARRAYS	----------------------------------
		//----------------------------------------------------------------------------------------

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

		while(result.next())
		{
			typename SpectrumType::PeakType p;
			p.setPosition(result.value(0).toDouble());
			p.setIntensity(result.value(1).toDouble());
			loadMetaInfo_(result.value(2).toInt(), p);
			spec.push_back(p);
			for (typename SpectrumType::FloatDataArrays::iterator mdarrays_it = spec.getFloatDataArrays().begin(); mdarrays_it != spec.getFloatDataArrays().end(); mdarrays_it++)
			{
				query.str("");
				query << "SELECT id FROM META_MetaInfoDescription WHERE Name='";
				query << mdarrays_it->getName() << "' AND fid_Spectrum=" << id;
				sub_result = db_con_.executeQuery(query.str(),true);
				query.str("");
				query << "SELECT Value FROM DATA_PeakMetaData WHERE fid_Peak=";
				query << result.value(3).toInt() << " AND fid_MetaInfoDescription=" << sub_result.value(0).toInt();
				sub_result = db_con_.executeQuery(query.str(),true);
				mdarrays_it->push_back(sub_result.value(0).toDouble());
			}
		}

		//id
		spec.setPersistenceId(id);
	}
}

#endif
