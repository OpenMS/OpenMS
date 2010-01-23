// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

/*

Implementation policy:

To keep the database in a consistent and clean state, the following policies are required to be fulfilled:
- As soon as an element is not referenced anymore (i.e. "not reachable by any path from any root element"),
  it must be removed from the database.
- Root elements of reference are:
  - a spectrum (loadSpectrum)
  - an experiment (loadExperiment / storeExperiment)

Problem cases appear in the following scenarios:
- an element does not need a certain property (optional properties)
- an element can have more than one property of a kind (1:n properties)

These cases are problematic because the implementation has to deal with the following procedure whenever an element is stored to the database:

Distinguish whether the entry is mandatory (1:1, ALWAYS EXACTLY 1 ENTRY) or optional (1:n, n entries)

If 1:1:
- if the element has never been stored before, store it.
- if the element has been stored before, overwrite it

If 1:n:
- delete all outdated elements AND THEIR CHILDREN and create current elements AND CHILDREN from scratch
  (otherwise we could run into big trouble if n changes)
  	
watch out: if an element, which is referenced in a 1:1 fashion, is not set, consider it mandatory anyway and store the empty element. This prevents a lot of case checking and makes implementation much more straightforward.


What do we do in the current implementation?

- storeExperiment:
	META_MSExperiment: if exp.getPersistenceId()==0: INSERT, else UPDATE
											storeMetaInfo_("META_MSExperiment", exp.getPersistenceId(), exp);
	DATA_Spectrum:		 if exp_it->getPersistenceId()==0: INSERT, else UPDATE
											storeFile_("DATA_Spectrum", exp_it->getPersistenceId(), exp_it->getSourceFile());
											storeMetaInfo_("DATA_Spectrum", exp_it->getPersistenceId(), *exp_it);
	META_MetaInfoDescription:
		foreach (exp_it->getMetaInfoDescriptions())
			if exists "SELECT id FROM META_MetaInfoDescription WHERE fid_Spectrum=exp_it->getPersistenceId() AND Name=desc_it->first: UPDATE; else INSERT
			storeFile_("META_MetaInfoDescription", parent_id, desc_it->second.getSourceFile());
			storeMetaInfo_("META_MetaInfoDescription", parent_id, desc_it->second);
	DATA_Precursor:		 if new_spectrum: INSERT, else UPDATE
											storeMetaInfo_("DATA_Precursor",parent_id, exp_it->getPrecursors()[0].);
	DATA_Peak: delete Meta, DELETE ALL; INSERT ALL; storeMetaInfo_;
	META_InstrumentSettings: if new_spectrum INSERT; else UPDATE;
														storeMetaInfo_("META_InstrumentSettings", parent_id, exp_it->getInstrumentSettings());
	META_AcquisitionInfo: if new_spectrum INSERT; else UPDATE;
	META_Acquisition: delete Meta, DELETE ALL; INSERT ALL;
	META_Instrument: if new_experiment INSERT, else UPDATE
	META_IonSource: if new_experiment INSERT, else UPDATE
	META_IonDetector: if new_experiment INSERT, else UPDATE
	META_MassAnalyzer: delete Meta, DELETE ALL, INSERT ALL;
  META_ContactPerson: delete Meta, DELETE ALL, INSERT ALL;
  META_HPLC: if new_experiment INSERT, else UPDATE
	
*/

//OpenMS includes
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>

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
		result = db_con_.executeQuery(query.str());
		while(result.next())
		{
			switch(result.value(0).toInt())
			{
				case 0: //string
					info.setMetaValue(result.value(1).toString(),result.value(2).toString());
					break;
				case 1: //double
					info.setMetaValue(result.value(1).toString(),result.value(2).toDouble());
					break;
				case 2: //int
					info.setMetaValue(result.value(1).toString(),result.value(2).toInt());
					break;
				default:
					throw Exception::BaseException(__FILE__, __LINE__, __PRETTY_FUNCTION__,"DBAdapter","Unknown META_TypeNameValue:type in DBAdapter!");
			}
		}
	}
	
	UID DBAdapter::storeMetaInfo_(const String& parent_table, UID parent_id, const MetaInfoInterface& info)
	{
		/* All MetaInfo data in META_TypeNameValue needs to be connected to entries in parent tables in 1:n fashion.
		   The additional table META_MetaInfo forces links between META_TypeNameValue and the parent table to get
		     a globally unique meta ID.
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
		bool debug=false;
		
		query << "SELECT fid_MetaInfo FROM " << parent_table << " WHERE id='" << parent_id << "' AND fid_MetaInfo IS NOT NULL";
		result = db_con_.executeQuery(query.str(),true);
		
		//metainfo present => delete values of it
		UID meta_id = 0;
		if (result.size()==1)
		{
			// information exists in DB -> clearing content to replace with more recent information
			if (debug) cout << "MetaInfo for entry '" << parent_id << "' in table '" << parent_table << "' already exists => deleting TypeNameValues..." << endl;
			meta_id = result.value(0).toInt();
			query.str("");
			query << "DELETE FROM META_TypeNameValue WHERE fid_MetaInfo='" << meta_id << "'";
			result = db_con_.executeQuery(query.str());
		}
		
		//connection between metainfo and object
		if (info.isMetaEmpty() && meta_id!=0)
		{
			// information exists in DB, but new information is empty -> reference to old information needs to be deleted
			if (debug) cout << "Nothing to save for entry '" << parent_id << "' in table '" << parent_table << "' => clearing reference..." << endl;
			query.str("");
			query << "DELETE FROM META_MetaInfo WHERE id=" << meta_id;
			result = db_con_.executeQuery(query.str());
			query.str("");
			query << "UPDATE " << parent_table << " SET fid_MetaInfo=NULL WHERE id=" << parent_id;
			result = db_con_.executeQuery(query.str());
		}
		
		if ((!info.isMetaEmpty()) && meta_id==0)
		{
			// information does not exist in DB, but there is new information to save -> create reference
			result = db_con_.executeQuery("INSERT INTO META_MetaInfo () VALUES ()");
			meta_id = db_con_.getAutoId();
			if (debug) cout << "Unsaved MetaInfo to save, saving " << meta_id << " to entry '" << parent_id << "' in table '" << parent_table << "' => creating reference..." << endl;
			query.str("");
			query << "UPDATE " << parent_table << " SET fid_MetaInfo=" << meta_id << " WHERE id=" << parent_id;
			if (debug) cout << query.str() << endl;
			result = db_con_.executeQuery(query.str());
		}

		if (!info.isMetaEmpty())
		{
			// there is new information to save -> create contents
			if (debug) cout << "writing meta values..." << endl;
			
			query.str("");
			query << "INSERT INTO META_TypeNameValue (fid_MetaInfo,Name,Type,Value) VALUES ";
			vector<String> keys;
			info.getKeys(keys);
			const DataValue* val = 0;
			vector<String>::iterator it = keys.begin();
			while (true)
			{
				query << "('" << meta_id << "','" << *it;
				val = &info.getMetaValue(*it);
				if (debug) cout << *it << "=" << *val << endl;
				//TODO Handling of int list, double list and string list
				switch (val->valueType())
				{
					case DataValue::STRING_VALUE:
						query << "','string','" << string(*val);
						break;
					case DataValue::INT_VALUE:
						query << "','int','" << Int(*val);
						break;
					case DataValue::DOUBLE_VALUE:
						query << "','double','" << DoubleReal(*val);
						break;
					default:
						throw Exception::BaseException(__FILE__, __LINE__, __PRETTY_FUNCTION__,"DBAdapter","Unknown DataValue type in DBAdapter!");
				}
				query << "')";
				++it;
				if (it==keys.end())
				{
					break;
				}
				query << ",";
			}
			result = db_con_.executeQuery(query.str());
		}
		
		return meta_id;
	}

	UID DBAdapter::storeMetaInfo_(const String& /*parent_table*/, UID /*parent_id*/, const Peak1D& /*peak*/)
	{
		return 0;
	}

	UID DBAdapter::storeMetaInfo_(const String& parent_table, UID parent_id, const RichPeak1D& peak)
	{
		return storeMetaInfo_(parent_table, parent_id, static_cast<const MetaInfoInterface&>(peak));
	}

  UID DBAdapter::storeMetaInfo_(const String& parent_table, UID parent_id, const Precursor& peak)
  {
    return storeMetaInfo_(parent_table, parent_id, static_cast<const MetaInfoInterface&>(peak));
  }
	
	void DBAdapter::deleteMetaInfo_(const String& parent_table, const String& condition)
	{
		//cout << "Loading Meta: " << id << endl;
		QSqlQuery result_select, result_delete;
		stringstream query;
		
		query.str("");
		query << "SELECT fid_MetaInfo FROM " << parent_table << " WHERE " << condition;
		result_select = db_con_.executeQuery(query.str());
		
		while (result_select.next())
		{
			query.str("");
			query << "DELETE FROM META_TypeNameValue WHERE fid_MetaInfo='" << result_select.value(0).toInt() << "'";
			result_delete = db_con_.executeQuery(query.str(),true);
			query.str("");
			query << "DELETE FROM META_MetaInfo WHERE id='" << result_select.value(0).toInt() << "'";
			result_delete = db_con_.executeQuery(query.str());
		}
	}
	
	void DBAdapter::loadFile_(UID id, SourceFile& file)
	{
		// if there is no SourceFile to load (=NULL in DB), set empty SourceFile and return
		if (id == 0)
		{
			SourceFile empty_sourcefile;
			file = empty_sourcefile;
			return;
		}
		
		stringstream query;
		query << "SELECT FileName, FilePath, Size, `Type`, sha1, ChecksumType-1,NativeIDType FROM META_File WHERE id='" << id << "'";
		QSqlQuery result = db_con_.executeQuery(query.str(),true);

		file.setNameOfFile(result.value(0).toString());
		file.setPathToFile(result.value(1).toString());
		file.setFileSize(result.value(2).toDouble());
		file.setFileType(result.value(3).toString());
		file.setChecksum(result.value(4).toString(), (SourceFile::ChecksumType)result.value(5).toInt());
		file.setNativeIDType(result.value(6).toString());
	}

	UID DBAdapter::storeFile_(const String& parent_table, UID parent_id, const SourceFile& file)
	{
		QSqlQuery result;
		stringstream query;
		String end;
		bool new_entry = true;
		bool debug = false;
		UID file_id = 0;
		
//		// If there is no file to save, set empty reference and return
//		if (file.isFileEmpty())
//		{
//			if (debug) cout << "Empty file for " << parent_table << " given, skipping..." << endl;
//			query.str("");
//			query << "UPDATE " << parent_table << " SET fid_File=NULL";
//			query << " WHERE id=" << String(parent_id);
//			result = db_con_.executeQuery(query.str());
//			return 0;
//		}
		
		if (debug) cout << "file given, saving to entry '" << parent_id << "' in table '" << parent_table << endl;
		query.str("");
		query << "SELECT fid_File FROM " << parent_table << " WHERE id=";
		query << String(parent_id) << " AND fid_File IS NOT NULL";
		if (debug) cout << query.str() << endl;
		result = db_con_.executeQuery(query.str(),true);
		query.str("");
		
		if (result.size() > 0)				// reference already exists
		{
			file_id = result.value(0).toInt();
			new_entry = false;
		}
				
		if (new_entry)								// no file information stored yet, storing new entry
		{
			query << "INSERT INTO META_File SET ";
			end = "";
		}
		else													// we already stored some file information, updating old entry
		{
			query << "UPDATE META_File SET ";
			end = " WHERE id=" + String(parent_id);
		}
		
		query << "FileName='" << file.getNameOfFile() << "',";
		query << "FilePath='" << file.getPathToFile() << "',";
		query << "Size=" << file.getFileSize() << ",";
		query << "sha1='" << file.getChecksum() << "',";
		query << "`Type`='" << file.getFileType() << "'";
		query << ",ChecksumType='" << 1u+file.getChecksumType() << "'";
		query << ",NativeIDType='" << file.getNativeIDType() << "'";
		query << end;

		if (debug) cout << query.str() << endl;
		result = db_con_.executeQuery(query.str());
		query.str("");
		
		if (new_entry)
		{
			file_id = db_con_.getAutoId();
			query << "UPDATE " << parent_table << " SET fid_File=";
			query << String(file_id) << " WHERE id=" << String(parent_id);
			if (debug) cout << query.str() << endl;
			result = db_con_.executeQuery(query.str());
		}

		return file_id;
	}

	UID DBAdapter::storeSample_(const Sample& sample, UID exp_id, UID parent_id)
	{
		//----------------------------------------------------------------------------------------
		//-------------------------------------- SAMPLE ------------------------------------------
		//----------------------------------------------------------------------------------------
		
		/*
		we don't know if this is the first call or if we're recursing.
		each subsample pwns foreign keys to
		- the supersample
		- the experiment
		
		fid_Sample points to the *parent* sample, so we need to store the parent first - otherwise we
		don't know the parent id.
		
		*/
		
		QSqlQuery result;
		stringstream query, treatment_query;
		UID meta_id;
		
		query.str("");		
		query << "INSERT INTO META_Sample SET ";
		query << "fid_MSExperiment=" << exp_id;
		if (parent_id > 0) query << ",fid_Sample=" << parent_id;
		query << ",Name='" << sample.getName() << "'";
		query << ",SampleID='" << sample.getNumber() << "'";
		query << ",Mass='" << sample.getMass() << "'";
		query << ",Volume='" << sample.getVolume() << "'";
		query << ",Concentration='" << sample.getConcentration() << "'";
		query << ",State=" << (1u+sample.getState());
		query << ",Organism='" << sample.getOrganism() << "'";
		query << ",Description='" << sample.getComment() << "'";
		result = db_con_.executeQuery(query.str());
		
		parent_id = db_con_.getAutoId();

		storeMetaInfo_("META_Sample", parent_id, sample);
		
		// recreate all treatments
		// a treatment can be a digestion or a modification
		const Digestion* digestion;
		const Modification* modification;
		const Tagging* tagging;
		
		for (int i = 0; i != sample.countTreatments(); i++)
		{
			// first save digestion / modification to get parent_ids for foreign keys in META_SampleTreatment
			
			query.str("");		
			treatment_query.str("");		
			query << "INSERT INTO META_SampleTreatment SET ";			
			query << "fid_Sample=" << parent_id;
			if (sample.getTreatment(i).getType() == "Digestion")
			{
				digestion = dynamic_cast<const Digestion*>(&sample.getTreatment(i));
				treatment_query << "INSERT INTO META_Digestion SET ";
				treatment_query << "Enzyme='" << digestion->getEnzyme() << "'";
				treatment_query << ",DigestionTime=" << digestion->getDigestionTime();
				treatment_query << ",Ph=" << digestion->getPh();
				treatment_query << ",Temperature=" << digestion->getTemperature();
				result = db_con_.executeQuery(treatment_query.str());
				
				query << ",fid_Digestion=" << db_con_.getAutoId();
				result = db_con_.executeQuery(query.str());
				meta_id = db_con_.getAutoId();
				
				storeMetaInfo_("META_SampleTreatment", meta_id, *digestion);
			}
			if (sample.getTreatment(i).getType() == "Modification")
			{
				modification = dynamic_cast<const Modification*>(&sample.getTreatment(i));
				treatment_query << "INSERT INTO META_Modification SET ";
				treatment_query << "ReagentName='" << modification->getReagentName() << "'";
				treatment_query << ",AffectedAminoAcids='" << modification->getAffectedAminoAcids() << "'";
				treatment_query << ",SpecificityType=" << (1u+modification->getSpecificityType());
				treatment_query << ",Mass=" << modification->getMass();
				result = db_con_.executeQuery(treatment_query.str());
				
				query << ",fid_Modification=" << db_con_.getAutoId();
				result = db_con_.executeQuery(query.str());
				meta_id = db_con_.getAutoId();
				
				storeMetaInfo_("META_SampleTreatment", meta_id, *modification);
			}
			if (sample.getTreatment(i).getType() == "Tagging")
			{
							tagging = dynamic_cast<const Tagging*>(&sample.getTreatment(i));
							// tagging goes into META_Modification
							treatment_query << "INSERT INTO META_Modification SET ";
							treatment_query << "ReagentName='" << tagging->getReagentName() << "'";
							treatment_query << ",AffectedAminoAcids='" << tagging->getAffectedAminoAcids() << "'";
							treatment_query << ",SpecificityType=" << (1u+tagging->getSpecificityType());
							treatment_query << ",MassShift=" << tagging->getMassShift();
							treatment_query << ",Variant=" << (1u+tagging->getVariant());
							treatment_query << ",Mass=" << tagging->getMass();
							result = db_con_.executeQuery(treatment_query.str());

							query << ",fid_Modification=" << db_con_.getAutoId();
							result = db_con_.executeQuery(query.str());
							meta_id = db_con_.getAutoId();

							storeMetaInfo_("META_SampleTreatment", meta_id, *tagging);
			}
		}
		
		// recursively save subsamples
		for (std::vector<Sample>::const_iterator sample_it = sample.getSubsamples().begin(); sample_it != sample.getSubsamples().end(); sample_it++)
		{
			storeSample_ (*sample_it, exp_id, parent_id);
		}
		
		return parent_id;
	}
	
	void DBAdapter::loadSample_(const UID id, Sample& sample)
	{
		QSqlQuery result, sub_result;
		stringstream query;
		Sample subsample;
		std::vector<Sample> subsamples;
		UID meta_id;
		
		query.str("");
		query << "SELECT Name,SampleID,Mass,Volume,Concentration,State-1,Organism,Description,fid_MetaInfo FROM META_Sample WHERE id=" << id;
		result = db_con_.executeQuery(query.str(),true);
		
		sample.setName(result.value(0).toString());
		sample.setNumber(result.value(1).toString());
		sample.setMass(result.value(2).toDouble());
		sample.setVolume(result.value(3).toDouble());
		sample.setConcentration(result.value(4).toDouble());
		sample.setState((Sample::SampleState) result.value(5).toInt());
		sample.setOrganism(result.value(6).toString());
		sample.setComment(result.value(7).toString());
		loadMetaInfo_(result.value(8).toInt(),sample);
		
		// loading Treatments
		Digestion digestion;
		Modification modification;
		Tagging tagging;
		
		query.str("");
		query << "SELECT id,fid_Digestion,fid_Modification,Description,fid_MetaInfo FROM META_SampleTreatment WHERE fid_Sample=" << id;
		result = db_con_.executeQuery(query.str());

		while (result.next())
		{
			meta_id = result.value(4).toInt();
			// we got a digestion
			if (result.value(1).toInt() > 0)
			{
				loadMetaInfo_(meta_id, digestion);
				
				query.str("");
				query << "SELECT Enzyme,DigestionTime,Ph,Temperature FROM META_Digestion WHERE id=" << result.value(1).toInt();
				sub_result = db_con_.executeQuery(query.str(),true);
				
				digestion.setEnzyme(sub_result.value(0).toString());
				digestion.setDigestionTime(sub_result.value(1).toDouble());
				digestion.setPh(sub_result.value(2).toDouble());
				digestion.setTemperature(sub_result.value(3).toDouble());
				sample.addTreatment(digestion);
			}
			// we got a modification OR a tagging
			else if (result.value(2).toInt() > 0)
			{
				// build query and boolean function to distinguish between tagging and modification (NULL values)
				query.str("");
				query << "SELECT ReagentName,AffectedAminoAcids,SpecificityType-1,Mass,MassShift,Variant-1,MassShift IS NOT NULL AND Variant IS NOT NULL FROM META_Modification WHERE id=" << result.value(2).toInt();
				sub_result = db_con_.executeQuery(query.str(),true);

				// distinguish whether we are dealing with a tagging
				if (sub_result.value(6).toInt() == 1)
				{
					loadMetaInfo_(meta_id, tagging);
					tagging.setReagentName(sub_result.value(0).toString());
					tagging.setAffectedAminoAcids(sub_result.value(1).toString());
					tagging.setSpecificityType((Modification::SpecificityType) sub_result.value(2).toInt());
					tagging.setMass(sub_result.value(3).toDouble());
					tagging.setMassShift(sub_result.value(4).toDouble());
					tagging.setVariant((Tagging::IsotopeVariant) sub_result.value(5).toInt());
					sample.addTreatment(tagging);
				}
				else
				// we have a real modification
				{
					loadMetaInfo_(meta_id, modification);
					modification.setReagentName(sub_result.value(0).toString());
					modification.setAffectedAminoAcids(sub_result.value(1).toString());
					modification.setSpecificityType((Modification::SpecificityType) sub_result.value(2).toInt());
					modification.setMass(sub_result.value(3).toDouble());
					sample.addTreatment(modification);
				}
			}
		}

		query.str("");
		query << "SELECT id FROM META_Sample WHERE fid_Sample='" << id << "'";
		result = db_con_.executeQuery(query.str());
		
		while (result.next())
		{
			subsample = Sample();
			loadSample_(result.value(0).toInt(), subsample);
			subsamples.push_back(subsample);
		}
		sample.setSubsamples(subsamples);
	}
	
	bool DBAdapter::checkDBVersion(bool warning)
	{
		QSqlQuery result;
		try
		{
			result = db_con_.executeQuery("SELECT Version FROM ADMIN_Version", true);
		}
		catch(DBConnection::InvalidQuery)
		{
			cerr << "Error: This is no OpenMS DB, as there is no table 'ADMIN_Version'!" << endl;
			return false;
		}
		
		if (result.size()==0)
		{
			if (warning)
			{
				cerr << "Error: There is no entry in the 'ADMIN_Version' table. This should not happen!" << endl;
			}
			return false;
		}
		else if (result.size()>1)
		{
			if (warning)
			{
				cerr << "Error: There are several entries in the 'ADMIN_Version' table. This should not happen!" << endl;
			}
			return false;
		}
		else
		{
			String db_version = result.value(0).toString();
			db_version = db_version.suffix(':');
			db_version = db_version.prefix('$');
			db_version.trim();
			
			String sql_path;
			try
			{
				sql_path = File::find("OpenMS_DB.sql");
			}
			catch(Exception::FileNotFound&)
			{
				cerr << "Warning: Could not verify DB version. Please set the environment variable OPENMS_DATA_PATH to the OpenMS data directory: $PREFIX/share/OpenMS/" << endl;
				return true;
			}
				
			TextFile sql(sql_path);
			// ReverseIterator instead or ConstReverseIterator as rend() it non-const and operator!= fails (until gcc 4.0.x)
			for (TextFile::ReverseIterator it = sql.rbegin(); it != sql.rend(); ++it)
			{
				if (it->hasSubstring("$Revision:"))
				{
					String file_version = *it;
					file_version = file_version.suffix(':');
					file_version = file_version.prefix('$');
					file_version.trim();
					if (file_version!=db_version)
					{
						if (warning)
						{
							cerr << "Error: The given DB (Rev: " << db_version  << ") has a different revision than OpenMS (Rev: " << file_version << ")!"<< endl;
						}
						return false;
					}
					break;
				}
			}
		}
		return true;
	}
	
	void DBAdapter::createDB()
	{
		String sql_path;
		try
		{
			sql_path = File::find("OpenMS_DB.sql");
		}
		catch(Exception::FileNotFound&)
		{
			cerr << "Error: Could not find the OpenMS DB declaration file. Please set the environment variable OPENMS_DATA_PATH to the OpenMS data directory: $PREFIX/share/OpenMS/" << endl;
			return;
		}

		// load sql queries
		TextFile sql(sql_path);

		/*
		// delete existing tables
		QSqlQuery result, dummy;
		result = db_con_.executeQuery("SHOW TABLES;");
		dummy = db_con_.executeQuery("SET FOREIGN_KEY_CHECKS=0;");
		while (result.next())
		{
			dummy = db_con_.executeQuery(String("DROP TABLE `") + result.value(0).toString() + "`;");
		}
		dummy = db_con_.executeQuery("SET FOREIGN_KEY_CHECKS=1;");
		*/

		// Conversion of phpMyAdmin output to required format
		// concatenate lines so that one line is one query
		vector<String> queries;
		String query, line;
		
		bool in_query = false;
		for (TextFile::ConstIterator it = sql.begin(); it != sql.end(); ++it)
		{
			line = *it;
			line.trim();
			
			if (in_query)
			{
				// line is empty or comment => query ends
				if (line=="" || (line.size()>=2 && line[0]=='-' && line[1]=='-'))
				{
					queries.push_back(query);
					in_query = false;
				}
				// query continues => append
				else
				{
					query = query + ' ' + line;
				}
			}
			else
			{
				// line is not empty and not a comment => query starts
				if (line!="" && (line.size()<2 || line[0]!='-' || line[1]!='-'))
				{
					query = line;
					in_query = true;
				}
			}
		}
		
		// execute queries
		db_con_.executeQueries(queries);
	}

	PeakFileOptions& DBAdapter::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& DBAdapter::getOptions() const
  {
  	return options_;
  }


	void DBAdapter::loadMetaInfo_(UID, Peak1D&)
	{
	}

	void DBAdapter::loadMetaInfo_(UID id, RichPeak1D& peak)
	{
		loadMetaInfo_(id, static_cast<MetaInfoInterface&>(peak));
	}

	void DBAdapter::loadMetaInfo_(UID id, Precursor& peak)
	{
		loadMetaInfo_(id, static_cast<MetaInfoInterface&>(peak));
	}
} //namespace
