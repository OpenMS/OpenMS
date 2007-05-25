// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectInfile.h>

using namespace std;

namespace OpenMS
{
	
	// default constructor
	InspectInfile::InspectInfile():
		mods_(-1),
		blind_(2),
		maxptmsize_(-1.0),
		precursor_mass_tolerance_(-1.0),
		peak_mass_tolerance_(-1.0),
		multicharge_(2),
		tag_count_(-1)
	{
	}
	
	// copy constructor
	InspectInfile::InspectInfile(
		const InspectInfile& inspect_infile):
		spectra_(inspect_infile.getSpectra()),
		enzyme_(inspect_infile.getEnzyme()),
		mod_(inspect_infile.getMod()),
		mods_(inspect_infile.getMods()),
		blind_(inspect_infile.getBlind()),
		maxptmsize_(inspect_infile.getMaxPTMsize()),
		precursor_mass_tolerance_(inspect_infile.getPrecursorMassTolerance()),
		peak_mass_tolerance_(inspect_infile.getPeakMassTolerance()),
		multicharge_(inspect_infile.getMulticharge()),
		instrument_(inspect_infile.getInstrument()),
		tag_count_(inspect_infile.getTagCount())
	{
	}
	
	// destructor
	InspectInfile::~InspectInfile()
	{
		mod_.clear();
	}
	
	// assignment operator
	InspectInfile& InspectInfile::operator= (const InspectInfile& inspect_infile)
	{
		if (this != &inspect_infile)
		{
			spectra_ = inspect_infile.getSpectra();
			enzyme_ = inspect_infile.getEnzyme();
			mod_ = inspect_infile.getMod();
			mods_ = inspect_infile.getMods();
			blind_ = inspect_infile.getBlind();
			maxptmsize_ = inspect_infile.getMaxPTMsize();
			precursor_mass_tolerance_ = inspect_infile.getPrecursorMassTolerance();
			peak_mass_tolerance_ = inspect_infile.getPeakMassTolerance();
			multicharge_ = inspect_infile.getMulticharge();
			instrument_ = inspect_infile.getInstrument();
			tag_count_ = inspect_infile.getTagCount();
		}
		return *this;
	}

	void
	InspectInfile::store(
		const String& filename)
	throw (
		Exception::UnableToCreateFile)
	{
		ofstream ofs( filename.c_str() );
		if ( !ofs ) throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);

		ofs << "spectra," << spectra_ << endl;

		if ( !db_.empty() ) ofs << "db," << db_ << endl;

		if ( !enzyme_.empty() ) ofs << "protease," << enzyme_ << endl;

		for ( vector< vector< String > >::const_iterator iter = mod_.begin(); iter != mod_.end(); ++iter )
		{
			ofs << "mod";
			for ( vector< String >::const_iterator i = iter->begin(); i != iter->end(); ++i)
			{
				if ( !i->empty()) ofs << "," << *i;
			}
			ofs << endl;
		}

		if ( mods_ > -1 ) ofs << "mods," << mods_ << endl;

		if ( blind_ != 2 ) ofs << "blind," << blind_ << endl;

		if ( maxptmsize_ >= 0) ofs << "maxptmsize," << maxptmsize_ << endl;

		if ( precursor_mass_tolerance_ >= 0 ) ofs << "PM_tolerance," << precursor_mass_tolerance_ << endl;

		if ( peak_mass_tolerance_ >= 0 ) ofs << "IonTolerance," << peak_mass_tolerance_ << endl;

		if ( multicharge_ != 2 ) ofs << "multicharge," << multicharge_ << endl;

		if ( !instrument_.empty() ) ofs << "instrument," << instrument_ << endl;

		if ( tag_count_ >= 0 ) ofs << "TagCount," << tag_count_ << endl;

		ofs.close();
		ofs.clear();
	}


	const string& InspectInfile::getSpectra() const {return spectra_;}
	void InspectInfile::setSpectra(const string& spectra) {spectra_ = spectra;}

	const String& InspectInfile::getDb() const {return db_;}
	void InspectInfile::setDb(const String& db) {db_ = db;}

	const String& InspectInfile::getEnzyme() const {return enzyme_;}
	void InspectInfile::setEnzyme(const String& enzyme) {enzyme_ = enzyme;}
	
	const vector< vector< String > >& InspectInfile::getMod() const {return mod_;}
	vector< vector< String > >& InspectInfile::getMod() {return mod_;}

	void InspectInfile::setMod(const vector< vector< String > >& mod) {mod_ = mod;}
	void InspectInfile::addMod(const vector< String >& mod) {mod_.push_back(mod);}

	const Int InspectInfile::getMods() const {return mods_;}
	void InspectInfile::setMods(Int mods) {mods_ = mods;}

	const UInt InspectInfile::getBlind() const {return blind_;}
	void InspectInfile::setBlind(UInt blind) {blind_ = blind;}

	const DoubleReal InspectInfile::getMaxPTMsize() const {return maxptmsize_;}
	void InspectInfile::setMaxPTMsize(DoubleReal maxptmsize) {maxptmsize_ = maxptmsize;}

	const DoubleReal InspectInfile::getPrecursorMassTolerance() const {return precursor_mass_tolerance_;}
	void InspectInfile::setPrecursorMassTolerance(DoubleReal precursor_mass_tolerance) {precursor_mass_tolerance_ = precursor_mass_tolerance;}

	const DoubleReal InspectInfile::getPeakMassTolerance() const {return peak_mass_tolerance_;}
	void InspectInfile::setPeakMassTolerance(DoubleReal ion_tolerance) {peak_mass_tolerance_ = ion_tolerance;}

	const UInt InspectInfile::getMulticharge() const {return multicharge_;}
	void InspectInfile::setMulticharge(UInt multicharge) {multicharge_ = multicharge;}

	const String& InspectInfile::getInstrument() const {return instrument_;}
	void InspectInfile::setInstrument(const String& instrument) {instrument_ = instrument;}

	const Int InspectInfile::getTagCount() const {return tag_count_;}
	void InspectInfile::setTagCount(Int tag_count) {tag_count_ = tag_count;}
} // namespace OpenMS
