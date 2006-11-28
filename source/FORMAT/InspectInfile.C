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
		PM_tolerance_(-1.0),
		ion_tolerance_(-1.0),
		multicharge_(2),
		tag_count_a_(-1),
		tag_count_b_(-1),
		twopass_(2)
	{
	}
	
	// copy constructor
	InspectInfile::InspectInfile(
		const InspectInfile& inspect_infile):
		spectra_(inspect_infile.getSpectra()),
		sequence_file_(inspect_infile.getSequenceFile()),
		protease_(inspect_infile.getProtease()),
		mod_(inspect_infile.getMod()),
		mods_(inspect_infile.getMods()),
		blind_(inspect_infile.getBlind()),
		maxptmsize_(inspect_infile.getMaxPTMsize()),
		PM_tolerance_(inspect_infile.getPMTolerance()),
		ion_tolerance_(inspect_infile.getIonTolerance()),
		jumpscores_(inspect_infile.getJumpscores()),
		multicharge_(inspect_infile.getMulticharge()),
		instrument_(inspect_infile.getInstrument()),
		tag_count_a_(inspect_infile.getTagCountA()),
		tag_count_b_(inspect_infile.getTagCountB()),
		twopass_(inspect_infile.getTwopass())
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
			sequence_file_ = inspect_infile.getSequenceFile();
			protease_ = inspect_infile.getProtease();
			mod_ = inspect_infile.getMod();
			mods_ = inspect_infile.getMods();
			blind_ = inspect_infile.getBlind();
			maxptmsize_ = inspect_infile.getMaxPTMsize();
			PM_tolerance_ = inspect_infile.getPMTolerance();
			ion_tolerance_ = inspect_infile.getIonTolerance();
			jumpscores_ = inspect_infile.getJumpscores();
			multicharge_ = inspect_infile.getMulticharge();
			instrument_ = inspect_infile.getInstrument();
			tag_count_a_ = inspect_infile.getTagCountA();
			tag_count_b_ = inspect_infile.getTagCountB();
			twopass_ = inspect_infile.getTwopass();
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

		if ( !sequence_file_.empty() ) ofs << "sequence_file," << sequence_file_ << endl;

		if ( !protease_.empty() ) ofs << "protease," << protease_ << endl;

		for ( vector< vector< String > >::const_iterator iter = mod_.begin(); iter != mod_.end(); ++iter )
		{
			ofs << "mod";
			for ( vector< String >::const_iterator i = iter->begin(); i != iter->end(); ++i)
			{
				ofs << "," << *i;
			}
			ofs << endl;
		}

		if ( mods_ > -1 ) ofs << "mods," << mods_ << endl;

		if ( blind_ != 2 ) ofs << "blind," << blind_ << endl;

		if ( maxptmsize_ >= 0) ofs << "maxptmsize," << maxptmsize_ << endl;

		if ( PM_tolerance_ >= 0 ) ofs << "PM_tolerance," << PM_tolerance_ << endl;

		if ( ion_tolerance_ >= 0 ) ofs << "IonTolerance," << ion_tolerance_ << endl;

		if ( !jumpscores_.empty() ) ofs << "jumpscores," << jumpscores_ << endl;

		if ( multicharge_ != 2 ) ofs << "multicharge," << multicharge_ << endl;

		if ( !instrument_.empty() ) ofs << "instrument," << instrument_ << endl;

		if ( tag_count_a_ >= 0 ) ofs << "TagCountA," << tag_count_a_ << endl;

		if ( tag_count_a_ >= 0 ) ofs << "TagCountB," << tag_count_b_ << endl;

		if ( twopass_ != 2 ) ofs << "twopass," << twopass_ << endl;

		ofs.close();
		ofs.clear();
	}


	const string& InspectInfile::getSpectra() const {return spectra_;}

	void InspectInfile::setSpectra(const string& spectra) {spectra_ = spectra;}

	const String& InspectInfile::getDb() const {return db_;}

	void InspectInfile::setDb(const String& db) {db_ = db;}

	const String& InspectInfile::getSequenceFile() const {return sequence_file_;}

	void InspectInfile::setSequenceFile(const String& sequence_file) {sequence_file_ = sequence_file;}

	const String& InspectInfile::getProtease() const {return protease_;}

	void InspectInfile::setProtease(const String& protease) {protease_ = protease;}
	
	const vector< vector< String > >& InspectInfile::getMod() const {return mod_;}
	
	vector< vector< String > >& InspectInfile::getMod() {return mod_;}

	void InspectInfile::setMod(const vector< vector< String > >& mod) {mod_ = mod;}
	
	void InspectInfile::addMod(const vector< String >& mod) {mod_.push_back(mod);}

	const SignedInt InspectInfile::getMods() const {return mods_;}

	void InspectInfile::setMods(SignedInt mods) {mods_ = mods;}

	const UnsignedInt InspectInfile::getBlind() const {return blind_;}

	void InspectInfile::setBlind(UnsignedInt blind) {blind_ = blind;}

	const DoubleReal InspectInfile::getMaxPTMsize() const {return maxptmsize_;}

	void InspectInfile::setMaxPTMsize(DoubleReal maxptmsize) {maxptmsize_ = maxptmsize;}

	const DoubleReal InspectInfile::getPMTolerance() const {return PM_tolerance_;}

	void InspectInfile::setPMTolerance(DoubleReal PM_tolerance) {PM_tolerance_ = PM_tolerance;}

	const DoubleReal InspectInfile::getIonTolerance() const {return ion_tolerance_;}

	void InspectInfile::setIonTolerance(DoubleReal ion_tolerance) {ion_tolerance_ = ion_tolerance;}

	const String& InspectInfile::getJumpscores() const {return jumpscores_;}

	void InspectInfile::setJumpscores(const String& jumpscores) {jumpscores_ = jumpscores;}

	const UnsignedInt InspectInfile::getMulticharge() const {return multicharge_;}

	void InspectInfile::setMulticharge(UnsignedInt multicharge) {multicharge_ = multicharge;}

	const String& InspectInfile::getInstrument() const {return instrument_;}

	void InspectInfile::setInstrument(const String& instrument) {instrument_ = instrument;}

	const SignedInt InspectInfile::getTagCountA() const {return tag_count_a_;}

	void InspectInfile::setTagCountA(SignedInt tag_count_a) {tag_count_a_ = tag_count_a;}

	const SignedInt InspectInfile::getTagCountB() const {return tag_count_b_;}

	void InspectInfile::setTagCountB(SignedInt tag_count_b) {tag_count_b_ = tag_count_b;}

	const UnsignedInt InspectInfile::getTwopass() const {return twopass_;}

	void InspectInfile::setTwopass(UnsignedInt twopass) {twopass_ = twopass;}
} // namespace OpenMS
