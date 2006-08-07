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
// $Id: InspectInfile.C,v 1.0 2006/07/25 13:46:15 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectInfile.h>

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
	InspectInfile::InspectInfile(const InspectInfile& inspect_infile):
		spectra_(inspect_infile.getSpectra()),
		db_(inspect_infile.getDb()),
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
			db_ = inspect_infile.getDb();
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
	
	void InspectInfile::store(const String& filename) throw (Exception::UnableToCreateFile)
	{
		std::ofstream ofs( filename.c_str() );
		if ( !ofs ) throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		
		// the spectrum
		ofs << "spectra," << spectra_ << std::endl;
			
		// the database
		if ( !db_.empty() ) ofs << "db," << db_ << std::endl;
		
		// the sequence file
		if ( !sequence_file_.empty() ) ofs << "sequence_file," << sequence_file_ << std::endl;
		
		// the protease ###
		if ( !protease_.empty() ) ofs << "protease," << protease_ << std::endl;
		
		// the modifications
		for ( std::vector< std::vector< String > >::const_iterator iter = mod_.begin(); iter != mod_.end(); ++iter )
		{
			ofs << "mod";
			for ( std::vector< String >::const_iterator i = iter->begin(); i != iter->end(); ++i)
			{
				ofs << "," << *i;
			}
			ofs << std::endl;
		}
		
		// number of modifications per peptide
		if ( mods_ > -1 ) ofs << "mods," << mods_ << std::endl;
		
		// whether to do a "blind" search
		if ( blind_ != 2 ) ofs << "blind," << blind_ << std::endl;

		// maximum modification size in a "blind" search
		if ( maxptmsize_ >= 0) ofs << "maxptmsize," << maxptmsize_ << std::endl;
		
		// parent mass tolerance
		if ( PM_tolerance_ >= 0 ) ofs << "PM_tolerance," << PM_tolerance_ << std::endl;
		
		// ion mass tolerance
		if ( ion_tolerance_ >= 0 ) ofs << "ion_tolerance," << ion_tolerance_ << std::endl;
		
		// jumpscores file
		if ( !jumpscores_.empty() ) ofs << "jumpscores," << jumpscores_ << std::endl;
		
		// whether precursor ion may be multiply charged
		if ( multicharge_ != 2 ) ofs << "multicharge," << multicharge_ << std::endl;
		
		// instrument type
		if ( !instrument_.empty() ) ofs << "instrument," << instrument_ << std::endl;
		
		// number of tags to generate for the first pass of a two-pass search
		if ( tag_count_a_ >= 0 ) ofs << "TagCountA," << tag_count_a_ << std::endl;
		
		// Number of tags to generate for the second pass of a two-pass search, OR the number of tags to use in a one-pass search.
		if ( tag_count_a_ >= 0 ) ofs << "TagCountB," << tag_count_b_ << std::endl;
		
		// whether to use two-pass search
		if ( twopass_ != 2 ) ofs << "twopass," << twopass_ << std::endl;
		
		ofs.close();
		ofs.clear();
	}


	const std::string& InspectInfile::getSpectra() const
	{
		return spectra_;
	}

	void InspectInfile::setSpectra(const std::string& spectra)
	{
		spectra_ = spectra;
	}

	const String& InspectInfile::getDb() const
	{
		return db_;
	}

	void InspectInfile::setDb(const String& db)
	{
		db_ = db;
	}

	const String& InspectInfile::getSequenceFile() const
	{
		return sequence_file_;
	}

	void InspectInfile::setSequenceFile(const String& sequence_file)
	{
		sequence_file_ = sequence_file;
	}

	const String& InspectInfile::getProtease() const
	{
		return protease_;
	}

	void InspectInfile::setProtease(const String& protease)
	{
		protease_ = protease;
	}
	
	const std::vector< std::vector< String > >& InspectInfile::getMod() const
	{
		return mod_;
	}
	
	std::vector< std::vector< String > >& InspectInfile::getMod()
	{
		return mod_;
	}

	void InspectInfile::setMod(const std::vector< std::vector< String > >& mod)
	{
		mod_ = mod;
	}
	
	void InspectInfile::addMod(const std::vector< String >& mod)
	{
		mod_.push_back(mod);
	}

	const int InspectInfile::getMods() const
	{
		return mods_;
	}

	void InspectInfile::setMods(int mods)
	{
		mods_ = mods;
	}

	const unsigned int InspectInfile::getBlind() const
	{
		return blind_;
	}

	void InspectInfile::setBlind(unsigned int blind)
	{
		blind_ = blind;
	}

	const double InspectInfile::getMaxPTMsize() const
	{
		return maxptmsize_;
	}

	void InspectInfile::setMaxPTMsize(double maxptmsize)
	{
		maxptmsize_ = maxptmsize;
	}

	const double InspectInfile::getPMTolerance() const
	{
		return PM_tolerance_;
	}

	void InspectInfile::setPMTolerance(double PM_tolerance)
	{
		PM_tolerance_ = PM_tolerance;
	}

	const double InspectInfile::getIonTolerance() const
	{
		return ion_tolerance_;
	}

	void InspectInfile::setIonTolerance(double ion_tolerance)
	{
		ion_tolerance_ = ion_tolerance;
	}

	const String& InspectInfile::getJumpscores() const
	{
		return jumpscores_;
	}

	void InspectInfile::setJumpscores(const String& jumpscores)
	{
		jumpscores_ = jumpscores;
	}

	const unsigned int InspectInfile::getMulticharge() const
	{
		return multicharge_;
	}

	void InspectInfile::setMulticharge(unsigned int multicharge)
	{
std::cout << multicharge << std::endl;
		multicharge_ = multicharge;
	}

	const String& InspectInfile::getInstrument() const
	{
		return instrument_;
	}

	void InspectInfile::setInstrument(const String& instrument)
	{
		instrument_ = instrument;
	}

	const int InspectInfile::getTagCountA() const
	{
		return tag_count_a_;
	}

	void InspectInfile::setTagCountA(int tag_count_a)
	{
		tag_count_a_ = tag_count_a;
	}

	const int InspectInfile::getTagCountB() const
	{
		return tag_count_b_;
	}

	void InspectInfile::setTagCountB(int tag_count_b)
	{
		tag_count_b_ = tag_count_b;
	}

	const unsigned int InspectInfile::getTwopass() const
	{
		return twopass_;
	}

	void InspectInfile::setTwopass(unsigned int twopass)
	{
		twopass_ = twopass;
	}
	
	void InspectInfile::generateSecondDatabase(const std::string& result_filename, const std::string& result_path, const std::string&  database_path, const std::string& database_filename, const double& cutoff_p_value, const double& cutoff_score_value, int min_annotated_spectra_per_protein, std::string second_database_filename, std::string second_index_filename, std::string second_database_path, std::string index_filename, std::string species) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument)
	{
		if ( (cutoff_p_value < 0) || (cutoff_p_value > 1) )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "cutoff_p_value is lower 0 or greater 1!");
		}
		std::string path_and_file = result_path;
		ensurePathChar(path_and_file);
		path_and_file.append(result_filename);
		std::ifstream result_file( path_and_file.c_str(), std::ios::in | std::ios::binary );
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		path_and_file = database_path;
		ensurePathChar(path_and_file);
		path_and_file.append(database_filename);
		std::ifstream database_file( path_and_file.c_str(), std::ios::in | std::ios::binary );
		if ( !database_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		
		std::string start_separator, buffer1;
		// get the start separator
		getLabels(path_and_file, buffer1, start_separator, buffer1, buffer1, buffer1);
		
		// processing the information from the result file
		enum columns
		{
			spectrum_file_column,
			scan_column,
			peptide_column,
			protein_column,
			charge_column,
			MQ_score_column,
			cut_score_column,
			intense_by_column,
			by_present_column,
			unused_column,
			p_value_column,
			delta_score_column,
			delta_score_other_column,
			record_number_column,
			DB_file_pos_column,
			spec_file_pos_column
		};
		unsigned int number_of_columns = 16;
		
		std::vector<String> substrings;
		String line;
		
		// map the proteins according to their record number in the result file (and count the number of annotated spectrum)
		std::map< unsigned int, unsigned int > record_map;
		std::map< unsigned int, unsigned int >::iterator i;
		
		// get the proteins whose p-value <= p_value_column and count the annotated spectrum (done via a set)
		unsigned int record_number;
		unsigned int max_record_number = 0;
		unsigned int line_number = 0; // used to report in which line an error occured
		bool missing_column;
		//std::set<std::string, string_less> spectrum_count;
		std::set< std::string > spectrum_count;
		//char buffer[10];
		
		// read out the whole result file
		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			++line_number;
			line.split('\t', substrings);
			// check whether the line has enough columns
			missing_column = ( substrings.size() == number_of_columns - 1 );
			if ( !missing_column && (substrings.size() < number_of_columns) ) continue;
			
			// if the version Inspect.20060620.zip is used, there is a header
			if ( substrings[0] == "#SpectrumFile" ) continue;
			
			spectrum_count.insert(substrings[spectrum_file_column]);
			
			// if the p_value of this record is lower or equal to the cutoff it is inserted or it's number of annotated spectrum is increased
			p_value = substrings[p_value_column - missing_column];
			
			if ( (p_value.length() >= 5) && (p_value.substr(0,5) == "0.000") && (substrings[p_value_column] == "nan") ) p_value = "nan";
			if ( ((p_value == "nan") && (atoi(substrings[MQ_score_column - missing_column].c_str()) < score_value_threshold)) || (atof(p_value.c_str()) > p_value_threshold) ) continue;
			
			record_number = atoi(substrings[record_number_column - missing_column].c_str()) - missing_column;
			max_record_number = std::max(max_record_number, record_number);
			// if the record has already been inserted it's number of annotated spectrum is increased  otherwise it is inserted
			++record_map[record_number];
		} // result file read
		result_file.close();
		result_file.clear();
		
		// if no protein has a p_value less equal the cutoff value return an empty database
		if ( record_map.empty() )
		{
			path_and_file = second_database_path;
			ensurePathChar(path_and_file);
			path_and_file.append(second_database_filename);
			FILE* file = fopen(path_and_file.c_str(), "w");
			if ( database_file == NULL )
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_database_filename);
			}
			fclose(file);
			
			path_and_file = second_database_path;
			ensurePathChar(path_and_file);
			path_and_file.append(second_index_filename);
			file = fopen(path_and_file.c_str(), "w");
			if ( database_file == NULL )
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_index_filename);
			}
			fclose(file);
			return;
		}
		
		// get the number of proteins
		unsigned int number_of_proteins = 0;
		unsigned int pos = 0;
		
		while ( getline(database_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			pos = line.find(start_separator, pos);
			while ( pos != std::string::npos )
			{
				++number_of_proteins;
				pos = line.find(start_separator, ++pos);
			}
		}
		// a trie datase has one more protein than delimiters
		if ( start_separator == std::string(1, trie_delimiter_) ) ++ number_of_proteins;
		
		// check whether the number of proteins in the database and the result file correspond
		if ( number_of_proteins < max_record_number+1 )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong database file ("+database_filename+") for this result file ("+result_filename+") (number of proteins in database is lower than a RecordNumber in the result file)!" ,result_filename);
		}
		
		// generate the new database using all proteins whose #annotated spectrum >= min_annotated_spectrum_per_peptid
		// compute min_annotated_spectra_per_protein if it's not given
		if ( min_annotated_spectra_per_protein < 0 )
		{
			min_annotated_spectra_per_protein = std::max(1, (signed) (2 * spectrum_count.size() / number_of_proteins));
		}
		spectrum_count.clear();
		
		// insert those records into a vector whose #annotated spectrum >= min_annotated_spectrum_per_peptid
		std::vector< unsigned int > wanted_records;
		for ( std::map< unsigned int, unsigned int >::const_iterator i = record_map.begin(); i != record_map.end(); ++i )
		{
			if ( i->second >= static_cast< unsigned int >(min_annotated_spectra_per_protein) ) wanted_records.push_back(i->first);
		}
		record_map.clear();
		database_file.close();
		result_file.close();
		
		// if no protein has a p_value less equal the cutoff value return an empty database
		if ( wanted_records.empty() )
		{
			path_and_file = second_database_path;
			ensurePathChar(path_and_file);
			path_and_file.append(second_database_filename);
			FILE* file = fopen(second_database_filename.c_str(), "w");
			if ( database_file == NULL )
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_database_filename);
			}
			fclose(file);
			
			path_and_file = second_database_path;
			ensurePathChar(path_and_file);
			path_and_file.append(second_index_filename);
			file = fopen(second_index_filename.c_str(), "w");
			if ( database_file == NULL )
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_index_filename);
			}
			fclose(file);
			return;
		}
		
		if (start_separator == std::string(1, trie_delimiter_)) compressTrieDB(database_filename, index_filename, database_path, wanted_records, second_database_filename, second_index_filename, second_database_path);
		else generateTrieDB(database_filename, database_path, second_database_path, wanted_records, second_database_filename, second_index_filename, false, species);
	}


} // namespace OpenMS
