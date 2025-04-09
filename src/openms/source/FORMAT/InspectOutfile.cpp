// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#if defined OPENMS_BIG_ENDIAN
#define OPENMS_IS_BIG_ENDIAN true
#else
#define OPENMS_IS_BIG_ENDIAN false
#endif

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <QtCore/QRegularExpression>

#include <fstream>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"

using namespace std;

namespace OpenMS
{
  InspectOutfile::InspectOutfile() = default;

  /// copy constructor
  InspectOutfile::InspectOutfile(const InspectOutfile&) = default;

  /// destructor
  InspectOutfile::~InspectOutfile() = default;

  /// assignment operator
  InspectOutfile& InspectOutfile::operator=(const InspectOutfile& inspect_outfile)
  {
    if (this == &inspect_outfile)
    {
      return *this;
    }
    return *this;
  }

  /// equality operator
  bool InspectOutfile::operator==(const InspectOutfile&) const
  {
    return true;
  }

  vector<Size> InspectOutfile::load(const String& result_filename, vector<PeptideIdentification>& peptide_identifications,
                                    ProteinIdentification& protein_identification, const double p_value_threshold, const String& database_filename)
  {
    // check whether the p_value is correct
    if ((p_value_threshold < 0) || (p_value_threshold > 1))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The parameters 'p_value_threshold' must be >= 0 and <=1 !");
    }

    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }

    String
      line,
      accession,
      accession_type,
      spectrum_file,
      identifier;

    Size
      record_number(0),
    scan_number(0),
    line_number(0),
    number_of_columns(0);

    vector<String> substrings;
    vector<Size> corrupted_lines;

    PeptideIdentification peptide_identification;

    if (!getline(result_file, line)) // the header is read in a special function, so it can be skipped
    {
      result_file.close();
      result_file.clear();
      throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }
    if (!line.empty() && (line[line.length() - 1] < 33))
      line.resize(line.length() - 1);
    line.trim();
    ++line_number;

    DateTime datetime = DateTime::now();
    if (protein_identification.getSearchEngine().empty())
    {
      identifier = "InsPecT_" + datetime.getDate();
    }
    else
    {
      protein_identification.getSearchEngine() + "_" + datetime.getDate();
    }
    // to get the precursor retention time and mz values later, save the filename and the numbers of the scans
    vector<pair<String, vector<pair<Size, Size> > > > files_and_peptide_identification_with_scan_number;
    // the record number is mapped to the position in the protein hits, to retrieve their sequences
    map<Size, Size> rn_position_map;

    // get the header
    Int
      spectrum_file_column(-1),
    scan_column(-1),
    peptide_column(-1),
    protein_column(-1),
    charge_column(-1),
    MQ_score_column(-1),
    p_value_column(-1),
    record_number_column(-1),
    DB_file_pos_column(-1),
    spec_file_pos_column(-1);

    String::size_type start(0), end(0);

    try
    {
      readOutHeader(result_filename, line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns);
    }
    catch (Exception::ParseError& p_e)
    {
      result_file.close();
      result_file.clear();
      OPENMS_LOG_WARN << "ParseError (" << p_e.what() << ") caught in " << __FILE__ << "\n";
      throw;
    }

    while (getline(result_file, line))
    {
      ++line_number;
      if (!line.empty() && (line[line.length() - 1] < 33))
        line.resize(line.length() - 1);
      line.trim();
      if (line.empty())
        continue;

      // check whether the line has enough columns
      line.split('\t', substrings);
      if (substrings.size() != number_of_columns)
      {
        corrupted_lines.push_back(line_number);
        continue;
      }

      // if the pvalue is too small, skip the line
      if (substrings[p_value_column].toFloat() > p_value_threshold)
      {
        continue;
      }
      // the protein
      ProteinHit protein_hit;
      // get accession number and type
      getACAndACType(substrings[protein_column], accession, accession_type);
      protein_hit.setAccession(accession);
//          protein_hit.setScore(0.0);

      // the database position of the protein (the i-th protein)
      record_number = substrings[record_number_column].toInt();

      // map the database position of the protein to its position in the
      // protein hits and insert it, if it's a new protein
      if (rn_position_map.find(record_number) == rn_position_map.end())
      {
        rn_position_map[record_number] = protein_identification.getHits().size();
        protein_identification.insertHit(protein_hit);
      }

      // if a new scan is found (new file or new scan), insert it into the
      // vector (the first time the condition is fulfilled because
      // spectrum_file is "")
      if ((substrings[spectrum_file_column] != spectrum_file) || ((Size) substrings[scan_column].toInt() != scan_number))
      {
        // if it's a new file, insert it into the vector (used to retrieve RT and MT later)
        if (substrings[spectrum_file_column] != spectrum_file)
        {
          // if it's the first file or if hits have been found in the file before, insert a new file
          if (files_and_peptide_identification_with_scan_number.empty() ||
              !files_and_peptide_identification_with_scan_number.back().second.empty())
          {
            files_and_peptide_identification_with_scan_number.emplace_back(substrings[spectrum_file_column],
                  vector<pair<Size, Size> >());
          }
          // otherwise change the name of the last file entry (the one without hits)
          else
            files_and_peptide_identification_with_scan_number.back().first = substrings[spectrum_file_column];
        }

        spectrum_file = substrings[spectrum_file_column];
        scan_number = substrings[scan_column].toInt();

        // if it's not the first scan and if hits have been found, insert the peptide identification
        if (!peptide_identification.empty() && !peptide_identification.getHits().empty())
        {
          files_and_peptide_identification_with_scan_number.back().second.emplace_back(peptide_identifications.size(), scan_number);
          peptide_identifications.push_back(peptide_identification);
        }
        peptide_identification = PeptideIdentification();

        peptide_identification.setIdentifier(identifier);
        peptide_identification.setSignificanceThreshold(p_value_threshold);
        peptide_identification.setScoreType(score_type_);
      }

      // get the peptide infos from the new peptide and insert it
      PeptideHit peptide_hit;
      peptide_hit.setCharge(substrings[charge_column].toInt());
      peptide_hit.setScore(substrings[MQ_score_column].toFloat());
      peptide_hit.setRank(0); // all ranks are set to zero and assigned later

      // get the sequence and the amino acid before and after
      String sequence, sequence_with_mods;
      sequence_with_mods = substrings[peptide_column];
      start = sequence_with_mods.find('.') + 1;
      end = sequence_with_mods.find_last_of('.');

      PeptideEvidence pe;

      if (start >= 2)
      {
        pe.setAABefore(sequence_with_mods[start - 2]);
      }

      if (end < sequence_with_mods.length() + 1)
      {
        pe.setAAAfter(sequence_with_mods[end + 1]);
      }

      //remove modifications (small characters and anything that's not in the alphabet)
      sequence_with_mods = substrings[peptide_column].substr(start, end - start);
      for (String::ConstIterator c_i = sequence_with_mods.begin(); c_i != sequence_with_mods.end(); ++c_i)
      {
        if ((bool) isalpha(*c_i) && (bool) isupper(*c_i))
          sequence.append(1, *c_i);
      }

      peptide_hit.setSequence(AASequence::fromString(sequence));
      pe.setProteinAccession(accession);
      peptide_hit.addPeptideEvidence(pe);

      peptide_identification.insertHit(peptide_hit);
    }

    // result file read
    result_file.close();
    result_file.clear();

    // if it's not the first scan and if hits have been found, insert the peptide identification
    if (!peptide_identification.empty() && !peptide_identification.getHits().empty())
    {
      files_and_peptide_identification_with_scan_number.back().second.emplace_back(peptide_identifications.size(), scan_number);
      peptide_identifications.push_back(peptide_identification);
    }

    // if the last file had no hits, delete it
    if (!files_and_peptide_identification_with_scan_number.empty() && files_and_peptide_identification_with_scan_number.back().second.empty())
    {
      files_and_peptide_identification_with_scan_number.pop_back();
    }

    if (!peptide_identifications.empty())
      peptide_identifications.back().assignRanks();

    // search the sequence of the proteins
    if (!protein_identification.getHits().empty() && !database_filename.empty())
    {
      vector<ProteinHit> protein_hits = protein_identification.getHits();
      vector<String> sequences;
      getSequences(database_filename, rn_position_map, sequences);

      // set the retrieved sequences
      vector<String>::const_iterator s_i = sequences.begin();
      for (map<Size, Size>::const_iterator rn_i = rn_position_map.begin(); rn_i != rn_position_map.end(); ++rn_i, ++s_i)
        protein_hits[rn_i->second].setSequence(*s_i);

      sequences.clear();
      rn_position_map.clear();
      protein_identification.setHits(protein_hits);
      protein_hits.clear();
    }

    // get the precursor retention times and mz values
    getPrecursorRTandMZ(files_and_peptide_identification_with_scan_number, peptide_identifications);
    protein_identification.setDateTime(datetime);
    protein_identification.setIdentifier(identifier);

    return corrupted_lines;
  }

  // < record number, number of protein in a vector >
  vector<Size>
  InspectOutfile::getSequences(
    const String& database_filename,
    const map<Size, Size>& wanted_records,
    vector<String>& sequences)
  {
    ifstream database(database_filename.c_str());
    if (!database)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, database_filename);
    }

    vector<Size> not_found;
    Size seen_records(0);
    stringbuf sequence;
    database.seekg(0, ios::end);
    streampos sp = database.tellg();
    database.seekg(0, ios::beg);

    for (map<Size, Size>::const_iterator wr_i = wanted_records.begin(); wr_i !=  wanted_records.end(); ++wr_i)
    {
      for (; seen_records < wr_i->first; ++seen_records)
      {
        database.ignore(sp, trie_delimiter_);
      }
      database.get(sequence, trie_delimiter_);
      sequences.emplace_back(sequence.str());
      if (sequences.back().empty())
      {
        not_found.push_back(wr_i->first);
      }
      sequence.str("");
    }

    // close the filestreams
    database.close();
    database.clear();

    return not_found;
  }

  void
  InspectOutfile::getACAndACType(
    String line,
    String& accession,
    String& accession_type)
  {
    String swissprot_prefixes = "JLOPQUX";
    /// @todo replace this by general FastA implementation? (Martin)
    accession.clear();
    accession_type.clear();

    // if it's a FASTA line
    if (line.hasPrefix(">"))
    {
      line.erase(0, 1);
    }
    if (!line.empty() && (line[line.length() - 1] < 33))
    {
      line.resize(line.length() - 1);
    }
    line.trim();

    // if it's a swissprot accession
    if (line.hasPrefix("tr") || line.hasPrefix("sp"))
    {
      accession = line.substr(3, line.find('|', 3) - 3);
      accession_type = "SwissProt";
    }
    else if (line.hasPrefix("gi"))
    {
      String::size_type snd(line.find('|', 3));
      String::size_type third(0);
      if (snd != String::npos)
      {
        third = line.find('|', ++snd) + 1;

        accession = line.substr(third, line.find('|', third) - third);
        accession_type = line.substr(snd, third - 1 - snd);
      }
      if (accession_type == "gb")
      {
        accession_type = "GenBank";
      }
      else if (accession_type == "emb")
      {
        accession_type = "EMBL";
      }
      else if (accession_type == "dbj")
      {
        accession_type = "DDBJ";
      }
      else if (accession_type == "ref")
      {
        accession_type = "NCBI";
      }
      else if ((accession_type == "sp") || (accession_type == "tr"))
      {
        accession_type = "SwissProt";
      }
      else if (accession_type == "gnl")
      {
        accession_type = accession;
        snd = line.find('|', third);
        third = line.find('|', ++snd);
        if (third != String::npos)
        {
          accession = line.substr(snd, third - snd);
        }
        else
        {
          third = line.find(' ', snd);
          if (third != String::npos)
          {
            accession = line.substr(snd, third - snd);
          }
          else
          {
            accession = line.substr(snd);
          }
        }
      }
      else
      {
        String::size_type pos1(line.find('(', 0));
        String::size_type pos2(0);
        if (pos1 != String::npos)
        {
          pos2 = line.find(')', ++pos1);
          if (pos2 != String::npos)
          {
            accession = line.substr(pos1, pos2 - pos1);
            if ((accession.size() == 6) && (String(swissprot_prefixes).find(accession[0], 0) != String::npos))
            {
              accession_type = "SwissProt";
            }
            else
            {
              accession.clear();
            }
          }
        }
        if (accession.empty())
        {
          accession_type = "gi";
          if (snd != String::npos)
          {
            accession = line.substr(3, snd - 4);
          }
          else
          {
            snd = line.find(' ', 3);
            if (snd != String::npos)
            {
              accession = line.substr(3, snd - 3);
            }
            else
            {
              accession = line.substr(3);
            }
          }
        }
      }
    }
    else if (line.hasPrefix("ref"))
    {
      accession = line.substr(4, line.find('|', 4) - 4);
      accession_type = "NCBI";
    }
    else if (line.hasPrefix("gnl"))
    {
      line.erase(0, 3);
      accession_type = line.substr(0, line.find('|', 0));
      accession = line.substr(accession_type.length() + 1);
    }
    else if (line.hasPrefix("lcl"))
    {
      line.erase(0, 4);
      accession_type = "lcl";
      accession = line;
    }
    else
    {
      String::size_type pos1(line.find('(', 0));
      String::size_type pos2(0);
      if (pos1 != String::npos)
      {
        pos2 = line.find(')', ++pos1);
        if (pos2 != String::npos)
        {
          accession = line.substr(pos1, pos2 - pos1);
          if ((accession.size() == 6) && (String(swissprot_prefixes).find(accession[0], 0) != String::npos))
          {
            accession_type = "SwissProt";
          }
          else
          {
            accession.clear();
          }
        }
      }
      if (accession.empty())
      {
        pos1 = line.find('|');
        accession = line.substr(0, pos1);
        if ((accession.size() == 6) && (String(swissprot_prefixes).find(accession[0], 0) != String::npos))
          accession_type = "SwissProt";
        else
        {
          pos1 = line.find(' ');
          accession = line.substr(0, pos1);
          if ((accession.size() == 6) && (String(swissprot_prefixes).find(accession[0], 0) != String::npos))
          {
            accession_type = "SwissProt";
          }
          else
          {
            accession = line.substr(0, 6);
            if (String(swissprot_prefixes).find(accession[0], 0) != String::npos)
            {
              accession_type = "SwissProt";
            }
            else
            {
              accession.clear();
            }
          }
        }
      }
    }
    if (accession.empty())
    {
      accession = line.trim();
      accession_type = "unknown";
    }
  }

  void
  InspectOutfile::getPrecursorRTandMZ(
    const vector<pair<String, vector<pair<Size, Size> > > >& files_and_peptide_identification_with_scan_number,
    vector<PeptideIdentification>& ids)
  {
    PeakMap experiment;
    String type;

    for (vector<pair<String, vector<pair<Size, Size> > > >::const_iterator fs_i = files_and_peptide_identification_with_scan_number.begin(); fs_i != files_and_peptide_identification_with_scan_number.end(); ++fs_i)
    {
      getExperiment(experiment, type, fs_i->first); // may throw an exception if the filetype could not be determined

      if (experiment.size() < fs_i->second.back().second)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough scans in file! (" + String(experiment.size()) + " available, should be at least " + String(fs_i->second.back().second) + ")", fs_i->first);
      }

      for (vector<pair<Size, Size> >::const_iterator pi_scan_i = fs_i->second.begin(); pi_scan_i != fs_i->second.end(); ++pi_scan_i)
      {
        ids[pi_scan_i->first].setMZ(experiment[pi_scan_i->second - 1].getPrecursors()[0].getMZ());
        ids[pi_scan_i->first].setRT(experiment[pi_scan_i->second - 1].getRT());
      }
    }
  }

  void
  InspectOutfile::compressTrieDB(
    const String& database_filename,
    const String& index_filename,
    vector<Size>& wanted_records,
    const String& snd_database_filename,
    const String& snd_index_filename,
    bool append)
  {
    if (database_filename == snd_database_filename)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Same filename can not be used for original and second database!", database_filename);
    }
    if (index_filename == snd_index_filename)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Same filename can not be used for original and second database!", index_filename);
    }
    ifstream database(database_filename.c_str());
    if (!database)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, database_filename);
    }

    ifstream index(index_filename.c_str(), ios::in | ios::binary);
    if (!index)
    {
      database.close();
      database.clear();
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index_filename);
    }

    // determine the length of the index file
    index.seekg(0, ios::end);
    streampos index_length = index.tellg();
    index.seekg(0, ios::beg);
    bool empty_records = wanted_records.empty();
    if (wanted_records.empty())
    {
      for (Size i = 0; i < index_length / record_length_; ++i)
        wanted_records.push_back(i);
    }

    // take the wanted records, copy their sequences to the new db and write the index file accordingly
    ofstream snd_database;
    if (append)
    {
      snd_database.open(snd_database_filename.c_str(), std::ios::out | std::ios::app);
    }
    else
    {
      snd_database.open(snd_database_filename.c_str(), std::ios::out | std::ios::trunc);
    }
    if (!snd_database)
    {
      database.close();
      database.clear();
      index.close();
      index.clear();
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, snd_database_filename);
    }

    ofstream snd_index;
    if (append)
    {
      snd_index.open(snd_index_filename.c_str(), std::ios::out | std::ios::binary | std::ios::app);
    }
    else
    {
      snd_index.open(snd_index_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    }
    if (!snd_index)
    {
      database.close();
      database.clear();
      index.close();
      index.clear();
      snd_database.close();
      snd_database.clear();
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, snd_index_filename);
    }

    char* index_record = new char[record_length_]; // to copy one record from the index file
    Size database_pos(0), snd_database_pos(0); // their sizes HAVE TO BE 4 bytes
    stringbuf sequence;

    for (vector<Size>::const_iterator wr_i = wanted_records.begin(); wr_i != wanted_records.end(); ++wr_i)
    {
      // get the according record in the index file
      if (index_length < Int((*wr_i + 1) * record_length_)) // if the file is too short
      {
        delete[] index_record;
        database.close();
        database.clear();
        index.close();
        index.clear();
        snd_database.close();
        snd_database.clear();
        snd_index.close();
        snd_index.clear();
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "index file is too short!", index_filename);
      }
      index.seekg((*wr_i) * record_length_);
      index.read(index_record, record_length_);

      // all but the first sequence are preceded by an asterisk
      if (append)
      {
        snd_database.put(trie_delimiter_);
      }
      append = true;

      // check if we have to reverse the database_pos part (which is saved in little endian)
      if (OPENMS_IS_BIG_ENDIAN)
      {
        for (Size i = 0; i < trie_db_pos_length_ / 2; i++)
        {
          char tmp = index_record[db_pos_length_ + i];
          index_record[db_pos_length_ + i] = index_record[db_pos_length_ + trie_db_pos_length_ - 1 - i];
          index_record[db_pos_length_ + trie_db_pos_length_ - 1 - i] = tmp;
        }
      }

      // go to the beginning of the sequence

      // whoever wrote this code - please don't ever do this again.
      // x86 does *not* have a monopoly, nor does little endian.
      memcpy(&database_pos, index_record + db_pos_length_, trie_db_pos_length_);
      database.seekg(database_pos);

      // store the corresponding index for the second database
      snd_database_pos = snd_database.tellp(); // get the position in the second database

      memcpy(index_record + db_pos_length_, &snd_database_pos, trie_db_pos_length_); // and copy to its place in the index record

      // fixing the above "suboptimal" code
      if (OPENMS_IS_BIG_ENDIAN)
      {
        for (Size i = 0; i < trie_db_pos_length_ / 2; i++)
        {
          char tmp = index_record[db_pos_length_ + i];
          index_record[db_pos_length_ + i] = index_record[db_pos_length_ + trie_db_pos_length_ - 1 - i];
          index_record[db_pos_length_ + trie_db_pos_length_ - 1 - i] = tmp;
        }
      }

      snd_index.write((char*) index_record, record_length_); // because only the trie-db position changed, not the position in the original database, nor the protein name

      // store the sequence
      database.get(sequence, trie_delimiter_);
      snd_database << sequence.str();
      sequence.str("");
    }


    if (empty_records)
    {
      wanted_records.clear();
    }
    delete[] index_record;
    database.close();
    database.clear();
    index.close();
    index.clear();
    snd_database.close();
    snd_database.clear();
    snd_index.close();
    snd_index.clear();
  }

  void InspectOutfile::generateTrieDB(
    const String& source_database_filename,
    const String& database_filename,
    const String& index_filename,
    bool append,
    const String& species)
  {
    ifstream source_database(source_database_filename.c_str());
    if (!source_database)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source_database_filename);
    }

    // get the labels
    String ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
    getLabels(source_database_filename, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);

    ofstream database;
    if (append)
    {
      database.open(database_filename.c_str(), ios::app | ios::out);
    }
    else
    {
      database.open(database_filename.c_str());
    }
    if (!database)
    {
      source_database.close();
      source_database.clear();
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, database_filename);
    }
    ofstream index;
    if (append)
    {
      index.open(index_filename.c_str(), ios::app | ios::out | ios::binary);
    }
    else
    {
      index.open(index_filename.c_str(), ios::out | ios::binary);
    }
    if (!index)
    {
      source_database.close();
      source_database.clear();
      database.close();
      database.clear();
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index_filename);
    }

    // using flags to mark what has already been read
    // the flags
    unsigned char ac_flag = 1;
    unsigned char species_flag = !species.empty() * 2; // if no species is given, take all proteins
    unsigned char sequence_flag = 4;
    // the value
    unsigned char record_flags = 0;

    String::size_type pos(0); // the position in a line
    unsigned long long source_database_pos = source_database.tellg(); // the start of a protein in the source database
    unsigned long long source_database_pos_buffer = 0; // because you don't know whether a new protein starts unless the line is read, the actual position is buffered before any new getline
    Size database_pos(0);
    String line, sequence, protein_name;
    char* record = new char[record_length_]; // a record in the index file
    char* protein_name_pos = record + db_pos_length_ + trie_db_pos_length_;

    while (getline(source_database, line))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      line.trim();

      // empty and comment lines are skipped
      if (line.empty() || line.hasPrefix(comment_label))
      {
        source_database_pos_buffer = source_database.tellg();
        continue;
      }

      // read the sequence if the accession and the species have been read already
      if (record_flags == (ac_flag | species_flag | sequence_flag))
      {
        if (!line.hasPrefix(sequence_end_label)) // if it is still the same protein, append the sequence
        {
          line.trim(); // erase all whitespaces from the sequence
          line.remove(trie_delimiter_);
          // save this part of the sequence
          sequence.append(line);
        }
        else // if a new protein is found, write down the old one
        {
          // if the sequence is not empty, the record has the correct form
          if (!sequence.empty())
          {
            // all but the first record in the database are preceded by an asterisk (if in append mode an asterisk has to be put at any time)
            if (append)
            {
              database.put('*');
            }
            database_pos = database.tellp();

            // write the record
            memcpy(record, &source_database_pos, db_pos_length_); // source database position
            if (OPENMS_IS_BIG_ENDIAN)
            {
              for (Size i = 0; i < db_pos_length_ / 2; i++)
              {
                char tmp = record[i];
                record[i] = record[db_pos_length_ - 1 - i];
                record[db_pos_length_ - 1 - i] = tmp;
              }
            }

            // whoever wrote this code - please don't ever do this again.
            // x86 does *not* have a monopoly, nor does little endian.
            memcpy(record + db_pos_length_, &database_pos, trie_db_pos_length_); // database position

            // fix the above "suboptimal" code
            if (OPENMS_IS_BIG_ENDIAN)
            {
              for (Size i = 0; i < trie_db_pos_length_ / 2; i++)
              {
                char tmp = record[db_pos_length_ + i];
                record[db_pos_length_ + i] = record[db_pos_length_ + trie_db_pos_length_ - 1 - i];
                record[db_pos_length_ + trie_db_pos_length_ - 1 - i] = tmp;
              }
            }

            index.write(record, record_length_);
            // protein name / accession has already been written
            database << sequence;
            source_database_pos = source_database_pos_buffer; // the position of the start of the new protein
            append = true;
          }
          sequence.clear();

          // set back the record flags for a new record
          record_flags = 0;
        }
      }

      // if not reading the sequence
      if (!(record_flags & sequence_flag))
      {
        if (line.hasPrefix(ac_label))
        {
          pos = ac_label.length(); // find the beginning of the accession

          while ((line.length() > pos) && (line[pos] < 33))
            ++pos; // discard the whitespaces after the label
          if (pos != line.length()) // if no accession is found, skip this protein
          {
            memset(protein_name_pos, 0, protein_name_length_); // clear the protein name
            // read at most protein_name_length_ characters from the record name and write them to the record
            protein_name = line.substr(pos, protein_name_length_);
            protein_name.substitute('>', '}');
            // cppcheck produces a false positive warning here -> ignore
            // cppcheck-suppress redundant copy
            memcpy(protein_name_pos, protein_name.c_str(), protein_name.length());

            record_flags |= ac_flag; // set the ac flag
          }
          else
            record_flags = 0;
        }
        // if a species line is found and an accession has already been found, check whether this record is from the wanted species, if not, skip it
        if (species_flag && line.hasPrefix(species_label) && (record_flags == ac_flag))
        {
          pos = species_label.length();
          if (line.find(species, pos) != String::npos)
          {
            record_flags |= species_flag;
          }
          else
          {
            record_flags = 0;
          }
        }
        // if the beginning of the sequence is found and accession and correct species have been found
        if (line.hasPrefix(sequence_start_label) && ((record_flags & (ac_flag | species_flag)) == (ac_flag | species_flag)))
        {
          record_flags |= sequence_flag;
        }
      }
      source_database_pos_buffer = source_database.tellg();
    }
    // source file read
    source_database.close();
    source_database.clear();

    // if the last record has no sequence end label, the sequence has to be appended nevertheless (e.g. FASTA)
    if (record_flags == (ac_flag | species_flag | sequence_flag) && !sequence.empty())
    {
      // all but the first record in the database are preceded by an asterisk (if in append mode an asterisk has to be put at any time)
      if (append)
      {
        database.put('*');
      }
      database_pos = database.tellp();

      // write the record
      // whoever wrote this code - please don't ever do this again.
      // x86 does *not* have a monopoly, nor does little endian.
      memcpy(record, &source_database_pos, db_pos_length_); // source database position
      if (OPENMS_IS_BIG_ENDIAN)
      {
        for (Size i = 0; i < db_pos_length_ / 2; i++)
        {
          char tmp = record[i];
          record[i] = record[db_pos_length_ - 1 - i];
          record[db_pos_length_ - 1 - i] = tmp;
        }
      }

      memcpy(record + db_pos_length_, &database_pos, trie_db_pos_length_); // database position

      // fix the above "suboptimal" code
      if (OPENMS_IS_BIG_ENDIAN)
      {
        for (Size i = 0; i < trie_db_pos_length_ / 2; i++)
        {
          char tmp = record[db_pos_length_ + i];
          record[db_pos_length_ + i] = record[db_pos_length_ + trie_db_pos_length_ - 1 - i];
          record[db_pos_length_ + trie_db_pos_length_ - 1 - i] = tmp;
        }
      }

      index.write(record, record_length_);
      // protein name / accession has already been written
      database << sequence;
      append = true;
    }

    delete[] record;

    // close the filestreams
    database.close();
    database.clear();
    index.close();
    index.clear();
  }

  void InspectOutfile::getLabels(
    const String& source_database_filename,
    String& ac_label,
    String& sequence_start_label,
    String& sequence_end_label,
    String& comment_label,
    String& species_label)
  {
    ac_label = sequence_start_label = sequence_end_label = comment_label = species_label = "";
    ifstream source_database(source_database_filename.c_str());
    if (!source_database)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source_database_filename);
    }

    String line;
    while (getline(source_database, line) && (sequence_start_label.empty()))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      if (line.trim().empty())
      {
        continue;
      }
      else if (line.hasPrefix(">"))
      {
        ac_label = ">";
        sequence_start_label = ">";
        sequence_end_label = ">";
        comment_label = ";";
        species_label = ">";
      }
      else if (line.hasPrefix("SQ"))
      {
        ac_label = "AC";
        sequence_start_label = "SQ";
        sequence_end_label = "//";
        comment_label = "CC";
        species_label = "OS";
      }
    }
    source_database.close();
    source_database.clear();

    // if no known start separator is found
    if (sequence_start_label.empty())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "database has unknown file format (neither trie nor FASTA nor swissprot)", source_database_filename);
    }
  }

  vector<Size> InspectOutfile::getWantedRecords(const String& result_filename, double p_value_threshold)
  {
    // check whether the p_value is correct
    if ((p_value_threshold < 0) || (p_value_threshold > 1))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the parameters 'p_value_threshold' must be >= 0 and <=1 !");
    }

    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }

    String line;
    vector<String> substrings;

    set<Size> wanted_records_set;

    vector<Size>
    wanted_records,
      corrupted_lines;

    Size line_number(0);

    // get the header
    Int
      spectrum_file_column(-1),
    scan_column(-1),
    peptide_column(-1),
    protein_column(-1),
    charge_column(-1),
    MQ_score_column(-1),
    p_value_column(-1),
    record_number_column(-1),
    DB_file_pos_column(-1),
    spec_file_pos_column(-1);

    Size number_of_columns(0);

    if (!getline(result_file, line))
    {
      result_file.close();
      result_file.clear();
      throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }
    ++line_number;
    readOutHeader(result_filename, line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns);

    while (getline(result_file, line))
    {
      ++line_number;
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      line.trim();
      if (line.empty())
      {
        continue;
      }
      line.split('\t', substrings);

      // check whether the line has enough columns
      if (substrings.size() != number_of_columns)
      {
        corrupted_lines.push_back(line_number);
        continue;
      }

      // check whether the line has enough columns
      if (substrings.size() != number_of_columns)
      {
        continue;
      }
      // take only those peptides whose p-value is less or equal the given threshold
      if (substrings[p_value_column].toFloat() > p_value_threshold)
      {
        continue;
      }
      wanted_records_set.insert(substrings[record_number_column].toInt());
    }

    result_file.close();
    result_file.clear();

    for (set<Size>::const_iterator rn_i = wanted_records_set.begin(); rn_i != wanted_records_set.end(); ++rn_i)
    {
      wanted_records.push_back(*rn_i);
    }

    return wanted_records;
  }

  bool
  InspectOutfile::getSearchEngineAndVersion(
    const String& cmd_output,
    ProteinIdentification& protein_identification)
  {
    protein_identification.setSearchEngine("InsPecT");
    protein_identification.setSearchEngineVersion("unknown");
    // searching for something like this: InsPecT version 20060907, InsPecT version 20100331
    QString response(cmd_output.toQString());
    QRegularExpression rx("InsPecT (version|vesrion) (\\d+)"); // older versions of InsPecT have typo...
    auto match = rx.match(response);
    if (!match.hasMatch())
    {
      return false;
    }
    protein_identification.setSearchEngineVersion(match.captured(2));
    return true;
  }

  void
  InspectOutfile::readOutHeader(
    const String& filename,
    const String& header_line,
    Int& spectrum_file_column,
    Int& scan_column,
    Int& peptide_column,
    Int& protein_column,
    Int& charge_column,
    Int& MQ_score_column,
    Int& p_value_column,
    Int& record_number_column,
    Int& DB_file_pos_column,
    Int& spec_file_pos_column,
    Size& number_of_columns)
  {
    spectrum_file_column = scan_column = peptide_column = protein_column = charge_column = MQ_score_column = p_value_column = record_number_column = DB_file_pos_column = spec_file_pos_column = -1;

    vector<String> substrings;
    header_line.split('\t', substrings);

    // #SpectrumFile Scan# Annotation Protein Charge MQScore Length TotalPRMScore MedianPRMScore FractionY FractionB Intensity NTT p-value F-Score DeltaScore DeltaScoreOther RecordNumber DBFilePos SpecFilePos
    for (vector<String>::const_iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i)
    {
      if ((*s_i) == "#SpectrumFile")
      {
        spectrum_file_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "Scan#")
      {
        scan_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "Annotation")
      {
        peptide_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "Protein")
      {
        protein_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "Charge")
      {
        charge_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "MQScore")
      {
        MQ_score_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "p-value")
      {
        p_value_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "RecordNumber")
      {
        record_number_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "DBFilePos")
      {
        DB_file_pos_column = s_i - substrings.begin();
      }
      else if ((*s_i) == "SpecFilePos")
      {
        spec_file_pos_column = s_i - substrings.begin();
      }
    }

    if ((spectrum_file_column == -1) || (scan_column == -1) || (peptide_column == -1) || (protein_column == -1) || (charge_column == -1) || (MQ_score_column == -1) || (p_value_column == -1) || (record_number_column == -1) || (DB_file_pos_column == -1) || (spec_file_pos_column == -1))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "at least one of the columns '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge', 'MQScore', 'p-value', 'RecordNumber', 'DBFilePos' or 'SpecFilePos' is missing!", filename);
    }
    number_of_columns = substrings.size();
  }

  const Size InspectOutfile::db_pos_length_ = 8;
  const Size InspectOutfile::trie_db_pos_length_ = 4;
  const Size InspectOutfile::protein_name_length_ = 80;
  const Size InspectOutfile::record_length_ = db_pos_length_ + trie_db_pos_length_ + protein_name_length_;
  const char InspectOutfile::trie_delimiter_ = '*';
  const String InspectOutfile::score_type_ = "Inspect";

} //namespace OpenMS

#pragma clang diagnostic pop

