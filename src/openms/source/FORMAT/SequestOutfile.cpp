// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <fstream>
#include <sstream>

using namespace std;

namespace OpenMS
{

#if 0 // useful for debugging
  template <typename ContainerType>
  void printContainer(std::ostream& os, ContainerType rhs, const String& separator = " ", const String& suffix = "\n", const String& prefix = "\n")
  {
    os << prefix;
    for (typename ContainerType::const_iterator cit = rhs.begin();; )
    {
      os << *cit;
      ++cit;
      if (cit == rhs.end())
      {
        break;
      }
      os << separator;
    }
    os << suffix;
  }

#endif

  SequestOutfile::SequestOutfile() = default;

  SequestOutfile::SequestOutfile(const SequestOutfile&) = default;

  SequestOutfile::~SequestOutfile() = default;

  SequestOutfile& SequestOutfile::operator=(const SequestOutfile& sequest_outfile)
  {
    if (this == &sequest_outfile)
      return *this;

    return *this;
  }

  bool SequestOutfile::operator==(const SequestOutfile&) const
  {
    return true;
  }

  void SequestOutfile::load(const String& result_filename,
                            vector<PeptideIdentification>& peptide_identifications,
                            ProteinIdentification& protein_identification,
                            const double p_value_threshold,
                            vector<double>& pvalues,
                            const String& database,
                            const bool ignore_proteins_per_peptide
                            )
  {

    // check whether the p_value is correct
    if ((p_value_threshold < 0) || (p_value_threshold > 1))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the parameters 'p_value_threshold' must be >= 0 and <=1 !");
    }

    // if no p_values were computed take all peptides
    bool no_pvalues = pvalues.empty();
    if (no_pvalues)
    {
      pvalues.push_back(0.0); // to make sure pvalues.end() is never reached
    }
    // generally used variables
    String
      line,
      buffer,
      sequest,
      sequest_version,
      database_type,
      identifier;

    vector<String> substrings;

    // map the protein hits according to their accession number in the result file
    map<String, Size> ac_position_map;

    // get the protein hits that have already been found in another out-file
    vector<ProteinHit> protein_hits = protein_identification.getHits();

    // and insert them Into the map
    for (const ProteinHit& phit_i : protein_hits)
    {
      ac_position_map.insert(make_pair(phit_i.getAccession(), ac_position_map.size()));
    }

    String accession, accession_type, score_type;

    DateTime datetime;
    double precursor_mz_value(0.0);
    Size
      precursor_mass_type(0),
      ion_mass_type(0),
      number_of_columns(0),
      displayed_peptides(0),
      proteins_per_peptide(0),
      line_number(0);

    Int
      charge(-1),
      number_column(-1),
      rank_sp_column(-1),
      id_column(-1),
      mh_column(-1),
      delta_cn_column(-1),
      xcorr_column(-1),
      sp_column(-1),
      sf_column(-1),
      ions_column(-1),
      reference_column(-1),
      peptide_column(-1),
      score_column(-1);

    String::size_type
      start(0),
      end(0);

    readOutHeader(result_filename, datetime, precursor_mz_value, charge,
        precursor_mass_type, ion_mass_type, displayed_peptides, sequest,
        sequest_version, database_type, number_column, rank_sp_column,
        id_column, mh_column, delta_cn_column, xcorr_column, sp_column,
        sf_column, ions_column, reference_column, peptide_column, score_column,
        number_of_columns);

    identifier = sequest + "_" + datetime.getDate();

    // set the search engine and its version and the score type
    protein_identification.setSearchEngine(sequest);
    protein_identification.setSearchEngineVersion(sequest_version);
    protein_identification.setIdentifier(identifier);
//      protein_identification.setScoreType("SEQUEST");

    // open the result
    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }

    while (getline(result_file, line)) // skip all lines until the one with '---'
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }

      line.trim();
      ++line_number;
      if (line.hasPrefix("---"))
      {
        break;
      }
    }

    PeptideIdentification peptide_identification;
    peptide_identification.setMZ(precursor_mz_value);
    peptide_identification.setIdentifier(identifier);
    peptide_identification.setSignificanceThreshold(p_value_threshold);

    vector<String> databases;
    databases.push_back(database);

    score_type = (sf_column == -1) ? "SEQUEST prelim." : "SEQUEST";
    peptide_identification.setScoreType(score_type);

    if (no_pvalues)
    {
      pvalues.insert(pvalues.end(), displayed_peptides, 0.0);
    }

    vector<double>::const_iterator p_value = pvalues.begin();

    for (Size viewed_peptides = 0; viewed_peptides < displayed_peptides; )
    {
      PeptideEvidence peptide_evidence;
      PeptideHit peptide_hit;
      ProteinHit protein_hit;

      ++line_number;
      // if less peptides were found than may be displayed, break
      if (!getline(result_file, line))
      {
        break;
      }
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      line.trim();
      if (line.empty())
      {
        continue; // skip empty lines
      }
      ++viewed_peptides;

      getColumns(line, substrings, number_of_columns, reference_column);

      // check whether the line has enough columns
      if (substrings.size() != number_of_columns)
      {
        stringstream error_message;
        error_message << "Wrong number of columns in line " << line_number << "! (" << substrings.size() << " present, should be " << number_of_columns << ")";
        result_file.close();
        result_file.clear();
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_message.str().c_str(), result_filename);
      }

      // check whether there are multiple proteins that belong to this peptide
      if (substrings[reference_column].find_last_of('+') != String::npos)
      {
        // save the number of multiple proteins
        proteins_per_peptide = substrings[reference_column].substr(substrings[reference_column].find_last_of('+')).toInt();
        // and remove this number from the String
        substrings[reference_column].resize(substrings[reference_column].find_last_of('+'));
      }
      else
      {
        proteins_per_peptide = 0;
      }

      // get the peptide information and insert it
      if (p_value != pvalues.end() && (*p_value) <= p_value_threshold)
      {
        peptide_hit.setScore(String(substrings[score_column].c_str()).toDouble());

        {
          if (rank_sp_column != (-1))
          {
            peptide_hit.setMetaValue("RankSp", substrings[rank_sp_column]);
          }
          if (id_column != (-1))
          {
            peptide_hit.setMetaValue("SequestId", atoi(substrings[id_column].c_str()));
          }
          if (mh_column != (-1))
          {
            peptide_hit.setMetaValue("MH", String(substrings[mh_column].c_str()).toDouble());
          }
          if (delta_cn_column != (-1))
          {
            peptide_hit.setMetaValue("DeltCn", String(substrings[delta_cn_column].c_str()).toDouble());
          }
          if (xcorr_column != (-1))
          {
            peptide_hit.setMetaValue("XCorr", String(substrings[xcorr_column].c_str()).toDouble());
          }
          if (sp_column != (-1))
          {
            peptide_hit.setMetaValue("Sp", String(substrings[sp_column].c_str()).toDouble());
          }
          if (sf_column != (-1))
          {
            peptide_hit.setMetaValue("Sf", String(substrings[sf_column].c_str()).toDouble());
          }
          if (ions_column != (-1))
          {
            peptide_hit.setMetaValue("Ions", substrings[ions_column]);
          }
        }

        peptide_hit.setCharge(charge);

        String sequence_with_mods = substrings[peptide_column];
        start = sequence_with_mods.find('.') + 1;
        end = sequence_with_mods.find_last_of('.');

        if (start >= 2)
        {
          peptide_evidence.setAABefore(sequence_with_mods[start - 2]);
        }

        if (end < sequence_with_mods.length() + 1)
        {
          peptide_evidence.setAAAfter(sequence_with_mods[end + 1]);
        }

        //remove modifications (small characters and everything that's not in the alphabet)
        String sequence;
        sequence_with_mods = substrings[peptide_column].substr(start, end - start);
        for (String::ConstIterator c_i = sequence_with_mods.begin(); c_i != sequence_with_mods.end(); ++c_i)
        {
          if ((bool) isalpha(*c_i) && (bool) isupper(*c_i))
          {
            sequence.append(1, *c_i);
          }
        }
        peptide_hit.setSequence(AASequence::fromString(sequence));

        peptide_hit.setRank(substrings[rank_sp_column].substr(0, substrings[rank_sp_column].find('/')).toInt());

        // get the protein information
        getACAndACType(substrings[reference_column], accession, accession_type);
        protein_hit.setAccession(accession);
//              protein_hit.setRank(ac_position_map.size());
        /// @todo simply sum up score? (Martin)

        if (ac_position_map.insert(make_pair(accession, protein_hits.size())).second)
        {
          protein_hits.push_back(protein_hit);
        }

        peptide_evidence.setProteinAccession(accession);
        peptide_hit.addPeptideEvidence(peptide_evidence);

        if (!ignore_proteins_per_peptide)
        {
          for (Size i = 0; i < proteins_per_peptide; ++i)
          {
            getline(result_file, line);
            if (!line.empty() && (line[line.length() - 1] < 33))
            {
              line.resize(line.length() - 1);
            }

            line.trim();
            // all these lines look like '0  accession', e.g. '0  gi|1584947|prf||2123446B gamma sar'
            /*if (!line.hasPrefix("0  ")) // if the line doesn't look like that
             {
             stringstream error_message;
             error_message << "Line " << line_number << " doesn't look like a line with additional found proteins! (Should look like this: 0  gi|1584947|prf||2123446B gamma sar)";
             result_file.close();
             result_file.clear();
             throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_message.str().c_str() , result_filename);
             }*/
            line.erase(0, 3);

            getACAndACType(line, accession, accession_type);
            protein_hit.setAccession(accession);
            // protein_hit.setRank(ac_position_map.size());
            // @todo simply add up score
            // Simply sum up score? (Martin)
            // protein_hit.setScore(0.0);

            if (ac_position_map.insert(make_pair(accession, protein_hits.size())).second)
            {
              protein_hits.push_back(protein_hit);
            }

            PeptideEvidence pe;
            pe.setProteinAccession(accession);
            peptide_hit.addPeptideEvidence(pe);
          }
        }

        peptide_identification.insertHit(peptide_hit);
      }
      else // if the pvalue is higher than allowed
      {
        if (!ignore_proteins_per_peptide)
        {
          for (Size i = 0; i < proteins_per_peptide; ++i)
          {
            getline(result_file, line);
          }
        }
      }

      ++p_value;
    }
    result_file.close();
    result_file.clear();

    if (no_pvalues)
    {
      pvalues.clear();
    }
    if (!peptide_identification.getHits().empty())
    {
      peptide_identifications.push_back(peptide_identification);
    }
    protein_identification.setHits(protein_hits);
    protein_identification.setDateTime(datetime);

    ac_position_map.clear();
  }

  // get the columns from a line
  bool
  SequestOutfile::getColumns(
    const String& line,
    vector<String>& substrings,
    Size number_of_columns,
    Size reference_column)
  {
    String buffer;

    if (line.empty())
    {
      return false;
    }
    line.split(' ', substrings);

    // remove any empty strings
    substrings.erase(remove(substrings.begin(), substrings.end(), ""), substrings.end());

    for (vector<String>::iterator s_i = substrings.begin(); s_i != substrings.end(); )
    {
      // if there are three columns, the middle one being a '/', they are merged
      if (s_i + 1 != substrings.end())
      {
        if (((*(s_i + 1)) == "/") && (s_i + 2 != substrings.end()))
        {
          s_i->append(*(s_i + 1));
          s_i->append(*(s_i + 2));
          substrings.erase(s_i + 2);
          substrings.erase(s_i + 1);
        }
        // if there are two columns, and the first ends with, or the second starts with a '/', they are merged
        else if ((*(s_i + 1))[0] == '/')
        {
          s_i->append(*(s_i + 1));
          substrings.erase(s_i + 1);
        }
        else if ((*(s_i))[s_i->length() - 1] == '/')
        {
          s_i->append(*(s_i + 1));
          substrings.erase(s_i + 1);
        }
        // if there are two columns and the second is a number preceded by a '+', they are merged
        else if ((*(s_i + 1))[0] == '+')
        {
          bool is_digit(true);
          for (Size i = 1; i < (s_i + 1)->length(); ++i)
          {
            is_digit &= (bool)isdigit((*(s_i + 1))[i]);
          }
          if (is_digit && ((s_i + 1)->length() - 1))
          {
            s_i->append(*(s_i + 1));
            substrings.erase(s_i + 1);
          }
          else
          {
            ++s_i;
          }
        }
        else
        {
          ++s_i;
        }
      }
      else
      {
        ++s_i;
      }
    }

    // if there are more columns than should be, there were spaces in the protein column
    for (vector<String>::iterator s_i = substrings.begin() + reference_column; substrings.size() > number_of_columns; )
    {
      s_i->append(" ");
      s_i->append(*(s_i + 1));
      substrings.erase(s_i + 1);
    }

    return true;
  }

  // retrieve the sequences
  void
  SequestOutfile::getSequences(
    const String& database_filename,
    const map<String, Size>& ac_position_map,
    vector<String>& sequences,
    vector<pair<String, Size> >& found,
    map<String, Size>& not_found)
  {
    ifstream database_file(database_filename.c_str());
    if (!database_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, database_filename);
    }

    String line, accession, accession_type, sequence;
    not_found = ac_position_map;
    map<String, Size>::iterator nf_i = not_found.end();
    while (getline(database_file, line) && !not_found.empty())
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }

      line.trim();

      // empty and comment lines are skipped
      if (line.empty() || line.hasPrefix(";"))
      {
        continue;
      }
      // the sequence belonging to the predecessing protein ('>') is stored, so
      // when a new protein ('>') is found, save the sequence of the old
      // protein
      if (line.hasPrefix(">"))
      {
        getACAndACType(line, accession, accession_type);
        if (nf_i != not_found.end())
        {
          sequences.push_back(sequence);
          found.emplace_back(*nf_i);
          not_found.erase(nf_i);
        }
        nf_i = not_found.find(accession); // for the first protein in the database, there's no predecessing protein
        sequence.clear();
      }
      else if (nf_i != not_found.end())
      {
        sequence.append(line);
      }
    }
    if (nf_i != not_found.end())
    {
      sequences.push_back(sequence);
      found.emplace_back(*nf_i);
      not_found.erase(nf_i);
    }

    database_file.close();
    database_file.clear();
  }

  void SequestOutfile::getACAndACType(String line, String& accession, String& accession_type)
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
        {
          accession_type = "SwissProt";
        }
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

  void SequestOutfile::readOutHeader(
    const String& result_filename,
    DateTime& datetime,
    double& precursor_mz_value,
    Int& charge,
    Size& precursor_mass_type,
    Size& ion_mass_type,
    Size& displayed_peptides,
    String& sequest,
    String& sequest_version,
    String& database_type,
    Int& number_column,
    Int& rank_sp_column,
    Int& id_column,
    Int& mh_column,
    Int& delta_cn_column,
    Int& xcorr_column,
    Int& sp_column,
    Int& sf_column,
// Int& P_column,
    Int& ions_column,
    Int& reference_column,
    Int& peptide_column,
    Int& score_column,
    Size& number_of_columns)
  {
    charge = 0;
    precursor_mz_value = 0.0;
    precursor_mass_type = ion_mass_type = 0;

    // open the result
    ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
    }

    String line, buffer;
    vector<String> substrings;

    // get the date and time
    DateTime datetime_empty;
    datetime.clear();
    datetime_empty.clear();

    while (getline(result_file, line))
    {
      if (!line.empty() && (line[line.length() - 1] < 33))
      {
        line.resize(line.length() - 1);
      }
      line.trim();
      line.split(',', substrings);

      if (line.hasSuffix(".out")) // next line is the sequest version
      {
        // \\bude\langwisc\temp\Inspect_Sequest.mzXML.13.1.d.out
        // TurboSEQUEST v.27 (rev. 12), (c) 1998-2005
        if (!getline(result_file, line))
        {
          result_file.close();
          result_file.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No Sequest version found!", result_filename);
        }
        if (!line.empty() && (line[line.length() - 1] < 33))
        {
          line.resize(line.length() - 1);
        }
        line.trim();

        if (line.hasSubstring(","))
        {
          line = line.substr(0, line.find(',', 0));
        }
        buffer = line;
        buffer.toUpper();
        if (!buffer.hasSubstring("SEQUEST"))
        {
          result_file.close();
          result_file.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No Sequest version found!", result_filename);
        }
        // the search engine
        String::size_type pos(0), pos1(0);
        pos1 = buffer.find("SEQUEST", 0) + strlen("SEQUEST");
        pos = line.find(' ', 0);
        if (pos == String::npos || pos >= pos1)
        {
          sequest = line.substr(0, pos1);
          pos = 0;
        }
        else
          sequest = line.substr(pos, pos1 - pos);
        // the version
        pos = line.find(' ', pos1);
        if (pos != String::npos)
        {
          pos1 = line.find(',', ++pos);
          if (pos1 == String::npos)
          {
              sequest_version = line.substr(pos);
          }
          else
          {
            sequest_version = line.substr(pos, pos1 - pos);
          }
        } // else no version was found
      }
      else if (line.hasPrefix("(M+H)+ mass = "))
      {
        line.erase(0, strlen("(M+H)+ mass = "));
        line.split(' ', substrings);
        precursor_mz_value = substrings[0].toFloat();
        charge = substrings[3].substr(1, 2).toInt();
        line = *(--substrings.end());
        line.split('/', substrings);
        substrings[0].toUpper();
        substrings[1].toUpper();
        precursor_mass_type = (substrings[0] == "MONO") ? 1 : 0;
        ion_mass_type = (substrings[1] == "MONO") ? 1 : 0;
      }
      else if ((!substrings.empty()) && (substrings[0].size() == 10) && (substrings[0][2] == '/') && (substrings[0][5] == '/'))
      {
        datetime.setDate(substrings[0]);
        buffer = substrings[1];
        buffer.trim();
        buffer.split(' ', substrings);
        // 12:00 = 12:00 PM; 24:00 = 12:00 AM
        Int hour = substrings[0].substr(0, 2).toInt();
        if ((hour == 12) && (substrings[1] == "AM"))
        {
          substrings[0].replace(0, 2, "00");
        }
        else if ((hour != 12) && (substrings[1] == "PM"))
        {
          hour += 12;
          substrings[0].replace(0, 2, String(hour));
        }
        substrings[0].append(":00");
        datetime.setTime(substrings[0]);
      }
      else if (line.hasPrefix("# bases"))
      {
          database_type = "bases";
      }
      else if (line.hasPrefix("# amino acids"))
      {
        database_type = "amino acids";
      }
      else if (line.hasPrefix("display top") && substrings[0].hasPrefix("display top")) // get the number of peptides displayed
      {
        displayed_peptides = strlen("display top ");
        displayed_peptides = substrings[0].substr(displayed_peptides, substrings[0].find('/', displayed_peptides) - displayed_peptides).toInt();
      }
      else if (line.hasPrefix("#"))
      {
        break; // the header is read
      }
    }

    if (datetime == datetime_empty)
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No time found!", result_filename);
    }
    if (sequest.empty())
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No Sequest version found!", result_filename);
    }
    if (sequest_version.empty())
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No Sequest version found!", result_filename);
    }
    if (!precursor_mz_value)
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No precursor mass found found!", result_filename);
    }
    if (!charge)
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No charge found!", result_filename);
    }

    if (!line.hasPrefix("#")) // check whether the header line was found
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not enough lines in file (no table header found)!", result_filename);
    }

    // retrieve the number of lines and get the needed column numbers
    number_column = -1;
    rank_sp_column = -1;
    id_column = -1;
    mh_column = -1;
    delta_cn_column = -1;
    xcorr_column = -1;
    sp_column = -1;
    sf_column = -1;
//      P_column = -1;
    ions_column = -1;
    reference_column = -1;
    peptide_column = -1;

    // get the single columns
    line.split(' ', substrings);

    // remove the empty strings
    for (vector<String>::iterator i = substrings.begin(); i != substrings.end(); )
    {
      i->trim();
      if (i->empty())
      {
          i = substrings.erase(i);
      }
      else
      {
          ++i;
      }
    }
    number_of_columns = substrings.size();

    // get the numbers of the columns
    for (vector<String>::const_iterator iter = substrings.begin(); iter != substrings.end(); ++iter)
    {
      if (!iter->compare("#"))
      {
        number_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Rank/Sp"))
      {
        rank_sp_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Id#"))
      {
        id_column = (iter - substrings.begin());
      }
      else if (!iter->compare("(M+H)+"))
      {
        mh_column = (iter - substrings.begin());
      }
      else if (!iter->compare("deltCn"))
      {
        delta_cn_column = (iter - substrings.begin());
      }
      else if (!iter->compare("XCorr"))
      {
        xcorr_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Sp"))
      {
        sp_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Sf"))
      {
        sf_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Ions"))
      {
        ions_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Reference"))
      {
        reference_column = (iter - substrings.begin());
      }
      else if (!iter->compare("Peptide"))
      {
        peptide_column = (iter - substrings.begin());
      }
    }
    // check whether the columns are available in the table header
    if ((number_column == -1) || (rank_sp_column == -1) || /* (id_column == -1) ||*/ (mh_column == -1) || (delta_cn_column == -1) || (xcorr_column == -1) || (sp_column == -1) || (ions_column == -1) || (reference_column == -1) || (peptide_column == -1))
    {
      result_file.close();
      result_file.clear();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "at least one of the columns '#', 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions', 'Reference' or 'Peptide' is missing!", result_filename);
    }

    score_column = (sf_column == -1) ? sp_column : sf_column;

    result_file.close();
    result_file.clear();
  }

//  void SequestOutfile::getPValuesFromOutFiles(vector< pair < String, vector< double > > >& out_filenames_and_pvalues)
//  throw (Exception::FileNotFound&, Exception::ParseError)
//  {
//      DateTime datetime;
//      double
//          precursor_mz_value(0),
//          discriminant_score,
//          xcorr,
//          rank_sp,
//          delta_mass;
//
//      Size
//          precursor_mass_type(0),
//          ion_mass_type(0),
//          number_of_columns(0),
//          displayed_peptides(0),
//          line_number(0),
//          proteins_per_peptide(0),
//          peptide_length(0);
//
//      Int
//          charge(0),
//          number_column(0),
//          rank_sp_column(0),
//          id_column(0),
//          mh_column(0),
//          delta_cn_column(0),
//          xcorr_column(0),
//          sp_column(0),
//          sf_column(0),
//          ions_column(0),
//          reference_column(0),
//          peptide_column(0),
//          score_column(0);
//
//      String
//          line,
//          sequence,
//          buffer,
//          out_filename,
//          database_type;
//
//      vector< String > substrings;
//      vector< double >
//          delta_cns,
//          current_discriminant_scores,
//          pvalues;
//
// //       map< String, vector< double > > out_filenames_and_discriminant_scores;
//      vector< vector< double > > discriminant_scores;
//      map< double, Size > discriminant_scores_histogram;
//
//      for ( vector< pair < String, vector< double > > >::const_iterator fp_i = out_filenames_and_pvalues.begin(); fp_i != out_filenames_and_pvalues.end(); ++fp_i )
//      {
//          current_discriminant_scores.clear();
//          readOutHeader(fp_i->first, datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, displayed_peptides, line, line, database_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns);
//
//          // the charge is allowed from 1 to 3 only
//          if ( charge < 0 ) charge *= -1;
//          if ( charge > 3 ) charge = 3;
//
//          // reopen the result file
//          ifstream out_file(out_filename.c_str());
//          if ( !out_file )
//          {
//              throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_filename);
//          }
//
//          while ( getline(out_file, line) ) // skip all lines until the one with '---'
//          {
//              if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
//              line.trim();
//              ++line_number;
//              if ( line.hasPrefix("---") ) break;
//          }
//
//          // needed: XCorr, peptide length, delta Cn, rankSp, delta Mass
//          for ( Size viewed_peptides = 0 ; viewed_peptides < displayed_peptides; )
//          {
//              if ( !getline(out_file, line) ) break; // if fewer peptides were found than may be displayed, break
//              ++line_number;
//              if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
//              line.trim();
//              if ( line.empty() ) continue;
//
//              getColumns(line, substrings, number_of_columns, reference_column);
//              ++viewed_peptides;
//
//              // check whether the line has enough columns
//              if (substrings.size() < number_of_columns )
//              {
//                  stringstream error_message;
//                  error_message << "Wrong number of columns in line " << line_number << "! (" << substrings.size() << " present, should be " << number_of_columns << ")";
//                  throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_message.str().c_str() , out_filename);
//              }
//              delta_cns.push_back(substrings[delta_cn_column].toFloat());
//              xcorr = substrings[xcorr_column].toFloat();
//              rank_sp = substrings[rank_sp_column].toFloat();
//              delta_mass = precursor_mz_value - substrings[mh_column].toFloat();
//              buffer = substrings[peptide_column].substr(2, substrings[peptide_column].length() - 4);
//              // remove all ptms
//              for ( String::ConstIterator c_i = buffer.begin(); c_i != buffer.end(); ++c_i )
//              {
//                  if ( (bool) isalpha(*c_i) && (bool) isupper(*c_i) ) sequence.append(1, *c_i);
//              }
//
//              // compute the discriminant score
//              peptide_length = min(max_pep_lens_[charge], sequence.length());
//              discriminant_score = xcorr_weights_[charge] * (log(xcorr) / log(peptide_length * num_frags_[charge]));
//              discriminant_score += rank_sp_weights_[charge] * log(rank_sp);
//              discriminant_score += delta_mass_weights_[charge] * abs(delta_mass);
//              discriminant_score += const_weights_[charge];
//              current_discriminant_scores.push_back(discriminant_score);
//
//              // if there are multiple proteins that belong to this peptide, skip these lines
//              if ( substrings[reference_column].find_last_of('+') != String::npos )
//              {
//                  proteins_per_peptide = substrings[reference_column].substr(substrings[reference_column].find_last_of('+')).toInt();
//                  for ( Size prot = 0; prot < proteins_per_peptide; ++prot ) getline(out_file, line);
//                  line_number += proteins_per_peptide;
//              }
//          }
//
//          // close and clear the stream for further use
//          out_file.close();
//          out_file.clear();
//
//          // if only one delta cn is found, it is set to 1
//
//          if ( delta_cns.size() == 1 ) current_discriminant_scores.back() += delta_cn_weights_[charge];
//          else if ( delta_cns.size() > 1 )
//          {
//              // the delta cns are recalculated and the discriminant scores are calculated correspondingly and added to the histogram
//              vector< double >::iterator ds_i = current_discriminant_scores.begin();
//              for ( vector< double >::const_iterator dcn_i = delta_cns.begin(); dcn_i != delta_cns.end(); ++dcn_i, ++ds_i )
//              {
//                  (*ds_i) += delta_cn_weights_[charge] * (delta_cns.back() - (*dcn_i));
//                  ++discriminant_scores_histogram[*ds_i]; // bucketing; not yet finished
//              }
//          }
//          // append the discriminant scores
//          discriminant_scores.push_back(current_discriminant_scores);
// //           out_filenames_and_discriminant_scores[out_filename] = current_discriminant_scores;
//      }
//
//      // now the p-values can be computed
//      // fit two normal distributions to the data
//      Math::BasicStatistics< >
//          correct,
//          incorrect;
//          // unfinished;
//      correct.setMean();
//      correct.setVariance();
//      incorrect.setMean();
//      incorrect.setVariance();
//
// //       for ( map< String, vector< double > >::const_iterator fnds_i = out_filenames_and_discriminant_scores.begin(); fnds_i != out_filenames_and_discriminant_scores.end(); ++fnds_i )
//      vector< vector< double >::const_iterator dss_i = discriminant_scores.begin();
//      for ( vector< pair < String, vector< double > > >::iterator fp_i = out_filenames_and_pvalues.begin(); fp_i != out_filenames_and_pvalues.end(); ++fp_i, ++dss_i )
//      {
//          pvalues.clear();
// //           for ( vector< double >::const_iterator ds_i = fnds_i->second.begin(); ds_i != fnds_i->second.begin(); ++ds_i )
//          for ( vector< double >::const_iterator ds_i = dss_i->begin(); ds_i != dss_i->end(); ++ds_i )
//          {
//              pvalues.push_back(correct.normalDensity(*ds_i) / (correct.normalDensity(*ds_i) + incorrect.normalDensity(*ds_i)));
// //               p_correct = exp(-0.5 * pow((*ds_i - mean_correct) / sd, 2)) / (sd_correct * sqrt(2 * pi) );
// //               p_incorrect = exp(-0.5 * pow((*ds_i - mean_incorrect) / sd, 2)) / (sd_incorrect * sqrt(2 * pi) );
// //               pvalues.push_back();
//          }
//          fp_i->second = pvalues;
//      }
//  }

  double SequestOutfile::const_weights_[] = {0.646f, -0.959f, -1.460f};
  double SequestOutfile::xcorr_weights_[] = {5.49f, 8.362f, 9.933f};
  double SequestOutfile::delta_cn_weights_[] = {4.643f, 7.386f, 11.149f};
  double SequestOutfile::rank_sp_weights_[] = {-0.455f, -0.194f, -0.201f};
  double SequestOutfile::delta_mass_weights_[] =  {-0.84f, -0.314f, -0.277f};
  Size SequestOutfile::max_pep_lens_[] = {100, 15, 25};
  Size SequestOutfile::num_frags_[] = {2, 2, 4};
} //namespace OpenMS
