// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectInfile.h>
#include "SearchEngineFileHelpers.h"
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <sstream>

using namespace std;

namespace OpenMS
{

  // default constructor
  InspectInfile::InspectInfile() :
    modifications_per_peptide_(-1),
    blind_(2),
    maxptmsize_(-1.0),
    precursor_mass_tolerance_(-1.0),
    peak_mass_tolerance_(-1.0),
    multicharge_(2),
    tag_count_(-1)
  {
  }

  // copy constructor
  InspectInfile::InspectInfile(const InspectInfile& inspect_infile) :
    spectra_(inspect_infile.getSpectra()),
    enzyme_(inspect_infile.getEnzyme()),
    modifications_per_peptide_(inspect_infile.getModificationsPerPeptide()),
    blind_(inspect_infile.getBlind()),
    maxptmsize_(inspect_infile.getMaxPTMsize()),
    precursor_mass_tolerance_(inspect_infile.getPrecursorMassTolerance()),
    peak_mass_tolerance_(inspect_infile.getPeakMassTolerance()),
    multicharge_(inspect_infile.getMulticharge()),
    instrument_(inspect_infile.getInstrument()),
    tag_count_(inspect_infile.getTagCount()),
    PTMname_residues_mass_type_(inspect_infile.getModifications())
  {
  }

  // destructor
  InspectInfile::~InspectInfile()
  {
    PTMname_residues_mass_type_.clear();
  }

  // assignment operator
  InspectInfile& InspectInfile::operator=(const InspectInfile& inspect_infile)
  {
    if (this != &inspect_infile)
    {
      spectra_ = inspect_infile.getSpectra();
      enzyme_ = inspect_infile.getEnzyme();
      modifications_per_peptide_ = inspect_infile.getModificationsPerPeptide();
      blind_ = inspect_infile.getBlind();
      maxptmsize_ = inspect_infile.getMaxPTMsize();
      precursor_mass_tolerance_ = inspect_infile.getPrecursorMassTolerance();
      peak_mass_tolerance_ = inspect_infile.getPeakMassTolerance();
      multicharge_ = inspect_infile.getMulticharge();
      instrument_ = inspect_infile.getInstrument();
      tag_count_ = inspect_infile.getTagCount();
      PTMname_residues_mass_type_ = inspect_infile.getModifications();
    }
    return *this;
  }

  // equality operator
  bool InspectInfile::operator==(const InspectInfile& inspect_infile) const
  {
    if (this != &inspect_infile)
    {
      bool equal = true;
      equal &= (spectra_ == inspect_infile.getSpectra());
      equal &= (enzyme_ == inspect_infile.getEnzyme());
      equal &= (modifications_per_peptide_ == inspect_infile.getModificationsPerPeptide());
      equal &= (blind_ == inspect_infile.getBlind());
      equal &= (maxptmsize_ == inspect_infile.getMaxPTMsize());
      equal &= (precursor_mass_tolerance_ == inspect_infile.getPrecursorMassTolerance());
      equal &= (peak_mass_tolerance_ == inspect_infile.getPeakMassTolerance());
      equal &= (multicharge_ == inspect_infile.getMulticharge());
      equal &= (instrument_ == inspect_infile.getInstrument());
      equal &= (tag_count_ == inspect_infile.getTagCount());
      equal &= (PTMname_residues_mass_type_ == inspect_infile.getModifications());
      return equal;
    }
    return true;
  }

  void InspectInfile::store(const std::string& filename)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::TSV))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::TSV) + "'");
    }

    ofstream ofs(filename.c_str());
    if (!ofs)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    stringstream file_content;

    file_content << "spectra," << spectra_ << "\n";

    if (!db_.empty())
    {
      file_content << "db," << db_ << "\n";
    }
    if (!enzyme_.empty())
    {
      file_content << "protease," << enzyme_ << "\n";
    }
    if (blind_ != 2)
    {
      file_content << "blind," << blind_ << "\n";
    }
    //mod,+57,C,fix,carbamidomethylation
    for (std::map<std::string, vector<std::string> >::iterator mods_i = PTMname_residues_mass_type_.begin(); mods_i != PTMname_residues_mass_type_.end(); ++mods_i)
    {
      // fix", "cterminal", "nterminal", and "opt
      StringUtils::toLower(mods_i->second[2]);
      if (StringUtils::hasSuffix(mods_i->second[2], "term"))
      {
        mods_i->second[2].append("inal");
      }
      file_content << "mod," << mods_i->second[1] << "," << mods_i->second[0] << "," << mods_i->second[2] << "," << mods_i->first << "\n";
    }

    if (modifications_per_peptide_ > -1)
    {
      file_content << "mods," << modifications_per_peptide_ << "\n";
    }
    if (maxptmsize_ >= 0)
    {
      file_content << "maxptmsize," << maxptmsize_ << "\n";
    }
    if (precursor_mass_tolerance_ >= 0)
    {
      file_content << "PM_tolerance," << precursor_mass_tolerance_ << "\n";
    }
    if (peak_mass_tolerance_ >= 0)
    {
      file_content << "IonTolerance," << peak_mass_tolerance_ << "\n";
    }
    if (multicharge_ != 2)
    {
      file_content << "multicharge," << multicharge_ << "\n";
    }
    if (!instrument_.empty())
    {
      file_content << "instrument," << instrument_ << "\n";
    }
    if (tag_count_ >= 0)
    {
      file_content << "TagCount," << tag_count_ << "\n";
    }
    ofs << file_content.str();

    ofs.close();
    ofs.clear();
  }

  void InspectInfile::handlePTMs(const std::string& modification_line, const std::string& modifications_filename, const bool monoisotopic)
  {
    Internal::SearchEngineFileHelpers::handlePTMs(modification_line, modifications_filename, monoisotopic, PTMname_residues_mass_type_, OPENMS_PRETTY_FUNCTION);
  }

  const std::map<std::string, vector<std::string> >& InspectInfile::getModifications() const
  {
    return PTMname_residues_mass_type_;
  }

  const std::string& InspectInfile::getSpectra() const
  {
    return spectra_;
  }

  void InspectInfile::setSpectra(const std::string& spectra)
  {
    spectra_ = spectra;
  }

  const std::string& InspectInfile::getDb() const
  {
    return db_;
  }

  void InspectInfile::setDb(const std::string& db)
  {
    db_ = db;
  }

  const std::string& InspectInfile::getEnzyme() const
  {
    return enzyme_;
  }

  void InspectInfile::setEnzyme(const std::string& enzyme)
  {
    enzyme_ = enzyme;
  }

  Int InspectInfile::getModificationsPerPeptide() const
  {
    return modifications_per_peptide_;
  }

  void InspectInfile::setModificationsPerPeptide(Int modifications_per_peptide)
  {
    modifications_per_peptide_ = modifications_per_peptide;
  }

  UInt InspectInfile::getBlind() const
  {
    return blind_;
  }

  void InspectInfile::setBlind(UInt blind)
  {
    blind_ = blind;
  }

  float InspectInfile::getMaxPTMsize() const
  {
    return maxptmsize_;
  }

  void InspectInfile::setMaxPTMsize(float maxptmsize)
  {
    maxptmsize_ = maxptmsize;
  }

  float InspectInfile::getPrecursorMassTolerance() const
  {
    return precursor_mass_tolerance_;
  }

  void InspectInfile::setPrecursorMassTolerance(float precursor_mass_tolerance)
  {
    precursor_mass_tolerance_ = precursor_mass_tolerance;
  }

  float InspectInfile::getPeakMassTolerance() const
  {
    return peak_mass_tolerance_;
  }

  void InspectInfile::setPeakMassTolerance(float ion_tolerance)
  {
    peak_mass_tolerance_ = ion_tolerance;
  }

  UInt InspectInfile::getMulticharge() const
  {
    return multicharge_;
  }

  void InspectInfile::setMulticharge(UInt multicharge)
  {
    multicharge_ = multicharge;
  }

  const std::string& InspectInfile::getInstrument() const
  {
    return instrument_;
  }

  void InspectInfile::setInstrument(const std::string& instrument)
  {
    instrument_ = instrument;
  }

  Int InspectInfile::getTagCount() const
  {
    return tag_count_;
  }

  void InspectInfile::setTagCount(Int tag_count)
  {
    tag_count_ = tag_count;
  }

} // namespace OpenMS
