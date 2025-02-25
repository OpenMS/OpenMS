// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>

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

  void InspectInfile::store(const String& filename)
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
    for (std::map<String, vector<String> >::iterator mods_i = PTMname_residues_mass_type_.begin(); mods_i != PTMname_residues_mass_type_.end(); ++mods_i)
    {
      // fix", "cterminal", "nterminal", and "opt
      mods_i->second[2].toLower();
      if (mods_i->second[2].hasSuffix("term"))
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

  void InspectInfile::handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic)
  {
    PTMname_residues_mass_type_.clear();
    // to store the information about modifications from the ptm xml file
    std::map<String, pair<String, String> > ptm_informations;
    if (!modification_line.empty()) // if modifications are used look whether whether composition and residues (and type and name) is given, the name (type) is used (then the modifications file is needed) or only the mass and residues (and type and name) is given
    {
      vector<String> modifications, mod_parts;
      modification_line.split(':', modifications); // get the single modifications
      if (modifications.empty())
      {
        modifications.push_back(modification_line);
      }
      // to get masses from a formula
      EmpiricalFormula add_formula, substract_formula;

      String types = "OPT#FIX#";
      String name, residues, mass, type;

      // 0 - mass; 1 - composition; 2 - ptm name
      Int mass_or_composition_or_name(-1);

      for (vector<String>::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i)
      {
        if (mod_i->empty())
        {
          continue;
        }
        // clear the formulae
        add_formula = substract_formula = EmpiricalFormula();
        name = residues = mass = type = "";

        // get the single parts of the modification string
        mod_i->split(',', mod_parts);
        if (mod_parts.empty())
          mod_parts.push_back(*mod_i);
        mass_or_composition_or_name = -1;

        // check whether the first part is a mass, composition or name

        // check whether it is a mass
        try
        {
          mass = mod_parts.front();
          // to check whether the first part is a mass, it is converted into a float and then back into a string and compared to the given string
          // remove + signs because they don't appear in a float
          if (mass.hasPrefix("+"))
          {
            mass.erase(0, 1);
          }
          if (mass.hasSuffix("+"))
          {
            mass.erase(mass.length() - 1, 1);
          }
          if (mass.hasSuffix("-")) // a - sign at the end will not be converted
          {
            mass.erase(mass.length() - 1, 1);
            mass.insert(0, "-");
          }
          // if it is a mass
          if (!String(mass.toFloat()).empty()) // just check if conversion does not throw, i.e. consumes the whole string
          {
            mass_or_composition_or_name = 0;
          }
        }
        catch (Exception::ConversionError& /*c_e*/)
        {
          mass_or_composition_or_name = -1;
        }

        // check whether it is a name (look it up in the corresponding file)
        if (mass_or_composition_or_name == -1)
        {
          if (ptm_informations.empty()) // if the ptm xml file has not been read yet, read it
          {
            if (!File::exists(modifications_filename))
            {
              throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, modifications_filename);
            }
            if (!File::readable(modifications_filename))
            {
              throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, modifications_filename);
            }

            // getting all available modifications from a file
            PTMXMLFile().load(modifications_filename, ptm_informations);
          }
          // if the modification cannot be found
          if (ptm_informations.find(mod_parts.front()) != ptm_informations.end())
          {
            mass = ptm_informations[mod_parts.front()].first; // composition
            residues = ptm_informations[mod_parts.front()].second; // residues
            name = mod_parts.front(); // name

            mass_or_composition_or_name = 2;
          }
        }

        // check whether it's an empirical formula / if a composition was given, get the mass
        if (mass_or_composition_or_name == -1)
        {
          mass = mod_parts.front();
        }
        if (mass_or_composition_or_name == -1 || mass_or_composition_or_name == 2)
        {
          // check whether there is a positive and a negative formula
          String::size_type pos = mass.find("-");
          try
          {
            if (pos != String::npos)
            {
              add_formula = EmpiricalFormula(mass.substr(0, pos));
              substract_formula = EmpiricalFormula(mass.substr(++pos));
            }
            else
            {
              add_formula = EmpiricalFormula(mass);
            }
            // sum up the masses
            if (monoisotopic)
            {
              mass = String(add_formula.getMonoWeight() - substract_formula.getMonoWeight());
            }
            else
            {
              mass = String(add_formula.getAverageWeight() - substract_formula.getAverageWeight());
            }
            if (mass_or_composition_or_name == -1)
            {
              mass_or_composition_or_name = 1;
            }
          }
          catch (Exception::ParseError& /*pe*/)
          {
            PTMname_residues_mass_type_.clear();
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, *mod_i, "There's something wrong with this modification. Aborting!");
          }
        }

        // now get the residues
        mod_parts.erase(mod_parts.begin());
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            PTMname_residues_mass_type_.clear();
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, *mod_i, "No residues for modification given. Aborting!");
          }

          // get the residues
          residues = mod_parts.front();
          residues.substitute('*', 'X');
          residues.toUpper();
          mod_parts.erase(mod_parts.begin());
        }

        // get the type
        if (mod_parts.empty())
        {
          type = "OPT";
        }
        else
        {
          type = mod_parts.front();
          type.toUpper();
          if (types.find(type) != String::npos)
          {
            mod_parts.erase(mod_parts.begin());
          }
          else
          {
            type = "OPT";
          }
        }

        if (mod_parts.size() > 1)
        {
          PTMname_residues_mass_type_.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, *mod_i, "There's something wrong with the type of this modification. Aborting!");
        }

        // get the name
        if (mass_or_composition_or_name < 2)
        {
          if (mod_parts.empty())
          {
            name = "PTM_" + String(PTMname_residues_mass_type_.size());
          }
          else
          {
            name = mod_parts.front();
          }
        }

        // insert the modification
        if (PTMname_residues_mass_type_.find(name) == PTMname_residues_mass_type_.end())
        {
          PTMname_residues_mass_type_[name] = vector<String>(3);
          PTMname_residues_mass_type_[name][0] = residues;
          // mass must not have more than 5 digits after the . (otherwise the test may fail)
          PTMname_residues_mass_type_[name][1] = mass.substr(0, mass.find(".") + 6);
          PTMname_residues_mass_type_[name][2] = type;
        }
        else
        {
          PTMname_residues_mass_type_.clear();
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, *mod_i, "There's already a modification with this name. Aborting!");
        }
      }
    }
  }

  const std::map<String, vector<String> >& InspectInfile::getModifications() const
  {
    return PTMname_residues_mass_type_;
  }

  const String& InspectInfile::getSpectra() const
  {
    return spectra_;
  }

  void InspectInfile::setSpectra(const String& spectra)
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

  const String& InspectInfile::getEnzyme() const
  {
    return enzyme_;
  }

  void InspectInfile::setEnzyme(const String& enzyme)
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

  const String& InspectInfile::getInstrument() const
  {
    return instrument_;
  }

  void InspectInfile::setInstrument(const String& instrument)
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
