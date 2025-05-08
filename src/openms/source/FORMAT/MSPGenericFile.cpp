// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <boost/regex.hpp>
#include <fstream>
#include <array>

namespace OpenMS
{
  MSPGenericFile::MSPGenericFile() :
    DefaultParamHandler("MSPGenericFile")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  MSPGenericFile::MSPGenericFile(const String& filename, MSExperiment& library) :
    DefaultParamHandler("MSPGenericFile")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
    load(filename, library);
  }

  void MSPGenericFile::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("synonyms_separator", "|", "The character that will separate the synonyms in the Synon metaValue.");
  }

  void MSPGenericFile::updateMembers_()
  {
    synonyms_separator_ = param_.getValue("synonyms_separator").toString();
  }

  void MSPGenericFile::load(const String& filename, MSExperiment& library)
  {
    loaded_spectra_names_.clear();
    synonyms_.clear();
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const Size BUFSIZE { 65536 };
    char line[BUFSIZE];
    library.clear(true);
    MSSpectrum spectrum;
    spectrum.setMetaValue("is_valid", 0); // to avoid adding invalid spectra to the library

    boost::cmatch m;
    boost::regex re_name("(?:^Name|^NAME): (.+)", boost::regex::no_mod_s);
    boost::regex re_retention_time("(?:^Retention Time|^RETENTIONTIME): (.+)", boost::regex::no_mod_s);
    boost::regex re_synon("^synon(?:yms?)?: (.+)", boost::regex::no_mod_s | boost::regex::icase);
    boost::regex re_points_line(R"(^\d)");
    boost::regex re_point(R"((\d+(?:\.\d+)?)[: \t](\d+(?:\.\d+)?);? ?)");
    boost::regex re_cas_nist(R"(^CAS#: ([\d-]+);  NIST#: (\d+))"); // specific to NIST db
    boost::regex re_precursor_mz("^PRECURSORMZ: (.+)");
    boost::regex re_num_peaks("(?:^Num Peaks|^Num peaks): (\\d+)");
    // regex for meta values required for MetaboliteSpectralMatcher
    boost::regex re_inchi("^INCHIKEY: (.+)");
    boost::regex re_smiles("^SMILES: (.+)");
    boost::regex re_sum_formula("^FORMULA: (.+)");
    boost::regex re_precursor_type("^PRECURSORTYPE: (.+)");
    // matches everything else
    boost::regex re_metadatum("^(.+): (.+)", boost::regex::no_mod_s);
    OPENMS_LOG_INFO << "\nLoading spectra from .msp file. Please wait." << std::endl;

    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      // Peaks
      if (boost::regex_search(line, m, re_points_line))
      {
        OPENMS_LOG_DEBUG << "re_points_line\n";
        boost::regex_search(line, m, re_point);
        do
        {
          OPENMS_LOG_DEBUG << "{" << m[1] << "} {" << m[2] << "}; ";
          const double position { std::stod(m[1]) };
          const double intensity { std::stod(m[2]) };
          spectrum.push_back( Peak1D(position, intensity) );
        } while ( boost::regex_search(m[0].second, m, re_point) );
      }
      // Synon
      else if (boost::regex_search(line, m, re_synon))
      {
        // OPENMS_LOG_DEBUG << "Synon: " << m[1] << "\n";
        synonyms_.emplace_back(m[1]);
      }
      // Name
      else if (boost::regex_search(line, m, re_name))
      {
        addSpectrumToLibrary(spectrum, library);
        // OPENMS_LOG_DEBUG << "\n\nName: " << m[1] << "\n";
        spectrum.clear(true);
        synonyms_.clear();
        spectrum.setName(String(m[1]));
        spectrum.setMetaValue(Constants::UserParam::MSM_METABOLITE_NAME, spectrum.getName());
        spectrum.setMetaValue("is_valid", 1);
      }
      // Number of Peaks
      else if (boost::regex_search(line, m, re_num_peaks))
      {
        spectrum.setMetaValue("Num Peaks", String(m[1]));
      }
      // Retention Time
      else if (boost::regex_search(line, m, re_retention_time))
      {
        spectrum.setRT(std::stod(m[1]));
      }
      // set Precursor MZ
      else if (boost::regex_search(line, m, re_precursor_mz))
      {
        std::vector<Precursor> precursors;
        Precursor p;
        p.setMZ(std::stod(m[1]));
        precursors.push_back(p);
        spectrum.setPrecursors(precursors);
      }
      //CAS# NIST#
      else if (boost::regex_search(line, m, re_cas_nist))
      {
        // OPENMS_LOG_DEBUG << "CAS#: " << m[1] << "; NIST#: " << m[2] << "\n";
        spectrum.setMetaValue(String("CAS#"), String(m[1]));
        spectrum.setMetaValue(String("NIST#"), String(m[2]));
      }
      // Meta values for MetaboliteSpectralMatcher
      else if (boost::regex_search(line, m, re_inchi))
      {
        spectrum.setMetaValue(Constants::UserParam::MSM_INCHI_STRING, String(m[1]));
      }
      else if (boost::regex_search(line, m, re_smiles))
      {
        spectrum.setMetaValue(Constants::UserParam::MSM_SMILES_STRING, String(m[1]));
      }
      else if (boost::regex_search(line, m, re_sum_formula))
      {
        spectrum.setMetaValue(Constants::UserParam::MSM_SUM_FORMULA, String(m[1]));
      }
      else if (boost::regex_search(line, m, re_precursor_type))
      {
        spectrum.setMetaValue(Constants::UserParam::MSM_PRECURSOR_ADDUCT, String(m[1]));
      }      
      // Other metadata, needs to be last, matches everything
      else if (boost::regex_search(line, m, re_metadatum))
      {
        // OPENMS_LOG_DEBUG << m[1] << m[2] << "\n";
        spectrum.setMetaValue(String(m[1]), String(m[2]));
      }
    }
    // To make sure a spectrum is added even if no empty line is present before EOF
    addSpectrumToLibrary(spectrum, library);
    OPENMS_LOG_INFO << "Loading spectra from .msp file completed." << std::endl;
  }

  void MSPGenericFile::store(const String& filename, const MSExperiment& library) const
  {
    std::ofstream output_file(filename.c_str());

    // checking if file is writable
    if (!File::writable(filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    for (const auto& spectrum : library.getSpectra())
    {
      if (spectrum.getName().empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "The current spectrum misses the Name information.");
      }
      if (spectrum.size())// we will not store spectrum with no peaks
      {
        output_file << "Name: " << spectrum.getName() << '\n';
        output_file << "Retention Time: " << spectrum.getRT() << '\n';
        const auto& synonyms = spectrum.getMetaValue("Synon");
        if (synonyms.valueType() == OpenMS::DataValue::DataType::STRING_VALUE)
        {
          StringList list;
          synonyms.toString().split(synonyms_separator_, list);
          for (const auto& syn : list)
          {
            output_file << "Synon: " << syn << '\n';
          }
        }
        if (spectrum.metaValueExists("CAS#") && spectrum.metaValueExists("NIST#"))
        {
          output_file << "CAS#: " << spectrum.getMetaValue("CAS#") << ";  NIST#: " << spectrum.getMetaValue("NIST#") << '\n';
        }
        // Other metadata
        static const std::array<std::string, 4> ignore_metadata = {"Synon", "CAS#", "NIST#", "Num Peaks"};
        std::vector<String> keys;
        spectrum.getKeys(keys);
        for (const auto& key : keys)
        {
          const auto& value = spectrum.getMetaValue(key);
          if (std::find(ignore_metadata.begin(), ignore_metadata.end(), key) == ignore_metadata.end())
          {
            output_file << key << ": " << value << '\n';
          }
        }
        // Peaks
        output_file << "Num Peaks: " << spectrum.size() << '\n';
        int peak_counter = 0;
        for (const auto& peak : spectrum)
        {
          output_file << peak.getPos() << ":" << peak.getIntensity() << " ";
          if ((++peak_counter % 5) == 0)
          {
            output_file << '\n';
          }
        }
        if ((peak_counter % 5) != 0)
        {
          output_file << '\n';
        }
        // separator
        output_file << '\n';
      }
    }

    output_file.close();
  }

  void MSPGenericFile::addSpectrumToLibrary(
    MSSpectrum& spectrum,
    MSExperiment& library
  )
  {
    if (static_cast<int>(spectrum.getMetaValue("is_valid")) == 0)
    {
      return;
    }
    // Check that required metadata (Name, Num Peaks) is present
    // Num Peaks is checked later in the code (when verifying for the number of points parsed)
    if (spectrum.getName().empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "The current spectrum misses the Name information.");
    }

    // Check that the spectrum is not a duplicate (i.e. already present in `library`)
    const Size name_found = loaded_spectra_names_.count(spectrum.getName());

    if (!name_found)
    {
      // Check that all expected points are parsed
      if (!spectrum.metaValueExists("Num Peaks"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The current spectrum misses the Num Peaks information.");
      }
      const String& num_peaks { spectrum.getMetaValue("Num Peaks") };
      if (spectrum.size() != std::stoul(num_peaks) )
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          num_peaks,
          "The number of points parsed does not coincide with `Num Peaks`.");
      }

      if (!synonyms_.empty())
      {
        String synon;
        for (const String& s : synonyms_)
        {
          synon += s + synonyms_separator_;
        }
        if (!synon.empty())
        {
          synon.pop_back();
        }
        spectrum.setMetaValue("Synon", synon);
      }

      spectrum.removeMetaValue("is_valid");
      if (spectrum.getRT() < 0) 
      { // set RT to spectrum index
        spectrum.setRT(library.getSpectra().size());
      }
      library.addSpectrum(spectrum);
      loaded_spectra_names_.insert(spectrum.getName());

      if (loaded_spectra_names_.size() % 20000 == 0)
      {
        OPENMS_LOG_INFO << "Loaded " << loaded_spectra_names_.size() << " spectra..." << std::endl;
      }
    }
    else
    {
      OPENMS_LOG_INFO << "DUPLICATE: " << spectrum.getName() << std::endl;
    }

    spectrum.setMetaValue("is_valid", 0);
  }
}
