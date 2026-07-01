// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <vector>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  /**
    @brief Read/write Mascot generic files (MGF).

    For details of the format, see http://www.matrixscience.com/help/data_file_help.html#GEN.

    In addition to the core MGF fields (TITLE, PEPMASS, CHARGE, RTINSECONDS, SCANS, MSLEVEL)
    and the GNPS/metabolomics extensions (NAME, COMPOUND_NAME, INCHI, SMILES, IONMODE, ...),
    this reader/writer also supports the sequence-query field SEQ. SEQ lines are exposed on
    the MSSpectrum via the "SEQ" meta value as a StringList (always, even for a single SEQ)
    and written back out as one SEQ= line per entry.

    @htmlinclude OpenMS_MascotGenericFile.parameters

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MascotGenericFile :
    public ProgressLogger,
    public DefaultParamHandler
  {
public:

    /// constructor
    MascotGenericFile();

    /// destructor
    ~MascotGenericFile() override;

    /// docu in base class
    void updateMembers_() override;

    /// stores the experiment data in a MascotGenericFile that can be used as input for MASCOT shell execution (optionally a compact format is used: no zero-intensity peaks, limited number of decimal places)
    void store(const std::string& filename, const PeakMap& experiment,
               bool compact = false);

    /// store the experiment data in a MascotGenericFile; the output is written to the given stream, the filename will be noted in the file (optionally a compact format is used: no zero-intensity peaks, limited number of decimal places)
    void store(std::ostream& os, const std::string& filename,
               const PeakMap& experiment, bool compact = false);

    /**
      @brief loads a Mascot Generic File into a PeakMap

      @param[in] filename file name which the map should be read from
      @param[out] exp the map which is filled with the data from the given file
      @throw FileNotFound is thrown if the given file could not be found
    */
    template <typename MapType>
    void load(const std::string& filename, MapType& exp)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      exp.reset();

      std::ifstream is(filename.c_str());
      // get size of file
      is.seekg(0, std::ios::end);
      startProgress(0, is.tellg(), "loading MGF");
      is.seekg(0, std::ios::beg);

      UInt spectrum_number(0);
      Size line_number(0); // carry line number for error messages within getNextSpectrum()

      typename MapType::SpectrumType spectrum;
      spectrum.setMSLevel(2);
      spectrum.getPrecursors().resize(1);
      spectrum.setType(SpectrumSettings::SpectrumType::CENTROID); // MGF is always centroided, by definition
      while (getNextSpectrum_(is, spectrum, line_number, spectrum_number))
      {
        exp.addSpectrum(spectrum);
        setProgress(is.tellg());
        ++spectrum_number;
      } // next spectrum
      exp.updateRanges();
      endProgress();
    }

    /**
      @brief enclosing Strings of the peak list body for HTTP submission

      Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
      adding the peaks (in whatever format, e.g. mzXML) enclosed in this body.
      The @p filename can later be found in the Mascot response.
    */
    std::pair<std::string, std::string> getHTTPPeakListEnclosure(const std::string& filename) const;

    /// writes a spectrum in MGF format to an ostream
    void writeSpectrum(std::ostream& os, const PeakSpectrum& spec, const std::string& filename, const std::string& native_id_type_accession);

protected:

    /// use a compact format for storing (no zero-intensity peaks, limited number of decimal places)?
    bool store_compact_;

    /// mapping of modifications with specificity groups, that have to be treated specially (e.g. "Deamidated (NQ)")
    std::map<std::string, std::string> mod_group_map_;

    /// writes a parameter header
    void writeParameterHeader_(const std::string& name, std::ostream& os);

    /// write a list of (fixed or variable) modifications
    void writeModifications_(const std::vector<std::string>& mods, std::ostream& os,
                             bool variable_mods = false);

     /// writes the full header
    void writeHeader_(std::ostream& os);

    /// writes the MSExperiment
    void writeMSExperiment_(std::ostream& os, const std::string& filename, const PeakMap& experiment);

    /// reads a spectrum block, the section between 'BEGIN IONS' and 'END IONS' of a MGF file
    template <typename SpectrumType>
    bool getNextSpectrum_(std::ifstream& is, SpectrumType& spectrum, Size& line_number, const Size& spectrum_number)
    {
      spectrum.resize(0);
      spectrum.setNativeID(std::string("index=") + (spectrum_number));

      if (spectrum.metaValueExists("TITLE"))
      {
        spectrum.removeMetaValue("TITLE");
      }
      if (spectrum.metaValueExists("SEQ"))
      {
        // SEQ is a per-query field; do not let it bleed across spectra
        spectrum.removeMetaValue("SEQ");
      }
      typename SpectrumType::PeakType p;

      std::string line;
      // seek to next peak list block
      while (getline(is, line, '\n'))
      {
        ++line_number;

        StringUtils::trim(line); // remove whitespaces, line-endings etc

        // found peak list block?
        if (line == "BEGIN IONS")
        {
          while (getline(is, line, '\n'))
          {
            ++line_number;
            StringUtils::trim(line); // remove whitespaces, line-endings etc

            if (line.empty()) continue;

            if (isdigit(line[0])) // actual data .. this comes first, since its the most common case
            {
              std::vector<std::string> split;
              do
              {
                if (line.empty())
                {
                  continue;
                }

                StringUtils::simplify(line); // merge double spaces (explicitly allowed by MGF), to prevent empty split() chunks and subsequent parse error
                StringUtils::substitute(line, '\t', ' '); // also accept Tab (strictly, only space(s) are allowed)
                if (StringUtils::split(line, ' ', split, false))
                {
                  try
                  {
                    p.setPosition(StringUtils::toDouble(split[0]));
                    p.setIntensity(StringUtils::toDouble(split[1]));
                  }
                  catch (Exception::ConversionError& /*e*/)
                  {
                    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The content '" + line + "' at line #" + StringUtils::toStr(line_number) + " could not be converted to a number! Expected two (m/z int) or three (m/z int charge) numbers separated by whitespace (space or tab).", "");
                  }
                  spectrum.push_back(p);
                }
                else
                {
                  throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The content '" + line + "' at line #" + StringUtils::toStr(line_number) + " does not contain m/z and intensity values separated by whitespace (space or tab)!", "");
                }
              }
              while (getline(is, line, '\n') && ++line_number && StringUtils::trim(line) != "END IONS"); // StringUtils::trim(line) is important here!

              if (line == "END IONS")
              {
                return true; // found end of spectrum
              }
              else
              {
                throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, R"(Reached end of file. Found "BEGIN IONS" but not the corresponding "END IONS"!)", "");
              }
            }
            else if (StringUtils::hasPrefix(line, "PEPMASS")) // parse precursor position
            {
              std::string tmp = StringUtils::substr(line, 8); // copy since we might need the original line for error reporting later
              StringUtils::substitute(tmp, '\t', ' ');
              std::vector<std::string> split;
              StringUtils::split(tmp, ' ', split);
              if (split.size() == 1)
              {
                spectrum.getPrecursors()[0].setMZ(StringUtils::toDouble(StringUtils::trim(split[0])));
              }
              else if (split.size() == 2)
              {
                spectrum.getPrecursors()[0].setMZ(StringUtils::toDouble(StringUtils::trim(split[0])));
                spectrum.getPrecursors()[0].setIntensity(StringUtils::toDouble(StringUtils::trim(split[1])));
              }
              else
              {
                throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot parse PEPMASS in '" + line + "' at line #" + StringUtils::toStr(line_number) + " (expected 1 or 2 entries, but " + StringUtils::toStr(split.size()) + " were present)!", "");
              }
            }
            else if (StringUtils::hasPrefix(line, "CHARGE"))
            {
              std::string tmp = StringUtils::substr(line, 7);
              StringUtils::remove(tmp, '+');
              spectrum.getPrecursors()[0].setCharge(StringUtils::toInt32(tmp));
            }
            else if (StringUtils::hasPrefix(line, "RTINSECONDS"))
            {
              std::string tmp = StringUtils::substr(line, 12);
              spectrum.setRT(StringUtils::toDouble(tmp));
            }
            else if (StringUtils::hasPrefix(line, "TITLE"))
            {
              // test if we have a line like "TITLE= Cmpd 1, +MSn(595.3), 10.9 min"
              if (StringUtils::hasSubstring(line, "min"))
              {
                try
                {
                  std::vector<std::string> split;
                  StringUtils::split(line, ',', split);
                  if (!split.empty())
                  {
                    for (Size i = 0; i != split.size(); ++i)
                    {
                      if (StringUtils::hasSubstring(split[i], "min"))
                      {
                        std::vector<std::string> split2;
                        StringUtils::trim(split[i]);
                        StringUtils::split(split[i], ' ', split2);
                        if (!split2.empty())
                        {
                          StringUtils::trim(split2[0]);
                          spectrum.setRT(StringUtils::toDouble(split2[0]) * 60.0);
                        }
                      }
                    }
                  }
                }
                catch (Exception::BaseException& /*e*/)
                {
                  // just do nothing and write the whole title to spec
                  std::vector<std::string> split;
                  if (StringUtils::split(line, '=', split))
                  {
                    if (!split[1].empty()) spectrum.setMetaValue("TITLE", split[1]);
                  }
                }
              }
              else // just write the title as metainfo to the spectrum and add native ID to make the titles unique
              {
                Size firstEqual = line.find('=', 4);
                if (firstEqual != std::string::npos)
                {
                  if (StringUtils::hasSubstring(StringUtils::toStr(spectrum.getMetaValue("TITLE")), spectrum.getNativeID()))
                  {
                    spectrum.setMetaValue("TITLE", StringUtils::substr(line, firstEqual + 1));
                  }
                  else
                  {
                    spectrum.setMetaValue("TITLE", StringUtils::substr(line, firstEqual + 1) + "_" + spectrum.getNativeID());
                  }
                }
              }
            }
            else if (StringUtils::hasPrefix(line, "NAME"))
            {
              std::string tmp = StringUtils::substr(line, 5);
              spectrum.setMetaValue(Constants::UserParam::MSM_METABOLITE_NAME, tmp);
            }
            else if (StringUtils::hasPrefix(line, "COMPOUND_NAME"))
            {
              std::string tmp = StringUtils::substr(line, 14);
              spectrum.setMetaValue(Constants::UserParam::MSM_METABOLITE_NAME, tmp);
            }
            else if (StringUtils::hasPrefix(line, "INCHI="))
            {
              std::string tmp = StringUtils::substr(line, 6);
              spectrum.setMetaValue(Constants::UserParam::MSM_INCHI_STRING, tmp);
            }
            else if (StringUtils::hasPrefix(line, "SMILES"))
            {
              std::string tmp = StringUtils::substr(line, 7);
              spectrum.setMetaValue(Constants::UserParam::MSM_SMILES_STRING, tmp);
            }
            else if (StringUtils::hasPrefix(line, "IONMODE"))
            {
              std::string tmp = StringUtils::substr(line, 8);
              spectrum.setMetaValue("IONMODE", tmp);
            }
            else if (StringUtils::hasPrefix(line, "MSLEVEL")) 
            {
             std::string tmp = StringUtils::substr(line, 8);
            try 
            {
            int ms_level = std::stoi(tmp); 
            spectrum.setMSLevel(ms_level);
            }
            catch (const std::invalid_argument& /*e*/)
            {
            // Default to MS2 if parsing fails
            spectrum.setMSLevel(2);
            spectrum.setMetaValue("MSLEVEL", "2");
            }
            catch (const std::out_of_range& /*e*/)
            {
                spectrum.setMSLevel(2);
            }
            }               
            else if (StringUtils::hasPrefix(line, "SOURCE_INSTRUMENT"))
            {
              std::string tmp = StringUtils::substr(line, 18);
              spectrum.setMetaValue("SOURCE_INSTRUMENT", tmp);
            }
            else if (StringUtils::hasPrefix(line, "ORGANISM"))
            {
              std::string tmp = StringUtils::substr(line, 9);
              spectrum.setMetaValue("ORGANISM", tmp);
            }
            else if (StringUtils::hasPrefix(line, "PI"))
            {
              std::string tmp = StringUtils::substr(line, 3);
              spectrum.setMetaValue("PI", tmp);
            }
            else if (StringUtils::hasPrefix(line, "DATACOLLECTOR"))
            {
              std::string tmp = StringUtils::substr(line, 14);
              spectrum.setMetaValue("DATACOLLECTOR", tmp);
            }
            else if (StringUtils::hasPrefix(line, "LIBRARYQUALITY"))
            {
              std::string tmp = StringUtils::substr(line, 15);
              spectrum.setMetaValue("LIBRARYQUALITY", tmp);
            }
            else if (StringUtils::hasPrefix(line, "SPECTRUMID"))
            {
              std::string tmp = StringUtils::substr(line, 11);
              spectrum.setMetaValue("GNPS_Spectrum_ID", tmp);
            }
            else if (StringUtils::hasPrefix(line, "SCANS="))
            {
              std::string tmp = StringUtils::substr(line, 6);
              spectrum.setMetaValue("Scan_ID", tmp);
            }
            else if (StringUtils::hasPrefix(line, "SEQ="))
            {
              // Mascot sequence-query field: peptide sequence in one-letter code.
              // Per spec, SEQ may appear multiple times per query (each entry is
              // an independent sequence filter). Always stored as a StringList
              // under the "SEQ" key for a stable interface regardless of how
              // many SEQ lines were present.
              std::string sequence = StringUtils::substr(line, 4);
              StringList sequences;
              if (spectrum.metaValueExists("SEQ"))
              {
                sequences = spectrum.getMetaValue("SEQ").toStringList();
              }
              sequences.push_back(sequence);
              spectrum.setMetaValue("SEQ", sequences);
            }
          }
        }
      }

      return false; // found end of file
    }

  };
} // namespace OpenMS
