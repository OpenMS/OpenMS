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
#include <OpenMS/DATASTRUCTURES/String.h>
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
    void store(const String& filename, const PeakMap& experiment,
               bool compact = false);

    /// store the experiment data in a MascotGenericFile; the output is written to the given stream, the filename will be noted in the file (optionally a compact format is used: no zero-intensity peaks, limited number of decimal places)
    void store(std::ostream& os, const String& filename,
               const PeakMap& experiment, bool compact = false);

    /**
      @brief loads a Mascot Generic File into a PeakMap

      @param filename file name which the map should be read from
      @param exp the map which is filled with the data from the given file
      @throw FileNotFound is thrown if the given file could not be found
    */
    template <typename MapType>
    void load(const String& filename, MapType& exp)
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

      endProgress();
    }

    /**
      @brief enclosing Strings of the peak list body for HTTP submission

      Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
      adding the peaks (in whatever format, e.g. mzXML) enclosed in this body.
      The @p filename can later be found in the Mascot response.
    */
    std::pair<String, String> getHTTPPeakListEnclosure(const String& filename) const;

    /// writes a spectrum in MGF format to an ostream
    void writeSpectrum(std::ostream& os, const PeakSpectrum& spec, const String& filename, const String& native_id_type_accession);

protected:

    /// use a compact format for storing (no zero-intensity peaks, limited number of decimal places)?
    bool store_compact_;

    /// mapping of modifications with specificity groups, that have to be treated specially (e.g. "Deamidated (NQ)")
    std::map<String, String> mod_group_map_;

    /// writes a parameter header
    void writeParameterHeader_(const String& name, std::ostream& os);

    /// write a list of (fixed or variable) modifications
    void writeModifications_(const std::vector<String>& mods, std::ostream& os,
                             bool variable_mods = false);

     /// writes the full header
    void writeHeader_(std::ostream& os);

    /// writes the MSExperiment
    void writeMSExperiment_(std::ostream& os, const String& filename, const PeakMap& experiment);

    /// reads a spectrum block, the section between 'BEGIN IONS' and 'END IONS' of a MGF file
    template <typename SpectrumType>
    bool getNextSpectrum_(std::ifstream& is, SpectrumType& spectrum, Size& line_number, const Size& spectrum_number)
    {
      spectrum.resize(0);
      spectrum.setNativeID(String("index=") + (spectrum_number));

      if (spectrum.metaValueExists("TITLE"))
      {
        spectrum.removeMetaValue("TITLE");
      }
      typename SpectrumType::PeakType p;

      String line;
      // seek to next peak list block
      while (getline(is, line, '\n'))
      {
        ++line_number;

        line.trim(); // remove whitespaces, line-endings etc

        // found peak list block?
        if (line == "BEGIN IONS")
        {
          while (getline(is, line, '\n'))
          {
            ++line_number;
            line.trim(); // remove whitespaces, line-endings etc

            if (line.empty()) continue;

            if (isdigit(line[0])) // actual data .. this comes first, since its the most common case
            {
              std::vector<String> split;
              do
              {
                if (line.empty())
                {
                  continue;
                }

                line.simplify(); // merge double spaces (explicitly allowed by MGF), to prevent empty split() chunks and subsequent parse error
                line.substitute('\t', ' '); // also accept Tab (strictly, only space(s) are allowed)
                if (line.split(' ', split, false))
                {
                  try
                  {
                    p.setPosition(split[0].toDouble());
                    p.setIntensity(split[1].toDouble());
                  }
                  catch (Exception::ConversionError& /*e*/)
                  {
                    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The content '" + line + "' at line #" + String(line_number) + " could not be converted to a number! Expected two (m/z int) or three (m/z int charge) numbers separated by whitespace (space or tab).", "");
                  }
                  spectrum.push_back(p);
                }
                else
                {
                  throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The content '" + line + "' at line #" + String(line_number) + " does not contain m/z and intensity values separated by whitespace (space or tab)!", "");
                }
              }
              while (getline(is, line, '\n') && ++line_number && line.trim() != "END IONS"); // line.trim() is important here!

              if (line == "END IONS")
              {
                return true; // found end of spectrum
              }
              else
              {
                throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, R"(Reached end of file. Found "BEGIN IONS" but not the corresponding "END IONS"!)", "");
              }
            }
            else if (line.hasPrefix("PEPMASS")) // parse precursor position
            {
              String tmp = line.substr(8); // copy since we might need the original line for error reporting later
              tmp.substitute('\t', ' ');
              std::vector<String> split;
              tmp.split(' ', split);
              if (split.size() == 1)
              {
                spectrum.getPrecursors()[0].setMZ(split[0].trim().toDouble());
              }
              else if (split.size() == 2)
              {
                spectrum.getPrecursors()[0].setMZ(split[0].trim().toDouble());
                spectrum.getPrecursors()[0].setIntensity(split[1].trim().toDouble());
              }
              else
              {
                throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot parse PEPMASS in '" + line + "' at line #" + String(line_number) + " (expected 1 or 2 entries, but " + String(split.size()) + " were present)!", "");
              }
            }
            else if (line.hasPrefix("CHARGE"))
            {
              String tmp = line.substr(7);
              tmp.remove('+');
              spectrum.getPrecursors()[0].setCharge(tmp.toInt());
            }
            else if (line.hasPrefix("RTINSECONDS"))
            {
              String tmp = line.substr(12);
              spectrum.setRT(tmp.toDouble());
            }
            else if (line.hasPrefix("TITLE"))
            {
              // test if we have a line like "TITLE= Cmpd 1, +MSn(595.3), 10.9 min"
              if (line.hasSubstring("min"))
              {
                try
                {
                  std::vector<String> split;
                  line.split(',', split);
                  if (!split.empty())
                  {
                    for (Size i = 0; i != split.size(); ++i)
                    {
                      if (split[i].hasSubstring("min"))
                      {
                        std::vector<String> split2;
                        split[i].trim().split(' ', split2);
                        if (!split2.empty())
                        {
                          spectrum.setRT(split2[0].trim().toDouble() * 60.0);
                        }
                      }
                    }
                  }
                }
                catch (Exception::BaseException& /*e*/)
                {
                  // just do nothing and write the whole title to spec
                  std::vector<String> split;
                  if (line.split('=', split))
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
                  if (String(spectrum.getMetaValue("TITLE")).hasSubstring(spectrum.getNativeID()))
                  {
                    spectrum.setMetaValue("TITLE", line.substr(firstEqual + 1));
                  }
                  else
                  {
                    spectrum.setMetaValue("TITLE", line.substr(firstEqual + 1) + "_" + spectrum.getNativeID());
                  }
                }
              }
            }
            else if (line.hasPrefix("NAME"))
            {
              String tmp = line.substr(5);
              spectrum.setMetaValue(Constants::UserParam::MSM_METABOLITE_NAME, tmp);
            }
            else if (line.hasPrefix("INCHI="))
            {
              String tmp = line.substr(6);
              spectrum.setMetaValue(Constants::UserParam::MSM_INCHI_STRING, tmp);
            }
            else if (line.hasPrefix("SMILES"))
            {
              String tmp = line.substr(7);
              spectrum.setMetaValue(Constants::UserParam::MSM_SMILES_STRING, tmp);
            }
            else if (line.hasPrefix("SPECTRUMID"))
            {
              String tmp = line.substr(11);
              spectrum.setMetaValue("GNPS_Spectrum_ID", tmp);
            }
            else if (line.hasPrefix("SCANS="))
            {
              String tmp = line.substr(6);
              spectrum.setMetaValue("Scan_ID", tmp);
            }
          }
        }
      }

      return false; // found end of file
    }

  };

} // namespace OpenMS
