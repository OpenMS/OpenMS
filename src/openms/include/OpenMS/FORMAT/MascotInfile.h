// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <fstream>

namespace OpenMS
{
  /**
    @brief Mascot input file adapter.

    @deprecated Use MascotGenericFile.h instead, which uses DefaultParamHandler, is more up to date and avoids weird behaviors of this class (especially in terms of
              MIME boundaries when writing MGF. This class is currently only used by MascotAdapter.cpp.

    Creates a file that can be used for Mascot search from a peak list or a whole experiment.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MascotInfile :
    public ProgressLogger
  {
public:

    /// constructor
    MascotInfile();

    /// constructor
    ~MascotInfile() override;

    /// stores the peak list in a MascotInfile that can be used as input for MASCOT shell execution
    void store(const String & filename, const PeakSpectrum & spec, double mz, double retention_time, String search_title);

    /// stores the experiment data in a MascotInfile that can be used as input for MASCOT shell execution
    void store(const String & filename, const PeakMap & experiment, String search_title);

    /** loads a Mascot Generic File into a PeakMap

            @param filename file name which the map should be read from
            @param exp the map which is filled with the data from the given file
            @throw FileNotFound is thrown if the given file could not be found
    */
    template <typename MapType>
    void load(const String & filename, MapType & exp)
    {
      exp.reset();
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      std::ifstream is(filename.c_str());
      std::vector<std::pair<double, double> > spec;
      UInt charge(0);
      double pre_mz(0), pre_int(0), rt(-1);
      String title;
      while (getNextSpectrum_(is, spec, charge, pre_mz, pre_int, rt, title))
      {
        typename MapType::SpectrumType spectrum;
        for (std::vector<std::pair<double, double> >::const_iterator it = spec.begin(); it != spec.end(); ++it)
        {
          typename MapType::PeakType p;
          p.setPosition(it->first);
          p.setIntensity(it->second);
          spectrum.push_back(p);
        }
        spectrum.setMSLevel(2);
        spectrum.getPrecursors().resize(1);
        spectrum.getPrecursors()[0].setMZ(pre_mz);
        spectrum.getPrecursors()[0].setIntensity(pre_int);
        spectrum.getPrecursors()[0].setCharge(charge);
        spectrum.setRT(rt);
        if (!title.empty())
        {
          spectrum.setMetaValue("TITLE", title);
          title = "";
        }

        exp.addSpectrum(spectrum);

        // clean up
        spec.clear();
        charge = 0;
        pre_mz = 0;
        pre_int = 0;
      }
    }

    /// returns the boundary used for the MIME format
    const String & getBoundary();
    /// sets the boundary used for the MIME format.<br>By default a 22 character random string is used
    void setBoundary(const String & boundary);

    /// returns the DB to use
    const String & getDB();
    /// sets the DB to use (default: MSDB). See &lt;mascot path&gt;/config/mascot.dat in "Databases" section for possible settings
    void setDB(const String & db);

    /// returns the search type
    const String & getSearchType();
    /// sets the search type (default: MIS). So far only MIS is supported!<br>Valid types are "MIS" (MS/MS Ion Search), "PMF" (Peptide Mass Fingerprint) , "SQ" (Sequence Query)
    void setSearchType(const String & search_type);

    /// returns the number of hits to report back
    const String & getHits();
    /// sets the number of hits to report back (default: 20)
    void setHits(const String & hits);

    /// returns the enzyme used for cleavage
    const String & getCleavage();
    /// sets the enzyme used for cleavage (default: Trypsin). <BR>See &lt;mascot path&gt;/config/enzymes for possible settings.
    void setCleavage(const String & cleavage);

    /// returns the used mass type ("Monoisotopic" or "Average")
    const String & getMassType();
    /// sets the used mass type "Monoisotopic" or "Average" (default: Monoisotopic)
    void setMassType(const String & mass_type);

    /// returns a vector containing the fixed modifications (default: none)
    const std::vector<String> & getModifications();
    /// sets the fixed modifications (default: none). <BR>See &lt;mascot path&gt;/config/mod_file for possible settings.
    void setModifications(const std::vector<String> & mods);

    /// returns a vector containing the variable modifications (default: none)
    const std::vector<String> & getVariableModifications();
    /// sets the fixed modifications (default: none). <BR>See &lt;mascot path&gt;/config/mod_file for possible settings.
    void setVariableModifications(const std::vector<String> & mods);

    /// returns the instrument type
    const String & getInstrument();
    /// sets the instrument type (Default: Default). <BR>Possible instruments: ESI-QUAD-TOF, MALDI-TOF-PSD, ESI-TRAP, ESI-QUAD, ESI-FTICR, MALDI-TOF-TOF, ESI-4SECTOR, FTMS-ECD, MALDI-QUAD-TOF, MALDI-QIT-TOF
    void setInstrument(const String & instrument);

    /// returns the number of allowed missed cleavages
    UInt getMissedCleavages() const;
    /// sets the number of allowed missed cleavages (default: 1)
    void setMissedCleavages(UInt missed_cleavages);

    /// returns the precursor mass tolerance
    float getPrecursorMassTolerance() const;
    /// sets the precursor mass tolerance in Da (default: 2.0)
    void setPrecursorMassTolerance(float precursor_mass_tolerance);

    /// returns the peak mass tolerance in Da
    float getPeakMassTolerance() const;
    /// sets the peak mass tolerance in Da (default: 1.0)
    void setPeakMassTolerance(float ion_mass_tolerance);

    /// returns the taxonomy
    const String & getTaxonomy();
    /// sets the taxonomy (default: All entries). <BR>See &lt;mascot path&gt;/config/taxonomy for possible settings.
    void setTaxonomy(const String & taxonomy);

    /// returns the Mascot form version
    const String & getFormVersion();
    /// sets the Mascot form version (default: 1.01)
    void setFormVersion(const String & form_version);

    /// returns the charges
    const String & getCharges();
    /// sets the charges (default: 1+, 2+ and 3+)
    void setCharges(std::vector<Int> & charges);

protected:
    /// parent mass
    double mz_;

    /// charge states to use
    String charges_;

    /// the search title of the mascot search
    String search_title_;

    /// the DB to search in
    String db_;

    /// search type: MIS, SQ or PMF
    String search_type_;

    /// number of hits to report
    String hits_;

    /// Enzyme used for cleavage
    String cleavage_;

    /// Monoisotopic/average mass
    String mass_type_;

    /// fixed Modifications
    std::vector<String> mods_;

    /// variable Modifications
    std::vector<String> variable_mods_;

    /// the used instrument
    String instrument_;

    /// number of missed cleavages
    UInt missed_cleavages_;

    /// precursor mass tolerance in Da
    float precursor_mass_tolerance_;

    /// m/z tolerance of ions  in Da
    float ion_mass_tolerance_;

    /// taxonomy
    String taxonomy_;

    /// form version
    String form_version_;

    /// the boundary used for the MIME format
    String boundary_;

    /// the retention time
    double retention_time_;

    /// writes a parameter header
    void writeParameterHeader_(const String & name, FILE * fp, bool line_break = true);

    /// writes the full header
    void writeHeader_(FILE * fp);

    /// writes the spectrum
    void writeSpectrum_(FILE * fp,
                        const String & filename,
                        const PeakSpectrum & peaks);

    /// writes the MSExperiment
    void writeMSExperiment_(FILE * fp,
                            const String & filename,
                            const PeakMap & experiment);

    bool getNextSpectrum_(std::istream & is, std::vector<std::pair<double, double> > & spectrum, UInt & charge, double & precursor_mz, double & precursor_int, double & rt, String & title);
  };

} // namespace OpenMS

