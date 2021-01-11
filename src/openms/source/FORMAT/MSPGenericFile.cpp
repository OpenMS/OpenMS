// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <boost/regex.hpp>
#include <fstream>

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
    synonyms_separator_ = (String)param_.getValue("synonyms_separator");
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
    boost::regex re_name("^Name: (.+)", boost::regex::no_mod_s);
    boost::regex re_synon("^synon(?:yms?)?: (.+)", boost::regex::no_mod_s | boost::regex::icase);
    boost::regex re_points_line("^\\d");
    boost::regex re_point("(\\d+(?:\\.\\d+)?)[: ](\\d+(?:\\.\\d+)?);? ?");
    boost::regex re_cas_nist("^CAS#: ([\\d-]+);  NIST#: (\\d+)"); // specific to NIST db
    boost::regex re_metadatum("^(.+): (.+)", boost::regex::no_mod_s);

    OPENMS_LOG_INFO << "\nLoading spectra from .msp file. Please wait." << std::endl;

    while (!ifs.eof())
    {
      ifs.getline(line, BUFSIZE);
      // Peaks
      if (boost::regex_search(line, m, re_points_line))
      {
        // OPENMS_LOG_DEBUG << "re_points_line\n";
        boost::regex_search(line, m, re_point);
        do
        {
          // OPENMS_LOG_DEBUG << "{" << m[1] << "} {" << m[2] << "}; ";
          const double position { std::stod(m[1]) };
          const double intensity { std::stod(m[2]) };
          spectrum.push_back( Peak1D(position, intensity) );
        } while ( boost::regex_search(m[0].second, m, re_point) );
      }
      // Synon
      else if (boost::regex_search(line, m, re_synon))
      {
        // OPENMS_LOG_DEBUG << "Synon: " << m[1] << "\n";
        synonyms_.push_back(String(m[1]));
      }
      // Name
      else if (boost::regex_search(line, m, re_name))
      {
        addSpectrumToLibrary(spectrum, library);
        // OPENMS_LOG_DEBUG << "\n\nName: " << m[1] << "\n";
        spectrum.clear(true);
        synonyms_.clear();
        spectrum.setName( String(m[1]) );
        spectrum.setMetaValue("is_valid", 1);
      }
      // Specific case of NIST's exported msp
      else if (boost::regex_search(line, m, re_cas_nist))
      {
        // OPENMS_LOG_DEBUG << "CAS#: " << m[1] << "; NIST#: " << m[2] << "\n";
        spectrum.setMetaValue(String("CAS#"), String(m[1]));
        spectrum.setMetaValue(String("NIST#"), String(m[2]));
      }
      // Other metadata
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

  void MSPGenericFile::addSpectrumToLibrary(
    MSSpectrum& spectrum,
    MSExperiment& library
  )
  {
    if (static_cast<int>(spectrum.getMetaValue("is_valid")) == 0) return;

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

      if (synonyms_.size())
      {
        String synon;
        for (const String& s : synonyms_)
        {
          synon += s + synonyms_separator_;
        }
        if (synon.size())
        {
          synon.pop_back();
        }
        spectrum.setMetaValue("Synon", synon);
      }

      spectrum.removeMetaValue("is_valid");
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
